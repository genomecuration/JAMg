#!/usr/bin/env perl

=pod

=head1 TODO

wishlist:
interesting: overlapping genes: group as genes before deleting and rename gene_id
Low: Would be nice to have a GenBank parser that doesn't use BioPerl also augusts' gff2gbSmallDNA.pl is just too ...

=head1 NAME

 prepare_golden_genes_for_predictors.pl

 Uses an existing exonerate output or just predicted proteins (e.g. from PASA or just a FASTA file) to prepare gene prediction inputs for Augustus, SNAP and geneid.
 Exonerate is run (enabled via AAT) if it is not provided.
 GTF file is produced to judge quality of annotation.
 Sorts out high quality alignments from those that don't meet the criteria.
 

=head1 USAGE

Mandatory:

  * -genome       :s    => Genome sequence FASTA file. Repeatmasked with Ns if softmasked genome is also provided. Otherwise not repeatmasked.

And one of:

1. If existing run with exonerate:

    -exonerate    :s    => Exonerate results file with special format (see perldoc)
    -cdna               => Exonerate was run with cDNA otherwise expect protein (no start/stop searched in query, just genome)

2. PASA: output of PASA/scripts/pasa_asmbls_to_training_set.dbi:

    -pasa_gff      :s => GFF output from pasa
    -pasa_peptides :s => Protein output from pasa
    -pasa_assembly :s => Contigs provided as input to pasa
    -pasa_cds      :s => CDS output from pasa
    -pasa_genome   :s => Genome alignment GFF

3. OR

    -peptides              :s => A FASTA file with proteins

4. OR

    -mrna                  :s => A FASTA file with cDNA sequence (not fully implemented yet)

Other options:

    -help                     => Show this menu
    -training        :i       => number of genes to go into a random training set (def. 66% of total or 5000, whichever is higher)
    -complete                 => only do full length (recommended)
    -softmasked      :s       => Genome that has been softmasked for repeats
    -nosingle                 => Don't allow single exon genes. Recommended, most are repeats even though some are legitimate genes (they will be detected downstream in JAMg)
    -flanks          :i       => Base pairs to add left and right when creating GenBank file for training
    -identical       :i       => Identity cutoff to mark an exonerate aln as golden, out of 100 (def 95)
    -similar         :i       => Similarity cutoff to mark an exonerate aln as golden, for exonerate out of 100 (def 98)
    -intron          :i       => Max size of intron (def 70000 bp)
    -threads|cpu     :i       => Number of threads (def 1)
    -stop_exonerate           => Stop after exonerate finishes.
    -stop_golden              => Stop after sorting out which alignments were very good. Use it to prevent checks for Augustus/snap/BioPerl etc
    -minorf          :i       => Minimum size of ORF (and therefore amino acids; def. 290 bp)
    -mismatch_cutoff :i       => Maximum number of mismatches allowed in exonerate (def 10)
    -same_species             => Contigs/proteins and genome is the same species
    -norefine                 => Don't use -refine for exonerate
    -norerun                  => Don't re-run exonerate (assume it already exists as [input file].exonerate.results
    -augustus                 => Directory where Augustus is installed (if not in your path)
    -extra_gff       :s       => Any extra GFF lines to consider (give a file)

=head1 DESCRIPTION

 Uses the following option in exonerate:

 --ryo "RYOAP_START\nRYOAP_STATS_D\tAlignment length\tidentical\tsimilar\tmismatch\tidentical%%\tsimilar%%\nRYOAP_STATS\t%s\t%et\t%ei\t%es\t%em\t%pi\t%ps\nRYOAP_CODING_QUERY\t%qi\t%qab\t%qae\n>%qi\n%qas\nRYOAP_CODING_GENOME\t%ti\t%tcb\t%tce\n>%ti\n%tcs\nRYOAP_END\n"

 To create softmasked FASTA, use RepeatMasker with the -gff option and then bedtools' maskFastaFromBed program.

=head1 Criteria for pass

 Methionine and * in protein sequence
 Start codon and stop codon in genome
 Fraction of identical residues > 98%
 Fraction of similar residues > 95%
 No more than 10 mismatches allowed
 Overlaps: For a particular genomic start/stop co-ordinate & length, only one gene is printed (one with best score); NO isoforms (alt. splicings) allowed.

=head1 AUTHORS

 Alexie Papanicolaou

        CSIRO Ecosystem Sciences
        alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

Copyright 2012-2014 the Commonwealth Scientific and Industrial Research Organization. 
This software is released under the Mozilla Public License v.2.

It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.mozilla.org/MPL/2.0.


=head1 BUGS & LIMITATIONS

Please report.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use POSIX qw(ceil);
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} =
"$RealBin:$RealBin/../3rd_party/bin/:$RealBin:$RealBin/../3rd_party/RepeatMasker/ncbi-blast/bin/"
  . $ENV{PATH};
use threads;
use Thread_helper;
use File::Basename;
use List::Util 'shuffle';

use Gene_obj;
use Gene_obj_indexer;
use CdbTools;
use GFF3_utils;
use GTF_utils;
use Nuc_translator;
use Fasta_reader;

# error handling
use Data::Dumper;

# debug
use Carp;

#options
my (
     $debug,              $contig_file,           $exonerate_file,
     $genome_file,        $is_cdna,               $no_single_exon,
     $overwrite,          $only_complete,         $peptide_file,
     $stop_after_golden,  $pasa_gff,              $pasa_cds,
     $pasa_assembly_file, $pasa_peptides,         $mrna_file,
     $softmasked_genome,  $stop_after_correction, $norefine,
     $nodataprint,        $no_gmap,               $no_exonerate,
     $pasa_genome_gff,    $extra_gff_file, $show_help
);

my $no_rerun_exonerate;
my $same_species             = '';
my $intron_size              = 70000;
my $training_set_size        = 4000;
my $aug_optimization_geneset = 450;
my $augustus_flank_region    = 4000;
my $identical_fraction       = 95;
my $minorf                   = 290;     #minimum orf size in bp
my $similar_fraction         = 98;
my $threads                  = 1;
my $mismatch_cutoff          = 10;

#globals
my $failed_cutoff = int(0);
my ( $cdbfasta_exec, $cdbyank_exec ) = &check_program( 'cdbfasta', 'cdbyank' );
my ( %get_id_seq_from_fasta_hash, $augustus_dir );

GetOptions(
            'help'               => \$show_help,
            'debug|verbose'      => \$debug,
            'exonerate:s'        => \$exonerate_file,
            'genome:s'           => \$genome_file,
            'softmasked:s'       => \$softmasked_genome,
            'pasa_gff:s'         => \$pasa_gff,
            'pasa_peptides:s'    => \$pasa_peptides,
            'pasa_cds:s'         => \$pasa_cds,
            'pasa_assembly:s'    => \$pasa_assembly_file,
            'pasa_genome:s'      => \$pasa_genome_gff,
            'peptides:s'         => \$peptide_file,
            'mrnas:s'            => \$mrna_file,
            'cdna'               => \$is_cdna,
            'training:i'         => \$training_set_size,
            'flanks:i'           => \$augustus_flank_region,
            'identical:i'        => \$identical_fraction,
            'similar:i'          => \$similar_fraction,
            'overwrite'          => \$overwrite,
            'intron:i'           => \$intron_size,
            'threads|cpu:i'      => \$threads,
            'nosingle|no_single' => \$no_single_exon,
            'complete'           => \$only_complete,
            'stop_exonerate'     => \$stop_after_correction,
            'stop_golden'        => \$stop_after_golden,
            'minorf:i'           => \$minorf,
            'mismatch_cutoff:i'  => \$mismatch_cutoff,
            'same_species'       => \$same_species,
            'norefine'           => \$norefine,
            'norerun'            => \$no_rerun_exonerate,
            'nodataprint'        => \$nodataprint,
            'augustus_dir:s'     => \$augustus_dir,
            'no_exonerate'       => \$no_exonerate,
            'no_gmap'            => \$no_gmap,
            'extra_gff:s'        => \$extra_gff_file
);

my ( $makeblastdb_exec, $tblastn_exec, $tblastx_exec ) =
  &check_program( 'makeblastdb', 'tblastn', 'tblastx' );
my ( $gmap_build_exec, $gmap_exec ) = &check_program( 'gmap_build', 'gmap' );

my ( $gff2gb_exec, $fathom_exec, $augustus_exec, $augustus_train_exec,
     $augustus_filterGenes_exec )
  = &check_program_optional(
                             'gff2gbSmallDNA.pl', 'fathom',
                             'augustus',          'etraining',
                             'filterGenes.pl'
  );
&check_for_options();

my $genome_sequence_file =
    $softmasked_genome
  ? $softmasked_genome
  : $genome_file;    # the full sequence, not repeatmasked
print "Indexing genome\n";
my $genome_sequence_file_dir  = dirname($genome_sequence_file);
my $genome_sequence_file_base = basename($genome_sequence_file);
my $genome_dir                = basename($genome_file) . '_dir';
my $scaffold_seq_hashref      = &read_fasta($genome_sequence_file);

&process_cmd(
"$makeblastdb_exec -in $genome_sequence_file -out $genome_sequence_file -hash_index -parse_seqids -dbtype nucl"
) unless -s "$genome_sequence_file.nin";
&process_cmd(
"$gmap_build_exec -D $genome_sequence_file_dir -d $genome_sequence_file_base.gmap $genome_file"
) unless -d "$genome_sequence_file_dir/$genome_sequence_file_base.gmap";

unlink("$genome_sequence_file.cidx");
&process_cmd("$cdbfasta_exec $genome_sequence_file");
my $uppercase_genome_files = &splitfasta( $genome_file, $genome_dir, 1 );

# We now have the option to use exonerate or another method.
# It turns out that exonerate is very good at getting rid of the false positives
# but aatpackage is not good enough as input when the species of the genome is the
# same as the other mRNA. Instead we can use BLAST (slower).
# also, in addition to exonerate, gmap might be of use. it will need post-processing
# but it might be better than exonerate (it is certainly MUCH faster).

my ( $gmap_gff, $gmap_passed, $exonerate_gff, $exonerate_passed,
     %passed_check, @to_evaluate );
unless ($no_gmap) {
 if ( $mrna_file && -s $mrna_file ) {
  ( $gmap_gff, $gmap_passed ) = &run_gmap($mrna_file);
 }
 elsif ( $pasa_cds && -s $pasa_cds ) {
  ( $gmap_gff, $gmap_passed ) = &run_gmap($pasa_cds);
 }
}
( $exonerate_gff, $exonerate_passed ) = &run_exonerate() unless $no_exonerate;

my $final_gff = 'final_golden_genes.gff3';
open( OUT1, ">$final_gff" );
open( OUT2, ">$final_gff.passed" );

# if hit in both, first is kept
# pasa/exonerate takes precedence to gmap
if ( $exonerate_gff && -s $exonerate_gff ) {
 open( IN, $exonerate_gff ) || die( "Cannot open $exonerate_gff " . $! );
 push(@to_evaluate,$exonerate_gff);
 while ( my $ln = <IN> ) { print OUT1 $ln; }
 close IN;
 open( IN, $exonerate_passed ) || die( "Cannot open $exonerate_passed " . $! );
 while ( my $ln = <IN> ) {
  if ( $ln =~ /^(\S+)/ ) {
   next if $passed_check{$1};
   $passed_check{$1}++;
   print OUT2 $ln;
  }
 }
 close IN;
}

if ( $gmap_gff && -s $gmap_gff ) {
 open( IN, $gmap_gff ) || die( "Cannot open $gmap_gff " . $! );
 push(@to_evaluate,$gmap_gff);
 while ( my $ln = <IN> ) {
  print OUT1 $ln;
 }
 close IN;
 open( IN, $gmap_passed ) || die( "Cannot open $gmap_passed " . $! );
 while ( my $ln = <IN> ) {
  if ( $ln =~ /^(\S+)/ ) {
   next if $passed_check{$1};
   $passed_check{$1}++;
   print OUT2 $ln;
  }
 }
 close IN;
}

if ( $extra_gff_file && -s $extra_gff_file ) {
 print "Processing user-provided GFF3 $extra_gff_file\n";
 &create_golden_gffs_from_gff( $extra_gff_file, $extra_gff_file . '.n' );
 &remove_overlapping_gff( $extra_gff_file . '.n', $extra_gff_file . '.nr' );
 $extra_gff_file .= '.nr';
 &create_golden_gffs_from_gff($extra_gff_file);
 push(@to_evaluate,$extra_gff_file.'.golden');
 open( IN, $extra_gff_file.'.golden' ) || die( "Cannot open $extra_gff_file.golden " . $! );
 while ( my $ln = <IN> ) { print OUT1 $ln; }
 open( IN, $extra_gff_file . 'passed' )
   || die( "Cannot open $extra_gff_file.passed " . $! );

 while ( my $ln = <IN> ) {
  if ( $ln =~ /^(\S+)/ ) {
   next if $passed_check{$1};
   $passed_check{$1}++;
   print OUT2 $ln;
  }
 }
 close IN;
}

close OUT1;
close OUT2;

if ($stop_after_golden) {
 print "User stop requested\n";
 exit(0);
}

if ( $final_gff && -s $final_gff ) {
 &remove_overlapping_gff( $final_gff, "$final_gff.nr",
                          undef,      "$final_gff.passed" );
 &process_for_gene_prediction( $final_gff . '.nr', "$final_gff.passed" );
 print "Done!\n";
}
else {
 warn "There was an issue producing the final GFF file. Exiting\n";
}
#################################################################################################################
######################################### Subroutines ###########################################################
#################################################################################################################
sub correct_exonerate_gff() {
 my $exonerate_file = shift;
 print
"Processing $exonerate_file using these cut-offs: identical:$identical_fraction similar:$similar_fraction;  No more than $mismatch_cutoff mismatches in total; Methionine and * in protein sequence; Start codon and stop codon in genome; No overlapping genes; Splice sites (GT..AG or GC..AG)\n";

 #pod2usage unless $exonerate_file && -s $exonerate_file;
 # IF WE NEED TO GET WHOLE SCAFFOLD SEQUENCE:
 pod2usage "Cannot find a file $exonerate_file or $genome_sequence_file\n"
   unless $exonerate_file
    && $genome_sequence_file
    && -s $exonerate_file
    && -s $genome_sequence_file;
 open( EXONERATE,               $exonerate_file );
 open( CORRECTED_EXONERATE_GFF, ">$exonerate_file.corrected.gff3" );

 my ( %details, $flag, %vulgar_data, %already_printed );
 my ( $gene_counter, $scounter, $ecounter ) = ( int(0), int(0), int(0) );

###########################################################################################
# EXONERATE protein has no UTR, the gene is actually the coding part of the mRNA
# get gene data
###########################################################################################
 my $header = <EXONERATE>;
 die
"It seems that the exonerate was run using protein2genome but -cdna was requested. Aborting\n"
   if ( $header =~ /protein2genome/ && $is_cdna );
 die
"It seems that the exonerate was run using cdna2genome but -cdna was not requested. Aborting\n"
   if ( ( $header =~ /cdna2genome/ || $header =~ /coding2genome/ )
        && !$is_cdna );
 while ( my $ln = <EXONERATE> ) {
  next if $ln =~ /^\s*$/;
  if ( $ln =~ /^vulgar:\s+(\S+)/ ) {
   $vulgar_data{$1} = $ln;
  }
  if ( $ln =~ /^# --- START OF GFF DUMP/ ) {
   my @exonerate_data;
   my @overlap_check;
   while ( my $ln2 = <EXONERATE> ) {
    last if ( $ln2 =~ /^# --- END OF GFF DUMP/ );
    next if ( $ln2 =~ /^#/ );
    my @data = split( "\t", $ln2 );
    if ( $data[8] ) {
     next if $data[2] eq 'similarity';
     my $offset = int(0);

# this is where the offset is taken from. hope that people don't run with references with this coords...
     if ( $data[0] =~ s/:(\d+)-(\d+)$// ) {
      $offset = $1 - 1;
      $data[3] += $offset;
      $data[4] += $offset;
     }
     push( @exonerate_data, join( "\t", @data ) );
     if ( $data[2] eq 'cds' ) {
      my %h = ( 'start' => $data[3], 'end' => $data[4] );
      push( @overlap_check, \%h );
     }
    }
   }
###########################################################################################
   # making other corrections. and fixing exonerate, snap etc bugs
###########################################################################################
   my ( $gene_start, $gene_end, $mRNA, $mRNA_id );
   my $exon_counter   = int(0);
   my $cds_counter    = int(0);
   my $splice_counter = int(0);
   $exonerate_data[0] =~ /sequence (\S+)/;
   my $gene_id = $1 if $1;

# if it has already been printed, skip it. someitmes there is a bug and this creates non-unique IDs
   if ( $already_printed{$gene_id} ) {
    warn "Warning: $gene_id was already printed before. skipping...\n";
    next;
   }
   $already_printed{$gene_id} = 1;

   for ( my $i = 0 ; $i < scalar(@exonerate_data) ; $i++ ) {
    my $gff_line = $exonerate_data[$i] || next;
    my @data = split( "\t", $gff_line );
    if ( $data[2] eq 'gene' ) {
     $data[8] =~ /sequence\s+(\S+)/;
     $gene_id    = $1;
     $gene_start = $data[3];
     $gene_end   = $data[4];
     undef($mRNA);    #reset
     undef($mRNA_id);

     #PASA HACK - review if there are problems in the future
     if ( $gene_id =~ s/^m\.(\d+)$/g.$1/ ) {
      $mRNA_id = "m.$1";
      $data[8] = "ID=$gene_id\n";
      $mRNA =
          $data[0] . "\t"
        . $data[1]
        . "\tmRNA\t"
        . $data[3] . "\t"
        . $data[4] . "\t.\t"
        . $data[6] . "\t.\t"
        . "ID=$mRNA_id;Parent=$gene_id\n";
     }

     # not PASA
     else {
      $mRNA_id = $gene_id . ".mRNA";
      $mRNA =
          $data[0] . "\t"
        . $data[1]
        . "\tmRNA\t"
        . $data[3] . "\t"
        . $data[4] . "\t.\t"
        . $data[6] . "\t.\t"
        . "ID=$mRNA_id;Parent=$gene_id\n";

      $data[8] = "ID=$gene_id\n";
     }
     $details{$gene_id}{'source'} = $data[1];
    }
    elsif ( $data[2] =~ /utr/ ) {
     if ( $data[2] =~ /5/ ) {
      $data[2] = 'five_prime_utr';
     }
     elsif ( $data[2] =~ /3/ ) {
      $data[2] = 'three_prime_utr';
     }
     else {
      $data[2] = 'UTR';
     }

     $data[8] = "ID=$gene_id." . $data[2] . ";Parent=$mRNA_id\n";
    }
    elsif ( $data[2] =~ /splice/ ) {
     $splice_counter++;
     if ( $data[8] =~ /splice_site "([A-Z]+)"/ ) {
      my $site = $1;
      if ( $data[2] eq 'splice5' ) {
       unless ( $site eq 'GT' || $site eq 'GC' ) {
        $details{$gene_id}{'non_canonical_splice_sites'}++;
       }
      }
      elsif ( $data[2] eq 'splice3' && $site ne 'AG' ) {
       $details{$gene_id}{'non_canonical_splice_sites'}++;
      }
     }
     $data[8] =
       "ID=$gene_id." . $data[2] . ":$splice_counter;Parent=$mRNA_id\n";
     $data[2] = 'splice_junction';

    }
    elsif ( $data[2] eq 'cds' ) {
     if ( $cds_counter > 0 ) {
      ### check if the next line is an exon. if it is not, check if there is an overlapping previous CDS (exonerate cdna2genome bug)
      # if there is we have to delete $exonerate_data[$i]; and move on
      my $previous = $overlap_check[ $cds_counter - 1 ];
      my @next_data = split( "\t", $exonerate_data[ $i + 1 ] )
        if $exonerate_data[ $i + 1 ];
      if (    $previous
           && $next_data[2]
           && $next_data[2] ne 'exon' )
      {
       if ( $data[6] eq '+' ) {

        # start must be higher than previous end
        if ( $data[3] < $previous->{'end'} ) {

#warn $data[8].":Start is smaller than previous end:\n".$data[3]." vs ".$previous->{'end'}."\n";
         delete( $exonerate_data[$i] );
         next;
        }
       }
       else {

        # end must be lower than previous start
        if ( $data[4] > $previous->{'start'} ) {
         delete( $exonerate_data[$i] );
         next;
        }
       }
      }
     }
     elsif ( scalar(@overlap_check) == 1 ) {

# if this is the first and only CDS, then it is a single exon gene. exonerate cdna2genome bug includes UTR co-ordidates. fix it
# by getting previous UTR line and adjust co-ordidates to by utr[end]+1 if + or utr[start]-1 if -.
      my @previous_data =
        split( "\t", $exonerate_data[ $i - 1 ] );
      if ( @previous_data && $previous_data[2] =~ /utr/i ) {
       if ( $data[6] eq '+' ) {
        if ( $data[3] < $previous_data[4] ) {
         $data[3] = $previous_data[4] + 1;
        }
       }
       else {
        if ( $data[4] > $previous_data[3] ) {
         $data[4] = $previous_data[3] - 1;
        }
       }
      }
     }
     $cds_counter++;
     $data[2] = 'CDS';
     $data[8] = "ID=cds.$gene_id;Parent=$mRNA_id\n";
    }
    elsif ( $data[2] eq 'exon' ) {
     $exon_counter++;
     $data[8] = "ID=$gene_id.exon.$exon_counter;Parent=$mRNA_id\n";
    }
    elsif ( $data[2] eq 'intron' ) {

     # we capture this from the GFF now;
     next;

     #$intron_counter++;
     #$data[8] = "ID=$gene_id.intron.$intron_counter;Parent=$mRNA_id\n";
    }
    $data[8] =~ s/\s+\;\s+/\;/g;
    $exonerate_data[$i] = join( "\t", @data );
   }
#################################################################################
   # now printing data
################################################################################
   for ( my $i = 0 ; $i < scalar(@exonerate_data) ; $i++ ) {
    my $gff_line = $exonerate_data[$i] || next;
    my @data = split( "\t", $gff_line );
    if ( $data[2] eq 'gene' ) {
     print CORRECTED_EXONERATE_GFF join( "\t", @data );
     print CORRECTED_EXONERATE_GFF $mRNA if $mRNA;
     $gene_counter++;
    }
    else {
     print CORRECTED_EXONERATE_GFF join( "\t", @data );
    }
   }

   # delimit end of gene model
   print CORRECTED_EXONERATE_GFF "##\n";

  }
###########################################################################################
  # for golden set:
  if ( $ln =~ /^RYOAP_END/ ) {
   $flag = 0;
   $ecounter++;
  }
  if ( $ln =~ /^RYOAP_START/ ) {
   $flag = 1;
   $scounter++;
   my $stats_header = <EXONERATE>;
   my $stats_str    = <EXONERATE>;
   chomp($stats_str);
   my @stats_data = split( "\t", $stats_str );
   my $query_data_str = <EXONERATE>;
   chomp($query_data_str);
   my @query_data = split( "\t", $query_data_str );
   my $query_id = $query_data[1];
   my $query_start =
       $query_data[2] <= $query_data[3]
     ? $query_data[2]
     : $query_data[3];
   my $query_end =
       $query_data[3] >= $query_data[2]
     ? $query_data[3]
     : $query_data[2];
   my $query_fasta = <EXONERATE>;
   my $query_seq;
   $ln = <EXONERATE>;

   while ( $ln !~ /^RYOAP/ ) {
    chomp($ln);
    $query_seq .= $ln;
    $ln = <EXONERATE>;
   }

   #       $query_fasta .= $query_seq;
   #TARGET
   my $target_data_str = $ln;
   chomp($target_data_str);
   my @target_data = split( "\t", $target_data_str );
   my $target_id = $target_data[1];
   my $target_start =
       $target_data[2] <= $target_data[3]
     ? $target_data[2]
     : $target_data[3];
   my $target_end =
       $target_data[3] >= $target_data[2]
     ? $target_data[3]
     : $target_data[2];
   my $alignment_strand = $target_data[2] <= $target_data[3] ? int(1) : int(-1);
   my $target_fasta = <EXONERATE>;
   my $target_seq;
   $ln = <EXONERATE>;

   while ( $ln !~ /^RYOAP/ ) {
    chomp($ln);
    $target_seq .= $ln;
    $ln = <EXONERATE>;
   }
   if ( $ln =~ /^RYOAP_END/ ) {
    $flag = 0;
    $ecounter++;
   }

   #       $target_fasta .= $target_seq;
   my $qlength = $query_data[3] - $query_data[2];
   my $tlength = abs( $target_data[3] - $target_data[2] );
   my $offset  = int(0);
   if ( $target_id =~ s/:(\d+)\-\d+$// ) {
    $offset = $1 - 1;
   }
   $details{$query_id}{'counter'}                    = $scounter;
   $details{$query_id}{'qlength'}                    = $qlength;
   $details{$query_id}{'tlength'}                    = $tlength + $offset;
   $details{$query_id}{'score'}                      = $stats_data[1];
   $details{$query_id}{'mismatch_number'}            = $stats_data[5];
   $details{$query_id}{'identical_frac'}             = $stats_data[6];
   $details{$query_id}{'similar_frac'}               = $stats_data[7];
   $details{$query_id}{'query_id'}                   = $query_id;
   $details{$query_id}{'query_seq'}                  = uc($query_seq);
   $details{$query_id}{'target_seq'}                 = uc($target_seq);
   $details{$query_id}{'query_start'}                = $query_start;
   $details{$query_id}{'query_end'}                  = $query_end;
   $details{$query_id}{'target_id'}                  = $target_id;
   $details{$query_id}{'target_start'}               = $target_start + $offset;
   $details{$query_id}{'target_end'}                 = $target_end + $offset;
   $details{$query_id}{'alignment_strand'}           = $alignment_strand;
   $details{$query_id}{'non_canonical_splice_sites'} = int(0)
     if !$details{$query_id}{'non_canonical_splice_sites'};
  }
 }
 print "Processed $gene_counter gene models.\n";
 close(EXONERATE);
 close CORRECTED_EXONERATE_GFF;
 die "Report seems incomplete!\n"
   unless $scounter == $ecounter;
 die "No genes have been found!\n" if $scounter == 0;

 if ($stop_after_correction) {
  print "User asked to stop after correction\n";
  exit();
 }
 &create_golden_gffs_from_exonerate( $exonerate_file . '.corrected.gff3',
                                     \%details, \%vulgar_data )
   unless -s "$exonerate_file.corrected.passed";
}

sub create_golden_gffs_from_exonerate() {
 my $exonerate_file      = shift;
 my $details_hashref     = shift;
 my $vulgar_data_hashref = shift;

 &gff3_fix_phase($exonerate_file) unless -s $exonerate_file . '.gtf';

 # it's getting a bit long...
 my $basename_exonerate = $exonerate_file;
 $basename_exonerate =~ s/\.gff3$//;
 print "Finding golden subset...\n";

 my $query_outfile =
   $is_cdna
   ? "$basename_exonerate.passed.cDNA"
   : "$basename_exonerate.passed.protein";
 open( OUT, ">$basename_exonerate.data" ) unless $nodataprint;
 open( VULGAR, ">$basename_exonerate.passed.vulgar" );
 open( LOG,    ">$basename_exonerate.golden.log" );
 open( PASSP,  ">$query_outfile" );

 #  open( PASSG,  ">$basename_exonerate.passed.genome" )  unless $nodataprint;
 open( PASSGC, ">$basename_exonerate.passed.genome.cds" ) unless $nodataprint;
 print OUT
"#Hit number\tQuery length\tTarget length\tScore\tNumber of mismatches\tFraction of identical\tFraction of similar\tNon-canonical splice_sites\tquery name\tquery start\tquery end\tquery sequence\tAlignment direction\ttarget name\ttarget start\ttarget end\ttarget sequence\n"
   unless $nodataprint;
 my $query_number = keys %{$details_hashref};
 my $pass_counter = int(0);
 my %pass_check;
 my %accepted;
 my $region_exists = int(0);
 my $counter       = int(0);

 foreach my $query_id (
  sort {
   $details_hashref->{$b}{'score'} <=> $details_hashref->{$a}{'score'}
  }
  keys %{$details_hashref}
   )
 {
  $counter++;
  print "\r$counter/$query_number    $pass_counter passed           ";
  print OUT "Hit "
    . $details_hashref->{$query_id}{'counter'} . "\t"
    . $details_hashref->{$query_id}{'qlength'} . "\t"
    . $details_hashref->{$query_id}{'tlength'} . "\t"
    . $details_hashref->{$query_id}{'score'} . "\t"
    . $details_hashref->{$query_id}{'mismatch_number'} . "\t"
    . $details_hashref->{$query_id}{'identical_frac'} . "\t"
    . $details_hashref->{$query_id}{'similar_frac'} . "\t"
    . $details_hashref->{$query_id}{'non_canonical_splice_sites'} . "\t"
    . $query_id . "\t"
    . $details_hashref->{$query_id}{'query_start'} . "\t"
    . $details_hashref->{$query_id}{'query_end'} . "\t"
    . $details_hashref->{$query_id}{'query_seq'} . "\t"
    . $details_hashref->{$query_id}{'alignment_strand'} . "\t"
    . $details_hashref->{$query_id}{'target_id'} . "\t"
    . $details_hashref->{$query_id}{'target_start'} . "\t"
    . $details_hashref->{$query_id}{'target_end'} . "\t"
    . $details_hashref->{$query_id}{'target_seq'} . "\n"
    unless $nodataprint;
  if (
       $pass_check{
            's'
          . $details_hashref->{$query_id}{'target_start'} . 'e'
          . $details_hashref->{$query_id}{'target_end'}
       }
    )
  {
   print LOG
     "$query_id rejected: higher scoring alignment with same co-ordinates ("
     . 's'
     . $details_hashref->{$query_id}{'target_start'} . 'e'
     . $details_hashref->{$query_id}{'target_end'} . "):"
     . $pass_check{ 's'
      . $details_hashref->{$query_id}{'target_start'} . 'e'
      . $details_hashref->{$query_id}{'target_end'} }
     . "\n";
   $region_exists++;
   next;
  }
  unless (
        $details_hashref->{$query_id}{'identical_frac'} >= $identical_fraction )
  {
   print LOG
     "$query_id rejected: identical_frac less than $identical_fraction: "
     . $details_hashref->{$query_id}{'identical_frac'} . "\n";
   $failed_cutoff++;
   next;
  }
  unless ( $details_hashref->{$query_id}{'similar_frac'} >= $similar_fraction )
  {
   print LOG "$query_id rejected: similar_frac less than $similar_fraction: "
     . $details_hashref->{$query_id}{'similar_frac'} . "\n";
   $failed_cutoff++;
   next;
  }
  unless (
          $details_hashref->{$query_id}{'mismatch_number'} <= $mismatch_cutoff )
  {
   print LOG "$query_id rejected: mismatches more than $mismatch_cutoff: "
     . $details_hashref->{$query_id}{'mismatch_number'} . "\n";
   $failed_cutoff++;
   next;
  }
  if ( $details_hashref->{$query_id}{'non_canonical_splice_sites'} > 0 ) {
   print LOG "$query_id rejected: there are non-canonical splice sites : "
     . $details_hashref->{$query_id}{'non_canonical_splice_sites'} . "\n";
   $failed_cutoff++;
   next;
  }

  unless ( $is_cdna
          || substr( $details_hashref->{$query_id}{'query_seq'}, 0, 1 ) eq 'M' )
  {
   print LOG "$query_id rejected: query does not start with M: "
     . substr( $details_hashref->{$query_id}{'query_seq'}, 0, 1 ) . "\n";
   $failed_cutoff++;
   next;
  }
  unless ( $is_cdna
         || substr( $details_hashref->{$query_id}{'query_seq'}, -1, 1 ) eq '*' )
  {
   print LOG "$query_id rejected: query does not end with stop codon (*): "
     . substr( $details_hashref->{$query_id}{'query_seq'}, -1, 1 ) . "\n";
   $failed_cutoff++;
   next;
  }
  unless (
          substr( $details_hashref->{$query_id}{'target_seq'}, 0, 3 ) eq 'ATG' )
  {
   print LOG "$query_id rejected: target does not start with ATG: "
     . substr( $details_hashref->{$query_id}{'target_seq'}, 0, 3 ) . "\n";
   $failed_cutoff++;
   next;
  }
  unless (
       (
           substr( $details_hashref->{$query_id}{'target_seq'}, -3, 3 ) eq 'TAG'
        || substr( $details_hashref->{$query_id}{'target_seq'}, -3, 3 ) eq 'TAA'
        || substr( $details_hashref->{$query_id}{'target_seq'}, -3, 3 ) eq 'TGA'
       )
    )
  {
   print LOG
"$query_id rejected: target does not end with stop codon (TAG, TAA or TGA): "
     . substr( $details_hashref->{$query_id}{'target_seq'}, -3, 3 ) . "\n";
   $failed_cutoff++;
   next;
  }
  my $scaffold_seq =
    $scaffold_seq_hashref->{ $details_hashref->{$query_id}{'target_id'} };
  if ( !$scaffold_seq ) {
   warn "Cannot find genome entry "
     . $details_hashref->{$query_id}{'target_id'} . "\n";
   next;
  }

  #bug check
  unless (int( $details_hashref->{$query_id}{'target_start'} )
       && int( $details_hashref->{$query_id}{'target_end'} )
       && $details_hashref->{$query_id}{'target_start'} > 0
       && $details_hashref->{$query_id}{'target_end'} > 0
       && $details_hashref->{$query_id}{'target_end'} <= length($scaffold_seq) )
  {
   print LOG "$query_id had to be skipped: has issues with target end ("
     . $details_hashref->{$query_id}{'target_id'}
     . ") being negative or larger than scaffold size - possible corruption of exonerate file! ("
     . $details_hashref->{$query_id}{'target_end'}
     . " is <0 or > "
     . length($scaffold_seq) . ")\n";
   next;
  }
  my $scaffold_subseq;

  #consider: adding stop codon as it is not reported at end??
  if ( $details_hashref->{$query_id}{'alignment_strand'} == -1 ) {
   $scaffold_subseq = &revcomp(
    substr(
            $scaffold_seq,
            $details_hashref->{$query_id}{'target_start'} - 2 - 1,
            $details_hashref->{$query_id}{'target_end'} - 1
      )    # starts from 0
   );
  }
  else {
   $scaffold_subseq = substr(
                              $scaffold_seq,
                              $details_hashref->{$query_id}{'target_start'},
                              $details_hashref->{$query_id}{'target_end'} + 2
   );
  }
  if ( $no_single_exon
       && length($scaffold_subseq) <=
       length( $details_hashref->{$query_id}{'target_seq'} ) )
  {
   $failed_cutoff++;
   print LOG "$query_id rejected: single coding exon\n";
   next;
  }

  my $number_of_ns = $scaffold_subseq =~ tr/N/N/;
  if ( $number_of_ns > ( length($scaffold_subseq) * 0.40 ) ) {
   $failed_cutoff++;
   print LOG
"$query_id rejected: genome sequence is more than 40% Ns ($number_of_ns)\n";
   next;
  }
  print LOG "$query_id accepted\n";
  print VULGAR $vulgar_data_hashref->{$query_id}
    if $vulgar_data_hashref->{$query_id};
  $pass_check{ 's'
     . $details_hashref->{$query_id}{'target_start'} . 'e'
     . $details_hashref->{$query_id}{'target_end'} } = $query_id;
  print PASSP ">$query_id\n"
    . $details_hashref->{$query_id}{'query_seq'} . "\n";

#    print PASSG ">$query_id "     . $details_hashref->{$query_id}{'target_id'} . " genome sequence\n$scaffold_subseq\n" unless $nodataprint;
  print PASSGC ">$query_id "
    . $details_hashref->{$query_id}{'target_id'}
    . " genome CDS\n"
    . $details_hashref->{$query_id}{'target_seq'} . "\n"
    unless $nodataprint;
  $accepted{$query_id} = $details_hashref->{$query_id}{'source'};
  $pass_counter++;
 }
 print "\r$counter/$query_number    $pass_counter passed           \n";
 close OUT   unless $nodataprint;
 close PASSP unless $nodataprint;

 #  close PASSG unless $nodataprint;
 close VULGAR;

 my $number_of_passing_genes = scalar( keys %accepted );
 open( PASS, ">$basename_exonerate.passed" );
 my @shuffled_genes = shuffle( keys %accepted );
 foreach my $gene (@shuffled_genes) {
  print PASS $gene . "\t" . $accepted{$gene} . "\n";
 }
 close PASS;

 my $orig_sep = $/;
 $/ = &get_gff_delimiter($exonerate_file);
 open( IN,  $exonerate_file );
 open( OUT, ">$basename_exonerate.golden" );

GENE: while ( my $record = <IN> ) {
  chomp($record);
  next if !$record || $record =~ /^\s*$/;
  my @record_data = split( "\n", $record );
  my $gene_line   = shift @record_data;
  my $mrna_line   = shift @record_data;
  my @gene_data   = split( "\t", $gene_line );
  my @mrna_data   = split( "\t", $mrna_line );
  my ( $gene_id, $mrna_id );

  if ( $gene_data[8] =~ /ID=([^;]+)/ ) {
   $gene_id = $1;
  }
  else {
   confess "Unkown GFF format for $record\n$gene_data[8]";
  }
  if ( $mrna_data[8] =~ /ID=([^;]+)/ ) {
   $mrna_id = $1;
  }
  else {
   confess "Unkown GFF format for $record\n$mrna_data[8]";
  }
  print OUT $record . $/ if $accepted{$gene_id} || $accepted{$mrna_id};
 }
 close OUT;
 $/ = $orig_sep;

 print LOG
"Processed $counter transcripts.\nFound $number_of_passing_genes sequences passing criteria.\n\n";
 close LOG;
 print
"Processed $counter transcripts.\nFound $number_of_passing_genes sequences passing criteria.\n\n";

 #cleanup
 unlink( $exonerate_file . '.pep' );
 unlink( $exonerate_file . '.cds' );
 unlink( $exonerate_file . '.gene' );
 unlink( $exonerate_file . '.original' );

 if ($stop_after_golden) {
  print "User stop requested\n";
  exit(0);
 }
}

sub process_for_gene_prediction() {
 print "Processing for gene prediction software...\n";
 my $gff_file  = shift;
 my $pass_file = shift;
 &order_fasta( $genome_sequence_file, $gff_file );

 my ( %accepted, %training_genes, %aug_optimization_genes );
 open( PASS, $pass_file ) || die "Cannot find pass file $pass_file";
 while ( my $ln = <PASS> ) {
  chomp($ln);
  my @data = split( "\t", $ln );
  next unless $data[1];
  $accepted{ $data[0] } = $data[1];
 }
 close PASS;

 print "Creating golden set using main GFF3 file $gff_file\n";

 # using fathom for a check 1
 &filter_gff( $gff_file, \%accepted, "$gff_file.fathom.gff3",
              "$gff_file.fathomnok.gff3" );
 &order_fasta( $genome_sequence_file, "$gff_file.fathom.gff3" );
 &gff2zff( "$gff_file.fathom.gff3", "$gff_file.fathom.zff" );
 my @failed =
`$fathom_exec $gff_file.fathom.zff $gff_file.fathom.gff3.fasta -validate 2>&1|grep 'error\b'`;
 if (@failed) {
  print "Fathom tells us that "
    . scalar(@failed)
    . " annotations failed... Taking them into account...\n";
  foreach my $fail (@failed) {
   $fail =~ /^(\S+)/;
   my $id = $1;
   delete( $accepted{$id} );
  }
  &filter_gff( $gff_file, \%accepted, "$gff_file.fathom.gff3",
               "$gff_file.fathomnok.gff3" );
  &order_fasta( $genome_sequence_file, "$gff_file.fathom.gff3" );
  &gff2zff( "$gff_file.fathom.gff3", "$gff_file.fathom.zff" );
 }

 # using fathom for a check 2
 my @warnings =
`$fathom_exec $gff_file.fathom.zff $gff_file.fathom.gff3.fasta -validate 2>&1|grep warnings|egrep -v 'GC\.\.AG|short'`;
 open( FATHOM, ">$gff_file.zff.warnings" );
 print FATHOM @warnings;
 close FATHOM;
 if ( scalar(@warnings) > 1 ) {
  pop(@warnings);
  print "Fathom gave us "
    . scalar(@warnings)
    . " warnings (possible non-canonical splice sites, incomplete CDSs (if -complete given) etc, see $gff_file.zff.warnings). Fixing...\n";
  $failed_cutoff += scalar(@warnings);
  foreach my $fail (@warnings) {
   $fail =~ /^\S+:\s+(\S+)/;
   my $id = $1;
   if ( !$only_complete && $fail =~ /cds:incomplete/ ) {
    next;
   }
   delete( $accepted{$id} );
  }
 }
 unlink("$gff_file.fathomnok.gff3");
 unlink("$gff_file.fathom.zff");
 unlink("$gff_file.fathom.zff.xdef");
 unlink("$gff_file.fathom.gff3");
 unlink("$gff_file.fathom.gff3.fasta");

 my $number_of_passing_genes      = scalar( keys %accepted );
 my $passed_genes_in_training     = int(0);
 my $passed_genes_in_optimization = int(0);

 $training_set_size = int( $number_of_passing_genes * 0.40 )
   if !$training_set_size
    || $training_set_size > int( $number_of_passing_genes * 0.40 );
 foreach my $gene ( shuffle( keys %accepted ) ) {
  $training_genes{$gene} = $accepted{$gene};
  $passed_genes_in_training++;
  last if $passed_genes_in_training >= $training_set_size;
 }
 foreach my $gene ( shuffle( keys %accepted ) ) {
  next if $training_genes{$gene};
  $aug_optimization_genes{$gene} = $accepted{$gene};
  $passed_genes_in_optimization++;
  last if $passed_genes_in_optimization >= $aug_optimization_geneset;
 }
 $aug_optimization_geneset = scalar( keys %aug_optimization_genes );
 $training_set_size        = scalar( keys %training_genes );
 die "Augustus optimization gene set is zero!\n"
   if $aug_optimization_geneset < 1;
 die "Augustus Training gene set is zero!\n" if $training_set_size < 1;

 print "'Golden': $number_of_passing_genes. Filtering GFFs...\n";
 &filter_gff( $gff_file, \%accepted, "$gff_file.golden.gff3" );
 push(@to_evaluate,"$gff_file.golden.gff3");
#    . "Training set: $training_set_size. Test: "
#    . ( $number_of_passing_genes - $training_set_size ) . ".\n"
#    . "Augustus optimization set: "
#    . $aug_optimization_geneset . ".\n"
#    . " (as the optimization set is taken from the test set,\n"
#    . " the Augustus test set is "
#    . ( $number_of_passing_genes - $training_set_size - $aug_optimization_geneset )
#    . ". more or less \n"

 print
"Processing for gene predictors... - Specific filters may reject some more genes during checking/conversion\n";
 print "\nProcessing for:\n";

 #augustus
 print "\tAugustus\n";

 #for training
 &gff2hints( "$gff_file.golden.gff3", 1 )
   if -s "$gff_file.golden.gff3";
 &filter_gff( "$gff_file.golden.gff3", \%training_genes,
              "$gff_file.golden.train.gff3" );
 &filter_gff( "$gff_file.golden.train.gff3.rest",
              \%aug_optimization_genes, "$gff_file.golden.optimization.gff3" );
 rename( "$gff_file.golden.optimization.gff3.rest",
         "$gff_file.golden.test.gff3" );

 if (    $gff2gb_exec
      && $augustus_train_exec
      && $augustus_filterGenes_exec )
 {

  system(
"$gff2gb_exec $gff_file.golden.train.gff3 $genome_sequence_file $augustus_flank_region $gff_file.golden.train.gb >/dev/null 2> /dev/null"
  ) if -s "$gff_file.golden.train.gff3";
  system(
"$augustus_train_exec --species=generic $gff_file.golden.train.gb 2>&1 | grep 'n sequence' | perl -pe 's/.*n sequence (\\S+):.*/\$1/' | sort -u > $gff_file.golden.train.gb.bad.lst 2>/dev/null"
  ) if -s "$gff_file.golden.train.gb";

  if (-s "$gff_file.golden.train.gb.bad.lst"){
    system("$augustus_filterGenes_exec $gff_file.golden.train.gb.bad.lst $gff_file.golden.train.gb > $gff_file.golden.train.good.gb 2>/dev/null");
  }else{
    symlink("$gff_file.golden.train.gb",$gff_file.golden.train.good.gb") if -s "$gff_file.golden.train.gb";
  }

  system(
"$gff2gb_exec $gff_file.golden.test.gff3 $genome_sequence_file $augustus_flank_region $gff_file.golden.test.gb  >/dev/null 2>/dev/null"
  ) if -s "$gff_file.golden.test.gff3";
  system(
"$augustus_train_exec --species=generic $gff_file.golden.test.gb 2>&1 | grep 'n sequence' | perl -pe 's/.*n sequence (\\S+):.*/\$1/' | sort -u > $gff_file.golden.test.gb.bad.lst 2>/dev/null"
  ) if -s "$gff_file.golden.test.gb";

  if (-s "$gff_file.golden.test.gb.bad.lst"){
    system("$augustus_filterGenes_exec $gff_file.golden.test.gb.bad.lst $gff_file.golden.test.gb > $gff_file.golden.test.good.gb 2>/dev/null");
  }else{
    symlink("$gff_file.golden.test.gb","$gff_file.golden.test.good.gb") if -s "$gff_file.golden.test.gb";
  }

  system(
"$gff2gb_exec $gff_file.golden.optimization.gff3 $genome_sequence_file $augustus_flank_region $gff_file.golden.optimization.gb  >/dev/null"
  ) if -s "$gff_file.golden.optimization.gff3";
  system(
"$augustus_train_exec --species=generic $gff_file.golden.optimization.gb 2>&1 | grep 'n sequence' | perl -pe 's/.*n sequence (\\S+):.*/\$1/' | sort -u > $gff_file.golden.optimization.gb.bad.lst 2>/dev/null"
  ) if -s "$gff_file.golden.optimization.gb";

  if (-s "$gff_file.golden.optimization.gb.bad.lst"){
     system("$augustus_filterGenes_exec $gff_file.golden.optimization.gb.bad.lst $gff_file.golden.optimization.gb > $gff_file.golden.optimization.good.gb 2>/dev/null");
  }else{
     symlink("$gff_file.golden.optimization.gb","$gff_file.golden.optimization.good.gb");
  }
 }

 #parse GB
 my $train_ref = &parse_gb("$gff_file.golden.train.good.gb")
   if -s "$gff_file.golden.train.good.gb";
 my $test_ref = &parse_gb("$gff_file.golden.test.good.gb")
   if -s "$gff_file.golden.test.good.gb";
 my $opt_ref = &parse_gb("$gff_file.golden.optimization.good.gb")
   if -s "$gff_file.golden.optimization.good.gb";

 #hints
 &gff2hints("$gff_file")
   if -s "$gff_file";
 &gff2hints("$gff_file.golden.gff3.rest")
   if -s "$gff_file.golden.gff3.rest";

 #geneid
 print "\tgeneid\n";
 my $geneid_gff_file = &create_geneid_from_gff($final_gff);
 if ( -s $geneid_gff_file ) {
  &confirm_geneid_gff( $geneid_gff_file, $geneid_gff_file . '.nr' )
    if $geneid_gff_file && -s $geneid_gff_file;
  if ( -s $geneid_gff_file . '.nr' ) {
   $geneid_gff_file .= '.nr';
   &filter_gff( $geneid_gff_file, \%accepted, "$geneid_gff_file.golden.gff3" );
   &filter_gff( "$geneid_gff_file.golden.gff3",
                \%training_genes, "$geneid_gff_file.golden.train.gff3" );
   rename( "$geneid_gff_file.golden.train.gff3.rest",
           "$geneid_gff_file.golden.test.gff3" );

   &order_fasta( $genome_sequence_file, "$geneid_gff_file.golden.gff3" );
   &order_fasta( $genome_sequence_file, "$geneid_gff_file.golden.train.gff3" );
   &order_fasta( $genome_sequence_file, "$geneid_gff_file.golden.test.gff3" );

#TODO run optimization script (already written)
# for i in {0..20}; do INDEX_PADDED=`printf "%02d" $i`; ParaFly -CPU 1 -v -c H.armigera_geneid.optimization.$INDEX_PADDED -failed_cmds H.armigera_geneid.optimization.$INDEX_PADDED.failed ; done
#TODO get intron/exon sizes with R
#data_intron <- read.table(file="exonerate.results.golden.intron.gff3",header=F,sep="\t")
# data_exon <- read.table(file="exonerate.results.golden.exon.gff3",header=F,sep="\t")
#dataInter<-read.table(file="intergenic.exon.gff",header=F,sep="\t");
#data_intron_length<-abs(data_intron$V5-data_intron$V4)
#data_exon_length<-abs(data_exon$V5-data_exon$V4)
#dataInter_length<-abs(dataInter$V5-dataInter$V4)
#quantile(dataInter_length,probs = c(0,0.25,0.5,0.95,0.96,0.97,0.98,0.99,1)) # pick 98%
# quantile(data_intron_length,probs = c(0,0.25,0.5,0.95,0.96,0.97,0.98,0.99,1)) # pick 98%
# INTRAgenic connections
# First+:Internal+                Internal+:Terminal+             15:2000 block
# Terminal-:Internal-             First-:Internal-                15:2000 blockr
# INTERgeneic connections
# aataaa+:Terminal+:Single+       Single+:First+:Promoter+        150:Infinity
# aataaa+:Terminal+:Single+       Single-:Terminal-:aataaa-       150:Infinity
# Promoter-:First-:Single-        Single+:First+:Promoter+        150:Infinity
# Promoter-:First-:Single-        Single-:Terminal-:aataaa-       150:Infinity
   &gb2geneid( $train_ref, "$geneid_gff_file.golden.train.good.gb.geneid" )
     if $train_ref;
   &gb2geneid( $test_ref, "$geneid_gff_file.golden.test.good.gb.geneid" )
     if $test_ref;
   &gb2geneid( $opt_ref, "$geneid_gff_file.golden.optimization.good.gb.geneid" )
     if $opt_ref;
   &order_fasta( "$gff_file.golden.optimization.good.gb.fasta",
                 "$geneid_gff_file.golden.optimization.good.gb.geneid" )
     if -s "$geneid_gff_file.golden.optimization.good.gb.geneid";
   &order_fasta( "$gff_file.golden.train.good.gb.fasta",
                 "$geneid_gff_file.golden.train.good.gb.geneid" )
     if -s "$geneid_gff_file.golden.train.good.gb.geneid";
   &order_fasta( "$gff_file.golden.test.good.gb.fasta",
                 "$geneid_gff_file.golden.test.good.gb.geneid" )
     if -s "$geneid_gff_file.golden.test.good.gb.geneid";
  }
 }

 print "\tsnap\n";

 #snap
 &gff2zff( "$gff_file.golden.gff3",       "$gff_file.golden.zff" );
 &gff2zff( "$gff_file.golden.train.gff3", "$gff_file.golden.train.zff" );
 &gff2zff( "$gff_file.golden.test.gff3",  "$gff_file.golden.test.zff" );
 &gff2zff( "$gff_file.golden.optimization.gff3",
           "$gff_file.golden.optimization.zff" );

 #glimmer
 print "\tglimmer\n";
 &gb2glimmer( $train_ref, "$gff_file.golden.train.good.gb.glimmer" )
   if $train_ref;
 &gb2glimmer( $test_ref, "$gff_file.golden.test.good.gb.glimmer" )
   if $test_ref;
 &gb2glimmer( $opt_ref, "$gff_file.golden.optimization.good.gb.glimmer" )
   if $opt_ref;

 #Different thresholds that can be chosen for the splice sites
 #can be consulted in: - false.acc  false.don  false.atg

 print "\tevaluation\n";
 foreach my $file (@to_evaluate){
  &gff_to_gtf($file);
 }
}

sub revcomp {
 my $dna     = shift;
 my $revcomp = reverse($dna);
 $revcomp =~ tr/ACGTacgt/TGCAtgca/;
 return $revcomp;
}

sub filter_gff() {
 my $gff      = shift;
 my $hash_ref = shift;

 return if !$gff || !$hash_ref || !-s $gff;

 my $filter_out = shift;
 $filter_out = "$gff.filtered" if !$filter_out;
 my $filter_out2 = shift;
 $filter_out2 = $filter_out . ".rest" if !$filter_out2;

 my $delimiter = &get_gff_delimiter($gff);

 my $orig_sep = $/;
 $/ = $delimiter;

 open( GFF, $gff ) || die;
 open( OUT, ">$filter_out" );
 open( OUT2, ">$filter_out2" );
 my $rejects;

GENE: while ( my $record = <GFF> ) {
  chomp($record);
  next if !$record || $record =~ /^\s*$/;
  my ( $gene_id, $mrna_id );
  my @record_data = split( "\n", $record );
  my @gene_data   = split( "\t", $record_data[0] );
  if ( $gene_data[8] && $gene_data[8] =~ /ID=([^;]+)/ ) {
   $gene_id = $1;
   my @mrna_data = split( "\t", $record_data[1] );
   if ( $mrna_data[8] && $mrna_data[8] =~ /ID=([^;]+)/ ) {
    $mrna_id = $1;
   }
  }
  else {

   # gene_id
   $gene_data[8] =~ s/\.mRNA$//;
   $gene_id = $gene_data[8];
  }
  confess "No gene ID for this record!\n$record\n" unless $gene_id;
  my $source = $gene_data[1];

  if ( ( $gene_id && $hash_ref->{$gene_id} && $hash_ref->{$gene_id} eq $source )
    || ( $mrna_id && $hash_ref->{$mrna_id} && $hash_ref->{$mrna_id} eq $source )
    )
  {
   print OUT $record . $delimiter;
  }
  else {
   print OUT2 $record . $delimiter;
   $rejects++;
  }
 }

 close GFF;
 close OUT;
 close OUT2;
 warn "Warning: filtered file $filter_out is empty...\n"  if !-s $filter_out;
 warn "Warning: filtered file $filter_out2 is empty...\n" if !-s $filter_out2;
 $/ = $orig_sep;

 return $rejects;
}

sub gff2zff() {

=cut

        None,       /* 0 A non-feature? */
        
        Inter,      /* 1 Intergenic */
        Int0,       /* 2 phase 0 intron */
        Int1,       /* 3 phase 1 intron */
        Int1T,      /* 4 phase 1 intron, previous exon has overhanging T */
        Int2,       /* 5 phase 2 intron */
        Int2TA,     /* 6 phase 2 intron, previous exon has overhanging TA */
        Int2TG,     /* 7 phase 2 intron, previous exon has overhanging TG */
        Intron,     /* 8 intron */

        UTR5,       /* 9 5' UTR */
        UTR3,       /* 10 3' UTR */
        
        Esngl,      /* 11 single exon gene */
        Einit,      /* 12 initial exon */
        Eterm,      /* 13 terminal exon */
        Exon,       /* 14 generic or internal exon */
        Coding,     /* 15 generically coding */
        Gene,       /* 16 geneircally gene */
        
        Acceptor,   /* 17 splice acceptor */
        Donor,      /* 18 splice donor */
        Start,      /* 19 ATG */
        Stop,       /* 20 stop codon */
        
        Repeat,     /* 21 Repetitive element */
        CNS,        /* 22 Conserved sequence */
        ORF,        /* 23 Open Reading Frame */
        
        PolyA,      /* 24 poly-A signal */
        Prom,       /* 25 promoter */
        BPS,        /* 26 branch point signal */
        TSS,        /* 27 Trans-splice site */
        
        Misc,       /* 28 miscellaneous feature */
        
        HSP_NN,     /* 29 high-scoring pair: nt-nt (BLASTN) */
        HSP_NA,     /* 20 high-scoring pair: nt-aa (BLASTX) */
        HSP_AN,     /* 31 high-scoring pair: aa-nn (TBLASTN) */
        HSP_AA      /* 32 high-scoring pair: aa-aa (BLASTP, TBLASTX) */

=cut

 ## NB snap uses the word exon to mean CDS... argh
 # don't use hardmasked.
 my $gff       = shift;
 my $out       = shift;
 my $out_hints = $out . ".xdef";
 open( GFF, $gff ) || die "Can't open $gff";
 open( OUT, ">$out" );
 open( XDEF, ">$out_hints" );
 my $orig_delim = $/;
 $/ = "\n\n";
 my $previous_ref;

 while ( my $record = <GFF> ) {
  my @lines = split( "\n", $record );
  chomp(@lines);
  my @gene_data = split( "\t", $lines[0] );
  my @mRNA_data = split( "\t", $lines[1] );
  my $ref       = $gene_data[0];
  my $zff_line;
  my $xdef_line;
  if ( !$previous_ref || $ref ne $previous_ref ) {
   $zff_line  .= ">$ref\n";
   $xdef_line .= ">$ref\n";
   $previous_ref = $ref;
  }
  $mRNA_data[8] =~ /ID=([^;]+)/;
  my $mRNA_id = $1 || die Dumper $record;
  my $cds_counter = 0;
  my @cds_data;
  foreach my $line (@lines) {
   next if $line =~ /^\s*$/ || $line =~ /^#/;
   my @data = split( "\t", $line );
   die "GFF $gff had a line that didn't begin with the reference $ref" . $line
     unless $data[0] eq $ref;
   next if $data[2] ne 'CDS';
   $cds_data[$cds_counter] = $line;
   $cds_counter++;
  }
  for ( my $i = 0 ; $i < $cds_counter ; $i++ ) {
   my $cds   = $cds_data[$i];
   my @data  = split( "\t", $cds );
   my $start = $data[6] eq '+' ? $data[3] : $data[4];
   my $end   = $data[6] eq '+' ? $data[4] : $data[3];
   my $type;
   if ( $cds_counter == 1 ) {
    $type = 'Esngl';

   }
   elsif ( $i == 0 && $data[6] eq '+' ) {
    $type = 'Einit';
   }
   elsif ( $i == ( $cds_counter - 1 ) && $data[6] eq '+' ) {
    $type = 'Eterm';
   }
   else {
    $type = 'Exon';
   }

   $zff_line .=
     join( "\t", ( $type, $start, $end, $mRNA_id ) )
     . "\n";    #unless $type eq 'Esngl';
   $xdef_line .= join(
                       "\t",
                       (
                         "Coding", $data[3], $data[4], $data[6],
                         '+100',   '.',      '.',      '.',
                         'ADJ'
                       )
   ) . "\n";    #unless $type eq 'Esngl';
  }
  print OUT $zff_line   if $zff_line;
  print XDEF $xdef_line if $xdef_line;
 }
 $/ = $orig_delim;
 close OUT;
 close XDEF;
 close GFF;
}

sub sort_gff3() {

 # find out why two empty lines still cause sort to crash
 my $gff       = shift;
 my $delimiter = shift;
 $delimiter = &get_gff_delimiter($gff) if !$delimiter;
 my $orig_sep = $/;
 open( GFF, $gff ) || confess "Cannot open $gff $!";
 open( OUT, ">$gff.sorted" );
 $/ = $delimiter;
 my @records;

 while ( my $rec = <GFF> ) {
  chomp($rec);
  next if !$rec || $rec =~ /^\s*$/;
  push( @records, $rec );
 }
 close GFF;

 foreach my $rec (
  sort {
   my @first  = split( "\n", $a );
   my @second = split( "\n", $b );

   my @split1 = split( "\t", $first[0] );
   my @split2 = split( "\t", $second[0] );
   return $split1[0] cmp $split2[0] || $split1[3] <=> $split2[3];
  } @records
   )
 {
  print OUT $rec . $/;
 }

 close OUT;
 $/ = $orig_sep;
 rename( "$gff.sorted", $gff );
}

sub order_fasta() {
 my $fasta = shift || die;
 my $gff   = shift || die;
 confess "Cannot find $fasta " . $! unless -s $fasta;
 confess "Cannot find $gff " . $!   unless -s $gff;
 my $fasta_hashref = &read_fasta($fasta);

 my @order = `cut -f 1 $gff|uniq `;
 open( OUT, ">$gff.fasta" );
 my %done;

 foreach my $ln (@order) {
  next if $ln =~ /^#/ || $ln =~ /^\s*$/;
  $ln =~ /^(\S+)/;
  my $id = $1;
  next unless $id;
  if ( !$done{$id} ) {
   my $seq = $fasta_hashref->{$id};
   unless ($seq) {
    die "Cannot get $id from $fasta for $gff\n";
   }
   print OUT ">$id\n$seq\n";
   $done{$id}++;
  }
 }
 close OUT;
 return "$gff.fasta";
}

sub confirm_geneid_gff() {
 my $file = shift;
 my $out  = shift;

 my %ignore;
 my $orig_sep = $/;
 $/ = '#$' . "\n";

 open( IN,   $file );
 open( OUT,  ">$out" );
 open( LOG2, ">$out.log" );

GENE: while ( my $record = <IN> ) {
  chomp($record);
  next if !$record || $record =~ /^\s*$/;
  my @record_data = split( "\n", $record );
  my @gene_data   = split( "\t", $record_data[0] );
  my $gene_id     = $gene_data[8];
  $gene_id =~ s/\.mRNA$//;
  next if $ignore{$gene_id};
  my %hash;
  $hash{'checkF'} = 0;
  $hash{'checkT'} = 0;

  foreach my $line (@record_data) {
   my @data = split( "\t", $line );
   if ( $data[2] eq 'Single' ) {
    print LOG2 "Excluding singleton $gene_id\n";
    $ignore{$gene_id} = 1;
    next GENE;
   }

# fixing various bugs with double first/last exons (and also more than one align per gene):
   $hash{'checkF'}++ if $data[2] eq 'First';
   $hash{'checkT'}++ if $data[2] eq 'Terminal';

  }

  if ( !$hash{'checkF'} || $hash{'checkF'} > 1 ) {
   print LOG2 "Gene $gene_id has 0 or more than 1 first exons ("
     . $hash{'checkF'}
     . "). Excluding gene\n";
   $ignore{$gene_id} = 1;
   next GENE;
  }
  if ( !$hash{'checkT'} || $hash{'checkT'} > 1 ) {
   print LOG2 "Gene $gene_id has 0 or more than 1 terminal exons("
     . $hash{'checkT'}
     . "). Excluding gene\n";
   $ignore{$gene_id} = 1;
   next GENE;
  }
  print OUT $record . '#$' . "\n";
 }

 close IN;
 close OUT;
 close LOG2;
 print "\tGeneID checked. Ignored "
   . scalar( keys %ignore )
   . " genes. See $out.log for details\n";

 $/ = $orig_sep;
}

sub remove_overlapping_gff() {

 # we allow genes on opposite strands
 my $file               = shift;
 my $out                = shift;
 my $delimiter          = shift;
 my $optional_pass_file = shift;

 $delimiter = &get_gff_delimiter($file) if !$delimiter;

 &sort_gff3( $file, $delimiter );

 my ( %overlap_found, %master_gene_list, %cds_sizes );
 my $skipped = int(0);
 my ( $previous_ref, $previous_mrna_id, %last_position_check );

 open( LOG2, ">$out.log" );
 my $orig_sep = $/;
 $/ = $delimiter;

 open( IN, $file );

GENE: while ( my $record = <IN> ) {
  chomp($record);
  my ( $gene_id, $mrna_id );
  next if !$record || $record =~ /^\s*$/;
  my @record_data = split( "\n", $record );
  my @gene_data   = split( "\t", $record_data[0] );
  confess "Weird GFF for $record" unless $gene_data[8];
  my $ref_id = $gene_data[0];
  my $strand = $gene_data[6];

  if ( $gene_data[8] =~ /ID=([^;]+)/ ) {
   $gene_id = $1;
  }
  else {

   # gene_id
   $gene_data[8] =~ s/\.mRNA$//;
   $gene_id = $gene_data[8];
  }
  my @mrna_data = split( "\t", $record_data[1] );
  if ( $mrna_data[8] && $mrna_data[8] =~ /ID=([^;]+)/ ) {
   $mrna_id = $1;
   $mrna_id =~ s/\.mRNA$//;
  }
  else {
   $mrna_id = $gene_id;
  }
  confess "No gene ID for this record!\n$record\n" unless $gene_id;

  unless ($mrna_id) {
   print LOG2 "No mRNA found for gene $mrna_id skipped\n";
   $skipped++;
   next GENE;
  }
  if ( $master_gene_list{$mrna_id} ) {
   print LOG2
"Gene $mrna_id found more than once. Skipping new data and keeping what I found first\n";
   $skipped++;
   next GENE;
  }

  foreach my $d (@record_data) {
   my @t = split( "\t", $d );
   if ( $t[2] eq 'CDS' ) {
    $cds_sizes{$mrna_id} += abs( $t[4] - $t[3] );
   }
  }
  if ( !$cds_sizes{$mrna_id} ) {
   print LOG2 "Non-coding transcript $mrna_id skipped\n";
   $skipped++;
   next GENE;
  }

  my ( $smallest_coord, $largest_coord ) =
    sort { $a <=> $b } ( $gene_data[3], $gene_data[4] );

  if (    $previous_mrna_id
       && $previous_ref
       && $last_position_check{$ref_id}{$strand}
       && $previous_ref eq $ref_id
       && $smallest_coord <= $last_position_check{$ref_id}{$strand} )
  {
   print LOG2 "Gene overlapping $previous_mrna_id found: $mrna_id\n";
   $overlap_found{$previous_mrna_id}{$mrna_id} = $record . $delimiter;
   next;
  }
  $master_gene_list{$mrna_id}            = $record . $delimiter;
  $previous_mrna_id                      = $mrna_id;
  $previous_ref                          = $ref_id;
  $last_position_check{$ref_id}{$strand} = $largest_coord;
 }
 $/ = $orig_sep;
 close IN;
 print LOG2 "\nDeciding on overlaps:\n";

 my %kept;

 open( OUT, ">$out" );
 foreach my $gene ( sort keys %master_gene_list ) {
  if ( $overlap_found{$gene} ) {
   my $master_size = $cds_sizes{$gene};
   my $gene_longest;
   foreach my $overlap_gene ( keys %{ $overlap_found{$gene} } ) {
    my $size = $cds_sizes{$overlap_gene};
    $gene_longest = $overlap_gene if $size > $master_size;
   }
   if ( !$gene_longest ) {
    print LOG2 "$gene overlaps: $gene kept as longest CDS\n";
    print OUT $master_gene_list{$gene};
    my @lines = split( "\n", $master_gene_list{$gene} );
    my @data  = split( "\t", $lines[0] );
    $kept{$gene} = $data[1];
   }
   else {
    print LOG2 "$gene overlaps: longer CDS $gene_longest found\n";
    print OUT $overlap_found{$gene}{$gene_longest};
    my @lines = split( "\n", $overlap_found{$gene}{$gene_longest} );
    my @data  = split( "\t", $lines[0] );
    $kept{$gene_longest} = $data[1];
   }
   $skipped += scalar( keys %{ $overlap_found{$gene} } );
  }
  else {

   # no overlaps
   print OUT $master_gene_list{$gene};
   my @lines = split( "\n", $master_gene_list{$gene} );
   my @data  = split( "\t", $lines[0] );
   $kept{$gene} = $data[1] if !$kept{$gene};    # only the first one
  }

 }
 close OUT;
 close LOG2;
 &sort_gff3($out);
 my %passed;
 if ( $optional_pass_file && -s $optional_pass_file ) {
  open( IN, $optional_pass_file );
  while ( my $ln = <IN> ) {
   chomp($ln);
   my @data = split( "\t", $ln );
   if ( $kept{ $data[0] } && $kept{ $data[0] } eq $data[1] ) {
    $passed{ $data[0] } = $data[1];
   }
  }
  close IN;

  open( OUT, ">$optional_pass_file" );
  foreach my $id ( keys %passed ) {
   print OUT $id . "\t" . $passed{$id} . "\n";
  }
  close OUT;

 }

 print "\tOverlaps checked. Skipped $skipped genes. See $out.log for details\n";
 return $skipped;
}

sub parse_gb() {
 use Bio::SeqIO;
 my $file      = shift;
 my $fasta_out = $file . ".fasta";
 my %hash;
 my $gb_obj = Bio::SeqIO->new( -file => $file, -format => 'genbank' );
 while ( my $seq_obj = $gb_obj->next_seq() ) {
  my ( $exon_count, $cds_strand );
  my $id  = $seq_obj->id();
  my $seq = $seq_obj->seq();
  foreach my $feat ( $seq_obj->get_SeqFeatures() ) {
   if (    $feat->location->isa('Bio::Location::SplitLocationI')
        && $feat->primary_tag eq 'CDS' )
   {
    $cds_strand = $feat->strand;
    foreach my $loc ( $feat->location->sub_Location ) {
     $exon_count++;
     push( @{ $hash{$id}{'exons'}{$exon_count} }, ( $loc->start, $loc->end ) );
    }
   }
  }
  $hash{$id}{'seq'}    = $seq;
  $hash{$id}{'strand'} = $cds_strand;
 }

 # have to because geneid demands it in order
 open( FSAOUT, ">$fasta_out" );
 foreach my $id ( sort keys %hash ) {
  print FSAOUT ">$id\n" . $hash{$id}{'seq'} . "\n";
 }
 close FSAOUT;
 &process_cmd("$cdbfasta_exec $fasta_out 2>/dev/null");
 return \%hash;
}

sub gb2gff3() {
 my $hash_ref = shift;
 my $out      = shift;
 open( GFF, ">$out" );
 foreach my $id ( sort keys %{$hash_ref} ) {
  next unless $hash_ref->{$id}{'strand'};
  my @sorted_exons =
    sort { $a <=> $b } keys %{ $hash_ref->{$id}{'exons'} };
  if ( $hash_ref->{$id}{'strand'} == 1 ) {
   my $strand = '+';
   my ( $gstart, $e ) = @{ $hash_ref->{$id}{'exons'}->{ $sorted_exons[0] } };
   my ( $s, $gend ) = @{ $hash_ref->{$id}{'exons'}->{ $sorted_exons[-1] } };
   print GFF $id
     . "\tGB\texon\t"
     . $gstart . "\t"
     . $gend
     . "\t.\t$strand\t.\t$id.1" . "\n";

   foreach my $exon_counter (@sorted_exons) {
    my ( $start, $stop ) = @{ $hash_ref->{$id}{'exons'}->{$exon_counter} };
    print GFF $id
      . "\tGB\texon\t"
      . $start . "\t"
      . $stop
      . "\t.\t$strand\t.\t$id.1" . "\n";
   }
  }
  elsif ( $hash_ref->{$id}{'strand'} == -1 ) {
   my $strand = '-';
   my ( $gstart, $e ) = @{ $hash_ref->{$id}{'exons'}->{ $sorted_exons[-1] } };
   my ( $s, $gend ) = @{ $hash_ref->{$id}{'exons'}->{ $sorted_exons[0] } };
   print GFF $id
     . "\tGB\texon\t"
     . $gstart . "\t"
     . $gend
     . "\t.\t$strand\t.\t$id.1" . "\n";
   foreach my $exon_counter (@sorted_exons) {
    my ( $start, $stop ) = @{ $hash_ref->{$id}{'exons'}->{$exon_counter} };
    print GFF $id
      . "\tGB\texon\t"
      . $start . "\t"
      . $stop
      . "\t.\t$strand\t.\t$id.1" . "\n";
   }
  }
  print GFF "\n\n";
 }
 close GFF;
 &sort_gff3($out);
}

sub gb2geneid() {
 my $hash_ref = shift;
 my $out      = shift;
 open( GENEID, ">$out" );
 foreach my $id ( sort keys %{$hash_ref} ) {
  next unless $hash_ref->{$id}{'strand'};
  if ( $hash_ref->{$id}{'strand'} == 1 ) {
   my $strand = '+';
   foreach my $exon_counter (
                              sort { $a <=> $b }
                              keys %{ $hash_ref->{$id}{'exons'} }
     )
   {
    my $exon_type;
    if ( $exon_counter == 1
         && ( scalar( keys %{ $hash_ref->{$id}{'exons'} } ) == 1 ) )
    {
     $exon_type = "Single";
    }
    elsif ( $exon_counter == 1 ) {
     $exon_type = "First";
    }
    elsif ( $exon_counter == ( scalar( keys %{ $hash_ref->{$id}{'exons'} } ) ) )
    {
     $exon_type = 'Terminal';
    }
    else {
     $exon_type = 'Internal';
    }
    my ( $start, $stop ) = @{ $hash_ref->{$id}{'exons'}->{$exon_counter} };
    print GENEID $id
      . "\tforGeneID\t$exon_type\t"
      . $start . "\t"
      . $stop
      . "\t.\t$strand\t.\t$id.1" . "\n";
   }
  }
  elsif ( $hash_ref->{$id}{'strand'} == -1 ) {
   my $strand = '-';
   foreach my $exon_counter (
                              sort { $a <=> $b }
                              keys %{ $hash_ref->{$id}{'exons'} }
     )
   {
    my $exon_type;
    if ( $exon_counter == 1
         && ( scalar( keys %{ $hash_ref->{$id}{'exons'} } ) == 1 ) )
    {
     $exon_type = "Single";
    }
    elsif ( $exon_counter == 1 ) {
     $exon_type = "Terminal";
    }
    elsif ( $exon_counter == ( scalar( keys %{ $hash_ref->{$id}{'exons'} } ) ) )
    {
     $exon_type = 'First';
    }
    else {
     $exon_type = 'Internal';
    }
    my ( $start, $stop ) = @{ $hash_ref->{$id}{'exons'}->{$exon_counter} };
    print GENEID $id
      . "\tforGeneID\t$exon_type\t"
      . $start . "\t"
      . $stop
      . "\t.\t$strand\t.\t$id.1" . "\n";
   }
  }
  print GENEID '#$' . "\n";
 }
 close GENEID;
 &sort_gff3($out);
}

sub gb2glimmer() {
 my $hash_ref = shift;
 my $out      = shift;
 open( GLIMMER, ">$out" );
 foreach my $id ( sort keys %{$hash_ref} ) {
  next unless $hash_ref->{$id}{'strand'};
  if ( $hash_ref->{$id}{'strand'} == 1 ) {
   foreach my $exon (
                      sort { $a <=> $b }
                      keys %{ $hash_ref->{$id}{'exons'} }
     )
   {
    my ( $start, $stop ) = @{ $hash_ref->{$id}{'exons'}->{$exon} };
    print GLIMMER $id . " " . $start . " " . $stop . "\n";
   }
  }
  elsif ( $hash_ref->{$id}{'strand'} == -1 ) {
   foreach my $exon (
                      sort { $b <=> $a }
                      keys %{ $hash_ref->{$id}{'exons'} }
     )
   {
    my ( $start, $stop ) = @{ $hash_ref->{$id}{'exons'}->{$exon} };
    print GLIMMER $id . " " . $stop . " " . $start . "\n";
   }
  }
  print GLIMMER "\n";
 }
 close GLIMMER;
}

sub evaluate_geneid() {

 #augustus gb file to fsa + gff_geneid flavour with bioperl
}

sub optimize_geneid() {

=cut
bash script to convert to perl
 IeWF='-4.5'                              # Initial Weigth First     value -7
deWF='0.5'                               # Delta Weigth First       value 0.5
FeWF='-1.5'                              # Final Weigth First       value -1.5
IeWFini=$IeWF                          # Initial value

IoWF='0.20'                              # Initial Weigth First     value 0.30
doWF='0.05'                             # Delta Weigth First       value 0.05
FoWF='0.45'                              # Final Weigth First       value 0.70


GENEID=$HOME/software/geneid/bin/
EVAL=$HOME/software/geneid/Evaluation/bin/

PARAM=$1     # Param file for geneid
SEQFILE=$2
CDSFILE=$3

if [[ ! $PARAM || ! $SEQFILE || ! $CDSFILE || ! -f $PARAM || ! -f $SEQFILE || ! -f $CDSFILE ]];then
        echo No input! Give: parameter_file genome_fasta gff_annotations
        exit 255
fi
 
 echo "
SN  = sensitivity nucleotide level
SP  = specificity nucleotide level
CC  = correlation SN/SP
SNe  = sensitivity exon level
SPe  = specificity exon level
SNSP  = correlation SNe/Spe
raME  = ratio missing exons
raWE  = ratio wrong exons
SNSPg = correlation SNg/SPg
raMG  = ratio missing genes
raWG  = ratio wrong genes

oWF     eWF     SN      SP      CC      SNe     SPe     SNSP    raME    raWE    raMG    raWG" > header

ef=`gawk "BEGIN{print ($IeWF <= $FeWF) ? 1 : 0}"`

of=`gawk "BEGIN{print ($IoWF <= $FoWF) ? 1 : 0}"`


while [ $of -ne 0 ]
do
        IeWF=$IeWFini
        ef=`gawk "BEGIN{print ($IeWF <= $FeWF) ? 1 : 0}"`

        while [ $ef -ne 0 ]
        do

echo "             gawk '{if (NR==27)  {print 1-$IoWF+0.2, 1-$IoWF-0.1, 1-$IoWF-0.1, 1-$IoWF+0.2;} else if (NR==30)  {print $IoWF, $IoWF, $IoWF, $IoWF;} else if (NR==36) {print $IeWF, $IeWF, $IeWF, $IeWF;} else  
print}' $PARAM > /tmp/tmp.$$.$IoWF.$IeWF   "

echo "           $GENEID/geneid -GP /tmp/tmp.$$.$IoWF.$IeWF $SEQFILE |egrep -v 'exon|^#' |sort +0 -1 +3n | gawk '{if (NR==1) ant=\$1; if (\$1!=ant) {print \"#\$\";ant=\$1}; print }' > /tmp/Predictions.$$.$IoWF.$I
eWF.gff"        

echo "       $EVAL/evaluation -sta /tmp/Predictions.$$.$IoWF.$IeWF.gff $CDSFILE | tail -2 | head -1 |  gawk '{printf \"%6.2f    %6.2f   %6.2f   %6.2f   %6.2f   %6.2f   %6.2f   %6.2f   %6.2f   %6.2f   %6.2f   %6.2
f\\n\", $IoWF, $IeWF, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$12, \$13}' > result_$IoWF.$IeWF.txt "

echo "              rm /tmp/Predictions.$$.$IoWF.$IeWF.gff "

          IeWF=$(echo $IeWF + $deWF|bc)
          ef=`gawk "BEGIN{print ($IeWF <= $FeWF) ? 1 : 0}"`
        done

    IoWF=$(echo $IoWF + $doWF|bc)
    of=`gawk "BEGIN{print ($IoWF <= $FoWF) ? 1 : 0}"`
  done
exit
rm /tmp/tmp.$$
{ echo "###"; echo "### Execution time for [$0] : $SECONDS secs";
  echo "$L$L$L$L";
  echo ""; } 1>&2;
#
exit 0
=cut

}

sub check_augustus() {
 pod2usage
"Can't find the Augustus directory, easiest way is to create a symlink of the augustus executable somewhere in your PATH (e.g. \$HOME/bin) or add the Augustus bin directory in your PATH\n"
   unless $augustus_dir && -d $augustus_dir;
 $gff2gb_exec = $augustus_dir . '/scripts/gff2gbSmallDNA.pl'
   if !$gff2gb_exec;
 $augustus_train_exec = $augustus_dir . '/bin/etraining'
   if !$augustus_train_exec;
 $augustus_filterGenes_exec = $augustus_dir . '/scripts/filterGenes.pl'
   if !$augustus_filterGenes_exec;
 pod2usage "Can't find etraining from Augustus\n"
   if !-s $augustus_train_exec;
 pod2usage "Can't find gff2gbSmallDNA.pl from Augustus\n"
   if !-s $gff2gb_exec;
 pod2usage "Can't find filterGenes.pl from Augustus\n"
   if !-s $augustus_filterGenes_exec;
}

sub shuffle_fasta() {

 # shuffling is important to allow aligners process the data
 # in roughly equal time per thread
 my $fasta    = shift;
 my $orig_sep = $/;
 $/ = '>';
 open( IN, $fasta );

 my %sequence_data;
 while ( my $record = <IN> ) {
  chomp($record);
  my @data = split( "\n", $record );
  next unless $record;
  my $id = shift @data;
  $sequence_data{$id} = join( '', @data );
 }

 open( OUT, ">$fasta.shuff" );
 foreach my $id ( shuffle( keys %sequence_data ) ) {
  print OUT $/ . $id . "\n" . &wrap_text( $sequence_data{$id} );
 }
 $/ = $orig_sep;
 close OUT;
 return "$fasta.shuff";
}

sub splitfasta() {
 my @files;
 my $file2split         = shift;
 my $outdir             = shift;
 my $how_many_in_a_file = shift;
 my $pattern            = shift;
 my $shuffle            = shift;

 return unless $file2split && -s $file2split && $outdir && $how_many_in_a_file;
 return if -d $outdir;

 mkdir($outdir) unless -d $outdir;
 my $filecount;
 my $seqcount = int(0);
 my $shuffled_file = $shuffle ? &shuffle_fasta($file2split) : $file2split;
 open( FILE, $shuffled_file ) || die;
 my $orig_sep = $/;
 $/ = ">";

 while ( my $record = <FILE> ) {
  chomp($record);
  next unless $record;
  my @lines = split( "\n", $record );
  my $id = shift @lines;
  next if ( $pattern && $id !~ /$pattern/ );
  $id =~ /^(\S+)/;
  $id = $1 || die "Cannot find ID for a sequence in $file2split";
  chomp(@lines);
  my $seq = join( '', @lines );
  $seq =~ s/\s+//g;
  $seq =~ s/\*$//;    # for protein stop codons
  next unless $seq;
  $seqcount++;

  if (   !$filecount
       || $how_many_in_a_file == 1
       || $seqcount > $how_many_in_a_file )
  {
   $seqcount = int(0);
   $filecount++;
   my $outfile =
     $how_many_in_a_file == 1 ? $id : $file2split . "_" . $filecount;
   $outfile = $outdir . '/' . $outfile;
   close(OUT);
   open( OUT, ">$outfile" ) || die("Cannot open $outfile");
   push( @files, $outfile );
  }
  print OUT $/ . $id . "\n" . &wrap_text($seq);
 }
 close(FILE);
 close(OUT);
 $/ = $orig_sep;
 unlink($shuffled_file) if $shuffle;
 return \@files;
}

sub gff2hints() {
 my $gff         = shift;
 my $golden      = shift;
 my $src         = $golden ? 'GLD' : 'XNT';
 my $type_suffix = $golden ? '' : 'part';
 my $priority    = $golden ? 7 : 5;
 my $delimiter   = "\n\n";
 my $orig_sep    = $/;
 open( GFF, $gff )          || die;
 open( OUT, ">$gff.hints" ) || die;
 $/ = $delimiter;

 while ( my $record = <GFF> ) {
  my @lines = split( "\n", $record );

  #my $gene_line = $lines[0];
  my $mRNA_line = $lines[1];
  $mRNA_line =~ /ID=([^;]+)/;
  my $mRNA_id = $1;
  $mRNA_id =~ s/[^\w\.]+/_/g;

  # CDS intron
  for ( my $i = 2 ; $i < scalar(@lines) ; $i++ ) {
   my @data = split( "\t", $lines[$i] );
   next unless $data[8];
   if ( $data[2] eq 'mRNA' ) {
    $data[2] = 'genicpart';
    print OUT join( "\t", @data[ 0 .. 7 ] )
      . "\tsrc=$src;pri=$priority;grp=$mRNA_id\n";
   }
   elsif (    $data[2] eq 'exon'
           || $data[2] eq 'CDS'
           || $data[2] eq 'intron' )
   {
    $data[2] .= $type_suffix;
    print OUT join( "\t", @data[ 0 .. 7 ] )
      . "\tsrc=$src;pri=$priority;grp=$mRNA_id\n";
   }
   elsif ( $data[2] =~ /splice/ ) {
    $data[2] = $data[2] eq 'three_prime_cis_splice_site' ? 'ass' : 'dss';
    print OUT join( "\t", @data[ 0 .. 7 ] )
      . "\tsrc=$src;pri=$priority;grp=$mRNA_id\n";
   }
   elsif ( $data[2] =~ /UTR/i ) {
    $data[2] = 'UTRpart';
    print OUT join( "\t", @data[ 0 .. 7 ] )
      . "\tsrc=$src;pri=$priority;grp=$mRNA_id\n";
   }
  }
 }
 close GFF;
 close OUT;
 $/ = $orig_sep;
}

sub process_cmd {
 my ($cmd) = @_;
 print "CMD: $cmd\n" if $debug;
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  die "Error, cmd died with ret $ret\n";
 }
 return $ret;
}

sub gff3_fix_phase() {
 my $gff3_file = shift;
 open( IN, $gff3_file ) || confess( "Cannot find $gff3_file " . $! );
 my $index_file = "$gff3_file.inx";
 my $gene_obj_indexer = new Gene_obj_indexer( { "create" => $index_file } );
 my $asmbl_id_to_gene_list_href =
   &GFF3_utils::index_GFF3_gene_objs( $gff3_file, $gene_obj_indexer );
 open( OUT,  ">$gff3_file.gff3" );
 open( PEP,  ">$gff3_file.pep" );
 open( CDS,  ">$gff3_file.cds" );
 open( GENE, ">$gff3_file.gene" );

 foreach my $asmbl_id (sort keys %$asmbl_id_to_gene_list_href ) {
  my $genome_seq = $scaffold_seq_hashref->{$asmbl_id};
  if ( !$genome_seq ) {
   warn "Cannot find sequence $asmbl_id from genome\n";
   next;
  }
  my @gene_ids = @{ $asmbl_id_to_gene_list_href->{$asmbl_id} };
  foreach my $gene_id (@gene_ids) {
   my %params;
   my %preferences;
   $preferences{'sequence_ref'} = \$genome_seq;
   $params{unspliced_transcript} = 1;

   my $gene_obj_ref = $gene_obj_indexer->get_gene($gene_id);

   $gene_obj_ref->create_all_sequence_types( \$genome_seq, %params );

   foreach
     my $isoform ( $gene_obj_ref, $gene_obj_ref->get_additional_isoforms() )
   {

    $isoform->delete_isoforms();
    my $isoform_id = $isoform->{Model_feat_name};
    my @model_span = $isoform->get_CDS_span();
    next
      if ( !$isoform->get_CDS_span()
           || abs( $model_span[0] - $model_span[1] ) < 3 );

    eval { $isoform->set_CDS_phases( \$genome_seq ); };

    # get sequences
    # CDS
    my $seq = $isoform->get_CDS_sequence();
    $seq =~ s/(\S{60})/$1\n/g if $seq;
    chomp $seq if $seq;
    if ( $seq && length($seq) >= $minorf ) {
     print CDS ">$isoform_id $gene_id\n$seq\n";
    }
    else {
     next;
    }

    # proteins
    $seq = $isoform->get_protein_sequence();
    $seq =~ s/(\S{60})/$1\n/g;
    chomp $seq;
    print PEP ">$isoform_id $gene_id\n$seq\n";

    # gene
    $seq = $isoform->get_gene_sequence();
    $seq =~ s/(\S{60})/$1\n/g;
    chomp $seq;
    print GENE ">$isoform_id $gene_id\n$seq\n";

    # GFF3
    print OUT $isoform->to_GFF3_format_extended(%preferences) . "\n";

   }
  }
 }
 close OUT;
 close PEP;
 close CDS;
 close GENE;
 rename( $gff3_file, "$gff3_file.original" );
 unlink $index_file;
 rename( "$gff3_file.gff3", $gff3_file );
 &sort_gff3("$gff3_file");
 return;
}

sub gff_to_gtf() {

 my $gff3_file = shift;
 open( OUT, ">$gff3_file.gtf" );
 my $inx_file = "$gff3_file.inx";
 my $gene_obj_indexer = new Gene_obj_indexer( { "create" => $inx_file } );

 my $asmbl_id_to_gene_list_href =
   &GFF3_utils::index_GFF3_gene_objs( $gff3_file, $gene_obj_indexer );

 foreach my $asmbl_id ( sort keys %$asmbl_id_to_gene_list_href ) {

  ## get the genome sequence
  my $genome_seq = $scaffold_seq_hashref->{$asmbl_id};
  if ( !$genome_seq ) {
   warn "Cannot find sequence $asmbl_id\n";
   next;
  }
  my @gene_ids = @{ $asmbl_id_to_gene_list_href->{$asmbl_id} };
  foreach my $gene_id (@gene_ids) {
   ## note models of isoforms are bundled into the same gene object.
   my $gene_obj_ref = $gene_obj_indexer->get_gene($gene_id);
   foreach
     my $gene_obj ( $gene_obj_ref, $gene_obj_ref->get_additional_isoforms() )
   {
    $gene_obj->delete_isoforms();    # unbundle the model object!
    my $gtf_text = "";
    eval { $gtf_text = $gene_obj->to_GTF_format( \$genome_seq ); };
    if ($@) {

     # do it in pseudogene mode - if not then UTR is printed as CDS... not good!
     $gene_obj->{is_pseudogene} = 1;
     $gtf_text = $gene_obj->to_GTF_format( \$genome_seq );
    }
    print OUT "$gtf_text\n";
   }
  }
 }
 close OUT;
 unlink $inx_file;
}

sub check_program() {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  die "Error, path to required $prog cannot be found\n"
    unless $path =~ /^\//;
  chomp($path);
  $path = readlink($path) if -l $path;
  push( @paths, $path );
 }
 return @paths;
}

sub check_program_optional() {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  warn
"Warning: path to optional $prog cannot be found in your path environment.\n"
    unless $path =~ /^\//;
  chomp($path);
  $path = readlink($path) if -l $path;
  push( @paths, $path );
 }
 return @paths;
}

sub get_id_seq_from_fasta {
 my ( $acc, $fasta_db ) = @_;
 my $seq;
 if ( !$get_id_seq_from_fasta_hash{$fasta_db}{$acc} ) {
  $seq = `$cdbyank_exec $fasta_db.cidx -a '$acc'`;
  chomp($seq);
  $get_id_seq_from_fasta_hash{$fasta_db}{$acc} = $seq;
 }
 else {
  $seq = $get_id_seq_from_fasta_hash{$fasta_db}{$acc};
 }
 unless ($seq) {

  #warn "WARNING: couldn't retrieve seq for $acc from $fasta_db\n";
  return;
 }
 my @x = split( /\n/, $seq );
 my $id = shift @x;
 $seq = join( "", @x );
 $seq = uc($seq);
 $seq =~ s/\s//g;
 return $seq;
}

sub parse_genome_gff() {
 my %hash;
 my ( $skipped, %allowed_data );
 my $gff_file = shift;
 die unless $gff_file && -s $gff_file;
 my $program = shift;
 $program = 'PASA' unless $program;
 my $number_overlapping =
   &remove_overlapping_gff( $gff_file, $gff_file . ".nr" );
 $gff_file .= '.nr';
 open( IN, $gff_file ) || die("Cannot find $gff_file\n");

 while ( my $ln = <IN> ) {
  next if $ln =~ /^#/ || $ln =~ /^\s*$/ || $ln !~ /\tmRNA\t/;
  chomp($ln);
  my @data = split( "\t", $ln );
  next unless $data[8];
  my $mrna_id;
  if ( $data[8] =~ /ID=([^;]+)/ ) {
   $mrna_id = $1;
   $hash{ $data[0] }{$mrna_id} = $ln;
  }
  if ($no_single_exon) {
   my $exon_counting = int(0);
   while ( my $ln2 = <IN> ) {
    last if $ln2 =~ /^#/ || $ln2 =~ /^\s*$/;
    $exon_counting++ if $ln2 =~ /\texon\t/;
   }
   if ( $exon_counting <= 1 ) {
    delete( $hash{ $data[0] }{$mrna_id} );
    $skipped++;
   }
  }
 }
 close IN;

 foreach my $scaffold ( keys %hash ) {
  open( OUT, ">$genome_dir/$scaffold.gff.filter" );
  print OUT "\t$scaffold\tProduced by $gff_file and $program\n";
  foreach my $mrna_id ( keys %{ $hash{$scaffold} } ) {
   my @data = split( "\t", $hash{$scaffold}{$mrna_id} );
   my $orient = $data[6] eq '+' ? 1 : 0;
   print OUT $data[3] . "\t"
     . $data[4]
     . "\t999\t0\t0\t$orient\t0\t0\t$mrna_id\n";
   $allowed_data{$mrna_id} = 1;
  }
  close OUT;
 }
 print "Skipped $skipped single exon genes\n" if $skipped;

 # update to non-overlapping file in case we need to use it again.
 $pasa_genome_gff = $gff_file;
 return \%allowed_data;
}

sub prepare_pasa_output() {

# the point of this exercise is to get the cDNA from pasa (not just the ORF) and an
# exonerate compatible annotations file.

 my $fasta_contigs = "$pasa_gff.contigs";
 return $fasta_contigs if -s $fasta_contigs;
 print "Processing pasa output to produce exonerate-compatible input\n";

 # will sort and remove overlapping genes. returns valid genes
 my $allowed_data_hashref = &parse_genome_gff( $pasa_genome_gff, 'PASA' );

 print "Indexing...\n";
 my $contig_seq_hashref = &read_fasta($pasa_assembly_file);

 my $index_file = "$pasa_gff.inx";
 my $gene_obj_indexer = new Gene_obj_indexer( { "create" => $index_file } );
 my $asmbl_id_to_gene_list_href =
   &GFF3_utils::index_GFF3_gene_objs( $pasa_gff, $gene_obj_indexer );

 print "Finding contigs from pasa...\n";
 my ( $contig_counter, $gene_counter, $mrna_counter ) =
   ( int(0), int(0), int(0) );
 my %contigs_used;
 open( ANNOT, ">$fasta_contigs.annotations" );
 open( TROUT, ">$fasta_contigs" );

 foreach my $asmbl_id ( keys %$asmbl_id_to_gene_list_href ) {
  my @gene_ids = @{ $asmbl_id_to_gene_list_href->{$asmbl_id} };
  $contig_counter++;
  foreach my $gene_id (@gene_ids) {
   my %params;
   my $gene_obj_ref = $gene_obj_indexer->get_gene($gene_id);
   $gene_obj_ref->create_all_sequence_types( \$contig_seq_hashref->{$asmbl_id},
                                             %params );

   $gene_counter++;
   foreach
     my $isoform ( $gene_obj_ref, $gene_obj_ref->get_additional_isoforms() )
   {
    my $com_name   = $isoform->{com_name};
    my $isoform_id = $isoform->{Model_feat_name};
    next unless $allowed_data_hashref->{$isoform_id};

    if ($only_complete) {
     next
       if $isoform->{is_5prime_partial}
        || $isoform->{is_3prime_partial}
        || $isoform->{is_pseudogene}
        || ( $com_name && $com_name =~ /prime_partial/ );
    }

    #$mrna_seq is now always given in + orientation
    #my $mrna_seq = $isoform->get_cDNA_sequence();
    # my $orientation = '+';
    my $orientation = $isoform->get_orientation();
    my $mrna_seq    = $contig_seq_hashref->{$asmbl_id};

    # set it to contig as there is no other way to get co-ords currently.

    my ( $model_lend, $model_rend ) =
      sort { $a <=> $b } $isoform->get_model_span();

    # die "Already exists: $isoform_id " if $contigs_used{$isoform_id};
    my $orf_seq_length = int( abs( $model_rend - $model_lend ) + 1 );
    next if $orf_seq_length < $minorf;

    if ( $orientation eq '+' ) {
     print ANNOT join(
                       ' ',
                       (
                         $isoform_id, $orientation,
                         $model_lend, $orf_seq_length
                       )
     ) . "\n";
    }
    else {
     print ANNOT join(
      ' ',
      (
       $isoform_id,
       $orientation,

       #the co-ordinates are revcomplemented
       ( length($mrna_seq) - $model_rend + 1 ),
       $orf_seq_length
      )
     ) . "\n";
    }
    print TROUT ">$isoform_id\n" . &wrap_text($mrna_seq);
    $contigs_used{$isoform_id} = 1;

   }

  }

 }

 close ANNOT;
 close TROUT;
 die "Could not create $fasta_contigs\n"
   unless -s $fasta_contigs && -s $fasta_contigs . '.annotations';
 print "Prepared $fasta_contigs with "
   . scalar( keys %contigs_used )
   . " sequences larger than $minorf b.p. from $mrna_counter transcripts,$gene_counter genes and $contig_counter contigs\n";
 return $fasta_contigs;
}

sub run_exonerate() {
 my $is_blast;
 my ($parafly_exec) = &check_program('ParaFly');
 $exonerate_file =
     $pasa_gff
   ? $pasa_gff . '.exonerate.results'
   : $peptide_file . '.exonerate.results';

 my $exonerate_command_file = "run_exonerate_commands.cmd";

 # check if it already has been processed
 $no_rerun_exonerate = 1
   if ( -s $exonerate_file
  && -s $exonerate_command_file
  && ( -s $exonerate_command_file == -s $exonerate_command_file . '.completed' )
   );
 unless ($no_rerun_exonerate) {
  &run_aat() if !$same_species;
  if (
   -s $exonerate_command_file
   && ( !-s $exonerate_command_file . '.completed'
    || (
     -s $exonerate_command_file != -s $exonerate_command_file . '.completed' ) )
    )
  {
   print "Re-processing with exonerate\n";
   &process_cmd(
"$parafly_exec -shuffle -CPU $threads -c $exonerate_command_file -failed_cmds $exonerate_command_file.failed -v "
   );
  }
  else {

   #TODO mrna_file only?
   if ($peptide_file) {

    if ($same_species) {
     $is_blast = 1;
     &run_blast( $peptide_file, 'protein' );
    }
    else {
     &run_aat( $peptide_file, 'protein' );
    }
    print "Running exonerate...\n";
    unlink($exonerate_file);
    my $exonerate_options =
" -minorf $minorf -protein -in $peptide_file -separate -filter $genome_dir/*filter -threads $threads -intron_max $intron_size  $same_species ";
    $exonerate_options .= " -softmask -ref $softmasked_genome "
      if ($softmasked_genome);
    $exonerate_options .= " -ref $genome_file "
      if ( !$softmasked_genome );
    $exonerate_options .= " -norefine "   if $norefine;
    $exonerate_options .= " -from_blast " if $is_blast;

    &process_cmd( 'run_exonerate.pl' . $exonerate_options );
    die "Exonerate run failed for some reason....\n"
      if ( !-d $peptide_file . "_queries" );
    system(   "cat "
            . $peptide_file
            . "_queries/*exonerate_results >> $exonerate_file" );
   }
   elsif ($mrna_file) {
    die "Not implemented yet - see GMAP output\n";
   }
   else {

    # we have pasa data
    my $fasta_contigs = &prepare_pasa_output();
    if ( !$same_species ) {
     &run_aat( $fasta_contigs, 'nucl' );
    }
    print "Running exonerate...\n";
    my $exonerate_options =
" -minorf $minorf -annotation $fasta_contigs.annotations -in $fasta_contigs -separate -filter $genome_dir/*filter "
      . " -threads $threads -intron_max $intron_size $same_species ";
    $exonerate_options .= " -softmask -ref $softmasked_genome "
      if ($softmasked_genome);
    $exonerate_options .= " -ref $genome_file "
      if ( !$softmasked_genome );
    $exonerate_options .= " -norefine "   if $norefine;
    $exonerate_options .= " -from_blast " if $is_blast;
    &process_cmd( 'run_exonerate.pl' . $exonerate_options );
    die "Exonerate run failed for some reason....\n"
      if ( !-d $fasta_contigs . "_queries" );
    unlink($exonerate_file);
    system(   "cat "
            . $fasta_contigs
            . "_queries/*exonerate_results >> $exonerate_file" );

   }

  }
  $no_rerun_exonerate = 1
    if ( -s $exonerate_file
   && -s $exonerate_command_file
   && -s $exonerate_command_file . '.completed'
   && (
    -s $exonerate_command_file == -s $exonerate_command_file . '.completed' ) );

 }
 die
"Have not completed with the current exonerate run. You can force this message to be ignored with -norerun \n"
   unless $no_rerun_exonerate;
 print "Completed processing exonerate file as $exonerate_file\n";

 &correct_exonerate_gff($exonerate_file)
   unless -s $exonerate_file . '.corrected.gff3'
    && $exonerate_file . '.corrected.passed';
 die "Cannot find $exonerate_file.corrected.gff3\n"
   unless -s $exonerate_file . '.corrected.gff3';
 if ($stop_after_correction) {
  print "User asked to stop after exonerate correction\n";
  exit();
 }

 return ( $exonerate_file . '.corrected.golden',
          $exonerate_file . '.corrected.passed' );
}

sub wrap_text() {
 my $string      = shift;
 my $wrap_length = shift;
 $wrap_length = 120 if !$wrap_length;
 $string =~ s/(.{0,$wrap_length})/$1\n/g;
 $string =~ s/\n{2,}/\n/;
 return $string;
}

sub count_seq() {
 my $fasta = shift;
 my $count = int(0);
 open( IN, $fasta );
 while ( my $ln = <IN> ) {
  next unless $ln =~ /^>/;
  $count++;
 }
 close IN;
 return $count;
}

sub run_aligner() {
 my $fasta_in = shift;
 my $aligner  = shift;
 my $type     = shift;
 print "$aligner: Preparing input sequences\n";
 my $seq_count = &count_seq($fasta_in);
 my $splits    = ceil( $seq_count / $threads );

 my $output_directory = "$fasta_in.$aligner";

 if ( -d $output_directory ) {
  warn
"$fasta_in.$aligner already exists. Will NOT overwrite and will skip existing output. Stop, delete it and restart otherwise\n";

  #sleep(3);
 }
 else {
  print "$aligner: Running alignment\n";
  if (    ( $pasa_cds && $fasta_in eq $pasa_cds )
       || ( $pasa_peptides && $fasta_in eq $pasa_peptides ) )
  {
   &splitfasta( $fasta_in, $output_directory, $splits, 'type:complete', 1 );
  }
  else {
   &splitfasta( $fasta_in, $output_directory, $splits, undef, 1 );
  }
 }
 my @fastas = glob("$output_directory/*");
 return unless @fastas && scalar(@fastas) > 0;

 my $thread_helper = new Thread_helper($threads);

 foreach my $fasta ( sort @fastas ) {
  next unless -f $fasta && -s $fasta;
  next unless $fasta =~ /\d$/;
  my $thread = threads->create( 'do_' . $aligner . '_cmd', $fasta,
                                $output_directory,         $type );
  $thread_helper->add_thread($thread);
 }

 $thread_helper->wait_for_all_threads_to_complete();

 my @failed_threads = $thread_helper->get_failed_threads();
 if (@failed_threads) {
  die "Error, " . scalar(@failed_threads) . " threads failed.\n";
  exit(1);
 }
 system("find $output_directory -empty -delete");
 my @aligner_out = glob("$output_directory/*");
 die "$aligner failed" unless @aligner_out && scalar(@aligner_out) > 0;
 print "$aligner: Completed\n";
 return \@aligner_out;
}

sub do_blast_cmd() {
 my $fasta            = shift;
 my $output_directory = shift;    # currently ignored.
 my $type             = shift;
 my $blast_opt = $type eq 'protein' ? $tblastn_exec : $tblastx_exec;
 $blast_opt .=
" -num_threads 1 -evalue 1e-20 -max_target_seqs 1 -outfmt 6 -db $genome_sequence_file";
 $blast_opt .= " -max_intron_length $intron_size" if $type eq 'protein';
 $blast_opt .= " -lcase_masking" if $softmasked_genome;
 &process_cmd("$blast_opt -query $fasta -out $fasta.blast");
 unlink($fasta);
}

sub do_gmap_cmd() {
 my $fasta            = shift;
 my $output_directory = shift;    # currently ignored.
 my $type             = shift;    # currently ignored.

 my $gmap_opt =
   "$gmap_exec -D $genome_sequence_file_dir -d $genome_sequence_file_base.gmap";
 my $identical_fraction_prop = sprintf( "%.2f", $identical_fraction / 100 );
 $gmap_opt .=
" -n 0 -p 3 --nofails -B 3 -t 1 -f gff3_gene --split-output=$fasta.gmap --min-trimmed-coverage=0.95 --min-identity=$identical_fraction_prop ";
 &process_cmd("$gmap_opt $fasta 2>/dev/null") unless -s "$fasta.gmap.uniq";
 unlink($fasta);

}

sub run_gmap() {
 my $fasta       = shift;
 my $gmap_output = $fasta . '_vs_genome.gmap.gff3';

 my $output_files_ref = &run_aligner( $fasta, 'gmap', 'nuc' );
 my $orig_sep = $/;

 # our gmap output has the gene as the main ID but we need it to
 # the mrna
 print "Post-processing GMAP output...\n";
 unless ( -s $gmap_output ) {
  $/ = "###\n";
  open( OUT, ">$gmap_output" );
  foreach my $file (@$output_files_ref) {
   next unless $file =~ /uniq$/;
   open( IN, $file );
   while ( my $record = <IN> ) {
    chomp($record);
    next if !$record || $record =~ /^\s*$/;
    my @record_data = split( "\n", $record );
    my @corrected_record_data;
    for ( my $i = 0 ; $i < scalar(@record_data) ; $i++ ) {
     next if $record_data[$i] =~ /^#/;
     my @data = split( "\t", $record_data[$i] );
     $data[8] =~ s/\.path\d+//g;
     $data[8] =~ s/\.mrna\d+//;
     $data[1] = 'GMAP';
     push( @corrected_record_data, join( "\t", @data ) );
    }

    my $gene_line = shift @corrected_record_data;
    my $mrna_line = shift @corrected_record_data;

    my @gene_data = split( "\t", $gene_line );
    my @mrna_data = split( "\t", $mrna_line );
    my ( $gene_id, $mrna_id );
    if ( $gene_data[8] =~ /ID=([^;]+)/ ) {
     $gene_id = $1;
    }
    else {
     confess "Unkown GFF format for $record\n$gene_data[8]";
    }
    if ( $mrna_data[8] =~ /ID=([^;]+)/ ) {
     $mrna_id = $1;
    }
    else {
     confess "Unkown GFF format for $record\n$mrna_data[8]";
    }

    # if PASA, convert to pasa format
    if ( $gene_id =~ s/^(asmbl_\d+\|)m(\.\d+)$/$1g$2/ ) {
     $gene_line =~ s/(asmbl_\d+\|)m(\.\d+)/$1g$2/;
    }
    else {
     $gene_id .= '.gene';
     $gene_line =~ s/(ID=[^;]+)/$1.gene/;
    }
    $mrna_line =~ s/Parent=[^;]+/Parent=$gene_id/;
    print OUT $gene_line . "\n";
    print OUT $mrna_line . "\n";

    foreach my $other_line (@corrected_record_data) {
     $other_line =~ s/\.path\d+//g;
     $other_line =~ s/Parent=[^;]+/Parent=$mrna_id/;
     print OUT $other_line . "\n";
    }
    print OUT "##\n";
   }
   close IN;
  }
  close OUT;
  $/ = $orig_sep;
 }

# consider searching for overlaps after we process golden. it will
# much slower ( a lot more genes) but will not discard golden overlaps
# i tried it and found 30 more genes (i.e 0.4%) and costs 10 minutes for 100,000 mrnas

 &create_golden_gffs_from_gff( $gmap_output, $gmap_output . '.n' )
   unless ( -s "$gmap_output.passed" );
 &remove_overlapping_gff( $gmap_output . '.n', $gmap_output . '.nr' )
   unless -s $gmap_output . '.nr';
 $gmap_output .= '.nr';
 &create_golden_gffs_from_gff($gmap_output) unless ( -s "$gmap_output.passed" );

 return ( $gmap_output . '.golden', $gmap_output . ".passed" );
}

sub run_blast() {
 my $fasta        = shift;
 my $type         = shift;
 my $blast_output = $fasta . '_vs_genome.blast';

 print "Preparing for BLAST...\n";
 if ( !-s $blast_output ) {
  my $output_files_ref = &run_aligner( $fasta, 'blast', $type );
  foreach my $file (@$output_files_ref) {
   &process_cmd(
         "sort -nk9,10 $file | sort -s -k1,1 | sort -s -k2,2 >> $blast_output");
  }
 }

 #produce a hit.filter file in $genome_dir
 my %print_hash;

 open( BLAST, $blast_output ) || die("Cannot find $blast_output $!");
 while ( my $ln = <BLAST> ) {
  next if $ln =~ /^\s*$/;
  chomp($ln);
  my @data = split( "\t", $ln );
  next unless $data[10];
  next if $data[10] > 1e-20;

# AAT ( $dstart, $dend,  $score, $astart, $aend,$orient, $zero1, $zero2, $query ) = ( $1, $2, $3, $4, $5, $6, $7, $8, $9 );
# BLAST qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
  my $g_start  = $data[8] < $data[9] ? $data[8] : $data[9];
  my $g_end    = $data[8] < $data[9] ? $data[9] : $data[8];
  my $g_strand = $data[8] < $data[9] ? '1'      : '0';
  my $score =
      $data[11] > 100
    ? $data[11]
    : 101
    ; # force to keep all alignments, even the very short ones if evalue is high enough
  $print_hash{ $data[1] } .=
      "$g_start\t$g_end\t" 
    . $score . "\t"
    . $data[6] . "\t"
    . $data[7]
    . "\t$g_strand\t0\t0\t"
    . $data[0] . "\n";

 }
 foreach my $hit ( keys %print_hash ) {
  open( OUT, ">$genome_dir/$hit.blast.filter" )
    || die("Cannot write to $genome_dir/$hit.blast.filter $!");
  print OUT "\t$hit\tProduced by BLAST\n";
  print OUT $print_hash{$hit};
  close OUT;
 }
 close BLAST;
 return $blast_output;
}

sub run_aat() {
 my $input_fasta = shift;
 my $type        = shift;

# this is different as it's parallelized versus genome.
# also it is not using the pasa_cds because that is not appropriate input for exonerate.
 my ( $aat_dir, $parafly_exec ) = &check_program( 'AAT.pl', 'ParaFly' );
 $aat_dir = dirname($aat_dir);
 my $aat_command_file = "./" . basename($genome_dir) . ".commands";

 # check if it already has been processed
 return
   if ( -s $aat_command_file
        && ( -s $aat_command_file == -s $aat_command_file . '.completed' ) );

 if ( -s $aat_command_file && !-s $aat_command_file . '.completed'
      || ( -s $aat_command_file != -s $aat_command_file . '.completed' ) )
 {
  print "Re-processing with AAT\n";
  &process_cmd(
"$parafly_exec -shuffle -CPU $threads -c $aat_command_file -failed_cmds $aat_command_file.failed -v "
  );
 }
 else {
  return unless $input_fasta && $type;
  my @commands;
  if ( $type eq 'protein' ) {
   my $aat_score = $same_species ? 400 : 80;
   my $aat_word =
     $same_species
     ? 5
     : 4;    # five is much faster than default and 3 is much slower.
   print "Preparing/running AAT...\n";
   my $matrix_file = "$aat_dir/matrices/BS";
   die "Cannot find AAT's matrices/BS as $matrix_file\n"
     unless -s $matrix_file;
   foreach my $genome_file (@$uppercase_genome_files) {
    next if $genome_file =~ /\.aat\./ || -d $genome_file;
    push( @commands,
"$aat_dir/dps $genome_file $input_fasta $matrix_file -c 300000 -f $aat_score -w $aat_word -i 30 -a $intron_size > $genome_file.aat.d ;"
       . "$aat_dir/ext $genome_file.aat.d -f $aat_score > $genome_file.aat.ext ;"
       . "$aat_dir/extCollapse.pl $genome_file.aat.ext > $genome_file.aat.extCol ;"
       . "$aat_dir/filter $genome_file.aat.extCol -c 1 > $genome_file.aat.filter ; rm -f $genome_file.aat.d $genome_file.aat.ext $genome_file.aat.extCol \n"
    ) unless -s "$genome_file.aat.filter";
   }
   open( CMD, ">$aat_command_file" );
   print CMD shuffle(@commands);
   close CMD;
   &process_cmd(
"$parafly_exec -shuffle -CPU $threads -c $aat_command_file -v -failed_cmds $aat_command_file.failed"
   );
   system("find $genome_dir -empty -delete");
  }
  else {
   my $aat_score = $same_species ? 200 : 80;
   my $aat_o     = $similar_fraction - 20;
   my $aat_p     = $identical_fraction - 20;
   print "Preparing/running AAT...\n";
   foreach my $genome_file (@$uppercase_genome_files) {
    next if $genome_file =~ /\.aat\./ || -d $genome_file;

    push( @commands,
"$aat_dir/dds $genome_file $input_fasta -o $aat_o -p $aat_p -c 300000 -f $aat_score -i 30 -a $intron_size > $genome_file.aat.d ;"
       . "$aat_dir/ext $genome_file.aat.d -f $aat_score > $genome_file.aat.ext ;"
       . "$aat_dir/extCollapse.pl $genome_file.aat.ext > $genome_file.aat.extCol ;"
       . "$aat_dir/filter $genome_file.aat.extCol -c 1 > $genome_file.aat.filter ; rm -f $genome_file.aat.d $genome_file.aat.ext $genome_file.aat.extCol \n"
    ) unless -s "$genome_file.aat.filter";
   }
   open( CMD, ">$aat_command_file" );
   print CMD shuffle(@commands);
   close CMD;
   &process_cmd(
"$parafly_exec -shuffle -CPU $threads -c $aat_command_file -v -failed_cmds $aat_command_file.failed"
   );
  }
  system("find $genome_dir -empty -delete");

 }

}

sub create_golden_gffs_from_gff() {
 my $gff_file = shift;
 my $output   = shift;
 $output = $gff_file . '.golden' if !$output;
 my $delimiter = shift;
 my $orig_sep  = $/;

 &gff3_fix_phase($gff_file);

 $delimiter = &get_gff_delimiter($gff_file) if !$delimiter;
 print "Finding golden subset...\n";

 #TODO :
 # if thomas changes the gmap format:
 # $mismatch_cutoff
 #DONE:
 # $identical_fraction
 # mrna: coverage=100.0;identity=100.0
 # pep: M and *
 # genome my $number_of_ns = $scaffold_subseq =~ tr/N/N/; (40%)
 # splice sites
 # acceptor_site=AG
 # donor_site=GC
 # donor_site=GT

 my ( %tocheck, %accepted );
 my ( $counter, $failed_cutoff ) = ( int(0), int(0), int(0) );
 my $gene_seqs_hashref = &read_fasta( $gff_file . '.gene' );
 my $pep_seqs_hashref  = &read_fasta( $gff_file . '.pep' );

 $/ = $delimiter;

 open( LOG, ">$output.log" );
 open( OUT, ">$output" );
 open( GFF, $gff_file );

GENE: while ( my $record = <GFF> ) {
  chomp($record);
  next if !$record || $record =~ /^\s*$/;
  my @record_data = split( "\n", $record );
  my $gene_line   = shift @record_data;
  my $mrna_line   = shift @record_data;
  my @gene_data   = split( "\t", $gene_line );
  my @mrna_data   = split( "\t", $mrna_line );
  my ( $gene_id, $mrna_id );

  if ( $gene_data[8] =~ /ID=([^;]+)/ ) {
   $gene_id = $1;
  }
  else {
   confess "Unkown GFF format for $record\n$gene_data[8]";
  }
  if ( $mrna_data[8] =~ /ID=([^;]+)/ ) {
   $mrna_id = $1;
  }
  else {
   confess "Unkown GFF format for $record\n$mrna_data[8]";
  }
  $counter++;

  if ( $mrna_line =~ /coverage=([\d\.]+)/ ) {
   my $coverage = $1;
   if ( $coverage < 90 ) {
    print LOG "$mrna_id rejected: coverage ($coverage) below 90%.\n";
    $failed_cutoff++;
    next;
   }
  }
  if ( $mrna_line =~ /identity=([\d\.]+)/ ) {
   my $identity = $1;
   if ( $identity < $identical_fraction ) {
    print LOG
      "$mrna_id rejected: identity ($identity) below $identical_fraction.\n";
    $failed_cutoff++;
    next;
   }
  }

  my $gene_seq = $gene_seqs_hashref->{$mrna_id};
  my $pep_seq  = $pep_seqs_hashref->{$mrna_id};
  unless ( $gene_seq && $pep_seq ) {
   print LOG
"$mrna_id rejected: cannot find protein or genome sequence. Possible CDS is unavailable.\n";
   $failed_cutoff++;
   next;
  }

  if ( $pep_seq !~ /^M/ ) {
   print LOG "$mrna_id rejected: query does not start with M.\n";
   $failed_cutoff++;
   next;
  }
  if ( $pep_seq !~ /\*$/ ) {
   print LOG "$mrna_id rejected: query does not start with M.\n";
   $failed_cutoff++;
   next;
  }

  my $number_of_ns = $gene_seq =~ tr/N/N/;
  if ( $number_of_ns > ( length($gene_seq) * 0.40 ) ) {
   $failed_cutoff++;
   print LOG
     "$mrna_id rejected: genome sequence is more than 40% Ns ($number_of_ns)\n";
   next;
  }

  foreach my $other_lines (@record_data) {
   next unless $other_lines =~ /splice/;
   if ( $other_lines =~ /acceptor_site=([A-Z]+)/ ) {
    my $site = $1;
    unless ( $site eq 'AG' ) {
     print LOG
       "$mrna_id rejected: non-canonical acceptor splice site ($site)\n";
     $failed_cutoff++;
     next GENE;
    }
   }
   elsif ( $other_lines =~ /donor_site=([A-Z]+)/ ) {
    my $site = $1;
    unless ( $site eq 'GC' || $site eq 'GT' ) {
     print LOG "$mrna_id rejected: non-canonical donor splice site ($site)\n";
     $failed_cutoff++;
     next GENE;
    }
   }
  }
  $accepted{$mrna_id} = $gene_data[1];
  print OUT $record . $delimiter;
 }
 close OUT;
 close GFF;
 $/ = $orig_sep;

 open( PASS, ">$gff_file.passed" );
 my @shuffled_genes = shuffle( keys %accepted );
 foreach my $gene (@shuffled_genes) {
  print PASS $gene . "\t" . $accepted{$gene} . "\n";
 }
 close PASS;

 my $number_of_passing_genes = scalar(@shuffled_genes);
 print LOG
"Processed $counter transcripts.\nFound $number_of_passing_genes sequences passing criteria.\n\n";
 close LOG;

 print
"Processed $counter transcripts.\nFound $number_of_passing_genes sequences passing criteria.\n\n";

 #cleanup
 unlink( $gff_file . '.pep' );
 unlink( $gff_file . '.gene' );
 unlink( $gff_file . '.cds' );

}

sub create_geneid_from_gff() {
 my $gff_file  = shift;
 my $delimiter = shift;
 my $out       = "$gff_file.geneid.gff3";
 my $orig_sep  = $/;
 $delimiter = &get_gff_delimiter($gff_file) if !$delimiter;

 $/ = $delimiter;
 open( IN, $gff_file ) || die( "Cannot open $gff_file " . $! );
 open( GENEIDGFF, ">$out" );

 while ( my $record = <IN> ) {
  chomp($record);
  next unless $record;
  next if $record =~ /^\s*$/;
  my @record_data = split( "\n", $record );

  # assume gene and then mrna - always.
  my $gene_line = shift(@record_data);
  my $mrna_line = shift(@record_data);
  my @gene_data = split( "\t", $gene_line );
  my @mrna_data = split( "\t", $mrna_line );
  die "Expected GFF $gff_file to have 'gene' instead of " . $gene_data[2]
    unless $gene_data[2] eq 'gene';
  die "Expected GFF $gff_file to have 'mRNA' instead of " . $mrna_data[2]
    unless $mrna_data[2] eq 'mRNA';

  my $gene_start = $gene_data[3];
  my $gene_end   = $gene_data[4];
  my $strand     = $mrna_data[6];
  my ( $mRNA_id, $gene_id, @geneid_data );
  if ( $mrna_data[8] =~ /ID=([^;]+)/ ) {
   $mRNA_id = $1;
  }
  if ( $gene_data[8] =~ /ID=([^;]+)/ ) {
   $gene_id = $1;
   $gene_id =~ s/\.mRNA$//;
  }

  my $geneid_mRNA_id   = $mRNA_id ? $mRNA_id : $gene_id;
  my $exon_counter     = int(0);
  my $max_exon_counter = int(0);
  my $max_cds_counter  = int(0);
  foreach my $other_features_line (@record_data) {
   my @data = split( "\t", $other_features_line );
   $max_exon_counter++ if $data[2] eq 'exon';
   $max_cds_counter++  if $data[2] eq 'CDS';
  }

  # no exons or non-coding
  next if $max_exon_counter == 0;
  next if $max_cds_counter == 0;

  foreach my $other_features_line (@record_data) {
   my @data = split( "\t", $other_features_line );

   # geneid asks for exons, not necessarely coding. right?
   # but ZFF asks for exons and checks they are all coding?! argh....
   if ( $data[2] eq 'exon' ) {
    my $exon_type;
    $exon_counter++;
    if (
         $max_exon_counter == 1 && (    $data[3] == $gene_start
                                     && $data[4] == $gene_end )
         || (    $data[4] == $gene_start
              && $data[3] == $gene_end )
      )
    {
     $exon_type = "Single";
    }
    elsif ( $exon_counter == 1 ) {
     $exon_type = "First";
    }
    elsif ( $data[6] eq '-' && $data[3] == $gene_start ) {
     $exon_type = 'Terminal';
    }
    elsif ( $data[6] eq '+' && $data[4] == $gene_end ) {
     $exon_type = 'Terminal';
    }
    else {
     $exon_type = 'Internal';
    }
    my $geneid_gffstr =
        $data[0] . "\t"
      . $data[1] . "\t"
      . $exon_type . "\t"
      . $data[3] . "\t"
      . $data[4] . "\t.\t"
      . $data[6]
      . "\t.\t$geneid_mRNA_id\n";
    $geneid_data[ ( $exon_counter - 1 ) ] = $geneid_gffstr
      if $geneid_gffstr;
   }
  }

  if (@geneid_data) {
   if ( $strand eq '+' ) {
    print GENEIDGFF join( '', @geneid_data );
   }
   else {
    print GENEIDGFF join( '', reverse(@geneid_data) );
   }
   print GENEIDGFF '#$' . "\n";
  }

 }

 close GENEIDGFF;

 $/ = $orig_sep;
 return $out;
}

sub predict_orfs() {

 # will not use PFAM
 $ENV{PATH} .= ":$RealBin/../3rd_party/transdecoder";
 my $fasta               = $mrna_file;
 my ($transdecoder_exec) = &check_program($fasta);
 my $cmd                 = $transdecoder_exec
   . " -t $fasta --workdir $fasta.transdecoder --CPU $threads ";
 &process_cmd($cmd) unless -s "$fasta.transdecoder.pep";
 die "Failed to produce $fasta.transdecoder.pep\n"
   unless -s "$fasta.transdecoder.pep";
 return "$fasta.transdecoder.pep";

}

sub read_fasta() {
 my $fasta = shift;
 my %hash;
 my $orig_sep = $/;
 $/ = '>';
 open( IN, $fasta ) || confess( "Cannot open $fasta : " . $! );
 while ( my $record = <IN> ) {
  chomp($record);
  next unless $record;
  my @lines = split( "\n", $record );
  my $id    = shift(@lines);
  my $seq   = join( '', @lines );
  $seq =~ s/\s+//g;
  if ( $id && $seq && $id =~ /^(\S+)/ ) {
   $hash{$1} = $seq;
  }
 }
 close IN;
 $/ = $orig_sep;
 return \%hash;
}

sub check_for_options() {
 pod2usage if $show_help;
 pod2usage "No genome found\n" unless ( $genome_file && -s $genome_file );
 pod2usage "Provide at least one set of the input files\n" unless ($exonerate_file || ($pasa_gff && $pasa_assembly_file && $pasa_peptides && $pasa_cds && $pasa_genome_gff) || $peptide_file || $mrna_file);
 
   unless (
            ( $exonerate_file && -s $exonerate_file )
            || (    $pasa_gff
                 && -s $pasa_gff
                 && $pasa_assembly_file
                 && -s $pasa_assembly_file
                 && $pasa_peptides
                 && -s $pasa_peptides
                 && $pasa_cds
                 && -s $pasa_cds
                 && $pasa_genome_gff
                 && -s $pasa_genome_gff )
            || ( $peptide_file && -s $peptide_file )
            || ( $mrna_file    && -s $mrna_file )
   ){
    warn "A required input file is missing\n";
    warn "\t$exonerate_file not found\n" if ($exonerate_file && !-s $exonerate_file);
    warn "\t$pasa_gff not found\n" if ($pasa_gff && !-s $pasa_gff);
    warn "\t$pasa_assembly_file not found\n" if ($pasa_assembly_file && !-s $pasa_assembly_file);
    warn "\t$pasa_peptides not found\n" if ($pasa_peptides && !-s $pasa_peptides);
    warn "\t$pasa_cds not found\n" if ($pasa_cds && !-s $pasa_cds);
    warn "\t$pasa_genome_gff not found\n" if ($pasa_gff && !-s $pasa_genome_gff);
    warn "\t$peptide_file not found\n" if ($peptide_file && !-s $peptide_file);
    warn "\t$mrna_file not found\n" if ($mrna_file && !-s $mrna_file);
    die "\n";
    
   }

 die "Max intron size (-intron) cannot be 0!\n"
   unless $intron_size && $intron_size > 0;
 die "CPUs (-cpu or -thread) cannot be 0!\n" unless $threads && $threads > 0;

 die "Cannot have both PASA and -peptide\n" if $pasa_gff && $peptide_file;
 die "Cannot have both PASA and -mrna\n"    if $pasa_gff && $mrna_file;

 die "Cannot have both -mRNA and -peptide" if $mrna_file && $peptide_file;
 die "cDNA mode is only used without peptides\n" if $peptide_file && $is_cdna;

 $is_cdna = 1 if $pasa_gff;

 die "Softmasked genome file does not exist\n"
   if $softmasked_genome && !-s $softmasked_genome;

 $augustus_dir =
   $augustus_dir ? readlink($augustus_dir) : dirname( dirname($augustus_exec) )
   if $augustus_exec && !$augustus_dir;
 &check_augustus() unless ( $stop_after_correction || $stop_after_golden );

 $same_species = '-same_species' if $same_species;

}

sub get_gff_delimiter() {
 my $delimiter;
 my $file = shift;
 return if !-s $file;
 open( IN, $file ) || confess("Can't find file $file");
 my $skip = <IN>;
 while ( my $ln = <IN> ) {
  if ( $ln =~ /^\s*$/ ) {

   # note that delimiter of an empty line is actually two \n\n
   $delimiter = $ln . "\n";
   last;
  }
  elsif ( $ln =~ /^#/ ) {
   $delimiter = $ln;
   last;
  }
 }
 close IN;
 confess "I don't know what delimiter to use for $file" if !$delimiter;
 return $delimiter;
}
