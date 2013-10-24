#!/usr/bin/env perl

=pod

=head1 TODO

smth wrong here (temporarely disabled)
38G .exonerate.results.passed.genome

=head1 NAME

 prepare_golden_genes_for_predictors.pl

 Uses an existing exonerate output or just predicted proteins (e.g. from Transdecoder or just a FASTA file) to prepare gene prediction inputs for Augustus, SNAP and geneid.
 Exonerate is run (enabled via AAT) if it is not provided.
 GTF file is produced to judge quality of annotation.
 Sorts out high quality alignments from those that don't meet the criteria.
 

=head1 USAGE

Mandatory:

  * -genome :s    => Genome sequence FASTA file. Repeatmasked with Ns if softmasked genome is also provided. Otherwise not repeatmasked.

And one of:

1. If existing run with exonerate:

    -exonerate :s => Exonerate results file with special format (see perldoc)
    -cdna         => Exonerate was run with cDNA otherwise expect protein (no start/stop searched in query, just genome)

2. Otherwise:

    -transdecoder_gff :s      => GFF output from transdecoder
    -transdecoder_peptides :s => Protein output from transdecoder
    -transdecoder_assembly :s  => Contigs provided as input to transdecoder

3. OR

    -peptides :s              => A file with proteins

Other options:

    -training :i              => number of genes to go into a random training set (def. 66% of total or 5000, whichever is higher)
    -complete                 => only do full length (recommended)
    -softmasked :s            => Genome that has been softmasked for repeats
    -nosingle                 => Don't allow single exon genes
    -flanks :i                => Base pairs to add left and right when creating GenBank file for training
    -identical :i             => Identity cutoff to mark an exonerate aln as golden, out of 100 (def 95)
    -similar :i               => Similarity cutoff to mark an exonerate aln as golden, for exonerate out of 100 (def 98)
    -intron :i                => Max size of intron (def 70000 bp)
    -threads|cpu :i           => Number of threads (def 1)
    -stop_exonerate           => Stop after exonerate finishes.
    -stop_golden              => Stop after sorting out which alignments were very good. Use it to prevent checks for Augustus/snap/BioPerl etc
    -minorf :i                => Minimum size of ORF (and therefore amino acids; def. 290 bp)
    -mismatch_cutoff :i       => Minimum number of mismatches in exonerate (def 10)
    -same_species             => Contigs/proteins and genome is the same species
    -norefine                 => Don't use -refine for exonerate
    -norerun                  => Don't re-run exonerate (assume it already exists as [input file].exonerate.results

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

        1 CSIRO Ecosystem Sciences, GPO 1700, Canberra 2601, Australia
        alexie@butterflybase.org
        
=head1 DISCLAIMER & LICENSE

This software is released under the GNU General Public License version 3 (GPLv3).
It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html.
Please note that incorporating the whole software or parts of its code in proprietary software
is prohibited under the current license.

=head1 BUGS & LIMITATIONS

Please report.

If you need the UTR, then you must use the cDNA model but exonerate has bugs when writing gff with cDNA.
If you don't need the UTR, then protein alignments are faster and the GFF doesn't have issues (i think).

Exonerate bugs:
Exonerate's GFF has issues with single exon genes: CDS includes UTR (if doing cdna2genome). fix by adjusting CDS by removing UTR coords
scaffold_99:447455-451806       exonerate:cdna2genome   gene    1999    2363    1338    +       .       gene_id 0 ; sequence m.36189 ; gene_orientation .
scaffold_99:447455-451806       exonerate:cdna2genome   utr5    1999    2228    .       +       .       
scaffold_99:447455-451806       exonerate:cdna2genome   cds     1999    2363    .       +       .       
scaffold_99:447455-451806       exonerate:cdna2genome   exon    1999    2363    .       +       .       insertions 0 ; deletions 0
scaffold_99:447455-451806       exonerate:cdna2genome   similarity      1999    2363    1338    +       .       alignment_id 0 ; Query m.36189 ; Align 1999 439 230 ; Align 2229 669 135
or
scaffold_99:717941-728003       exonerate:cdna2genome   gene    2001    8063    25933   +       .       gene_id 0 ; sequence m.107900 ; gene_orientation +
scaffold_99:717941-728003       exonerate:cdna2genome   utr5    2001    2080    .       +       .       
scaffold_99:717941-728003       exonerate:cdna2genome   exon    2001    2080    .       +       .       insertions 0 ; deletions 0
scaffold_99:717941-728003       exonerate:cdna2genome   splice5 2081    2082    .       +       .       intron_id 1 ; splice_site "GT"
scaffold_99:717941-728003       exonerate:cdna2genome   intron  2081    2742    .       +       .       intron_id 1
scaffold_99:717941-728003       exonerate:cdna2genome   splice3 2741    2742    .       +       .       intron_id 0 ; splice_site "AG"
scaffold_99:717941-728003       exonerate:cdna2genome   utr5    2743    7443    .       +       .       
scaffold_99:717941-728003       exonerate:cdna2genome   cds     2743    7764    .       +       .       
scaffold_99:717941-728003       exonerate:cdna2genome   utr3b   7765    8063    .       +       .       
scaffold_99:717941-728003       exonerate:cdna2genome   exon    2743    8063    .       +       .       insertions 0 ; deletions 0
scaffold_99:717941-728003       exonerate:cdna2genome   similarity      2001    8063    25933   +       .       alignment_id 0 ; Query m.107900 ; Align 2001 1 80 ; Align 2743 81 4701 ; Align 7444 4782 321 ; Align 7765 5103 299
or
scaffold_99:779380-785309       exonerate:cdna2genome   gene    1570    3930    8824    -       .       gene_id 0 ; sequence m.108043 ; gene_orientation .
scaffold_99:779380-785309       exonerate:cdna2genome   utr5    3327    3930    .       -       .       
scaffold_99:779380-785309       exonerate:cdna2genome   cds     2442    3930    .       -       .       
scaffold_99:779380-785309       exonerate:cdna2genome   utr3b   1570    2441    .       -       .       
scaffold_99:779380-785309       exonerate:cdna2genome   exon    1570    3930    .       -       .       insertions 0 ; deletions 0
scaffold_99:779380-785309       exonerate:cdna2genome   similarity      1570    3930    8824    -       .       alignment_id 0 ; Query m.108043 ; Align 3931 2362 604 ; Align 3327 1758 885 ; Align 2442 873 872


also it causes some weird overlapping CDSs
scaffold_620:11428-19100        exonerate:cdna2genome   cds     4265    4386    .       +       .       
scaffold_620:11428-19100        exonerate:cdna2genome   exon    4265    4386    .       +       .       insertions 0 ; deletions 0
scaffold_620:11428-19100        exonerate:cdna2genome   splice5 4387    4388    .       +       .       intron_id 4 ; splice_site "GT"
scaffold_620:11428-19100        exonerate:cdna2genome   intron  4387    5038    .       +       .       intron_id 4
scaffold_620:11428-19100        exonerate:cdna2genome   splice3 5037    5038    .       +       .       intron_id 3 ; splice_site "AG"
scaffold_620:11428-19100        exonerate:cdna2genome   cds     4265    5038    .       +       .       
scaffold_620:11428-19100        exonerate:cdna2genome   utr3b   5039    5673    .       +       .       
scaffold_620:11428-19100        exonerate:cdna2genome   exon    5039    5673    .       +       .       insertions 0 ; deletions 0
basically last CDS is not supported by an exon as well as overlapping previous CDS 
funny enough, if this makes the gene single exon, the UTR bug above doesn't happen

=head1 TODO

fix bugs
splice site recognition using exonerate (splice +- 40bp); nb gc/ag is fine as a splice site.
add zff in first pass. use fathom to prune bad ones before golden set (some cds incomplete for some reason).
fathom zff fasta -validate 2>&1|grep error # parse first word: ^(\S+)
fathom zff fasta -validate 2>&1|grep warnings|egrep -v 'GC\.\.AG|short' # pass first word after colon: ^\S:\s(\S+)
5311 genes, 5245 OK, 413 warnings, 66 errors

=cut
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin";

use Digest::SHA qw/sha1_hex/;
use Data::Dumper;
use File::Basename;
use List::Util 'shuffle';

use Gene_obj;
use Gene_obj_indexer;
use CdbTools;
use GFF3_utils;
use GTF_utils;



#options
my (
     $contig_file,              $exonerate_file,  
     $genome_file,              $is_cdna,
     $no_single_exon,           $overwrite,
     $only_complete,            $peptide_file,
     $stop_after_golden,        $transdecoder_gff,
     $transdecoder_assembly_file, $transdecoder_peptides, $mrna_file,
     $softmasked_genome,        $stop_after_correction, $norefine, $nodataprint
);

my $no_rerun_exonerate;
my $same_species = '';
my $intron_size              = 70000;
my $training_set_size        = 4000;
my $aug_optimization_geneset = 450;
my $augustus_flank_region    = 4000;
my $identical_fraction       = 95;
my $minorf                   = 290;     #minimum orf size in bp
my $similar_fraction         = 98;
my $threads                  = 1;
my $mismatch_cutoff = 10;

#globals
my $failed_cutoff = int(0);
my ($cdbfasta_exec,$cdbyank_exec) = &check_program('cdbfasta','cdbyank');
my (%get_id_seq_from_fasta_hash);

GetOptions(
            'exonerate:s'             => \$exonerate_file,
            'genome:s'                => \$genome_file,
            'softmasked:s'            => \$softmasked_genome,
            'transdecoder_gff:s'      => \$transdecoder_gff,
            'transdecoder_peptides:s' => \$transdecoder_peptides,
            'transdecoder_assembly:s' => \$transdecoder_assembly_file,
            'peptides:s'              => \$peptide_file,
            'mrnas:s'                 => \$mrna_file,
            'cdna'                    => \$is_cdna,
            'training:i'              => \$training_set_size,
            'flanks:i'                => \$augustus_flank_region,
            'identical:i'             => \$identical_fraction,
            'similar:i'               => \$similar_fraction,
            'overwrite'               => \$overwrite,
            'intron:i'                => \$intron_size,
            'threads|cpu:i'           => \$threads,
            'complete'                => \$only_complete,
            'stop_exonerate'          => \$stop_after_correction,
            'stop_golden'             => \$stop_after_golden,
            'minorf:i'                  => \$minorf,
            'mismatch_cutoff:i'  => \$mismatch_cutoff,
	    'same_species' => \$same_species,
            'norefine'  => \$norefine,
            'norerun' => \$no_rerun_exonerate,
	    'nodataprint' => \$nodataprint
);
pod2usage "A required file is missing"
  unless (
           ( $exonerate_file && -s $exonerate_file )
           || (    $transdecoder_gff
                && -s $transdecoder_gff
                && $transdecoder_assembly_file
                && -s $transdecoder_assembly_file
                && $transdecoder_peptides
                && -s $transdecoder_peptides )
           || ($peptide_file && -s $peptide_file) || ($mrna_file && -s $mrna_file)
  )
  && ($genome_file
  && -s $genome_file);
die "Cannot have both TransDecoder and peptides\n" if $transdecoder_gff && $peptide_file;
die "cDNA mode is only used without peptides\n" if $peptide_file && $is_cdna;
die "Softmasked genome file does not exist\n" if $softmasked_genome && !-s $softmasked_genome;

my ($gff2gb_exec,$fathom_exec,$augustus_exec,$augustus_train_exec,$augustus_filterGenes_exec) = &check_program_optional('gff2gbSmallDNA.pl','fathom','augustus','etraining','filterGenes.pl');
my $augustus_dir= dirname(dirname($augustus_exec)) if $augustus_exec;
&check_augustus() unless ($stop_after_correction || $stop_after_golden);

$same_species = '-same_species' if  $same_species;

my $genome_sequence_file = $softmasked_genome ? $softmasked_genome : $genome_file; # the full sequence, not repeatmasked
unlink("$genome_sequence_file.cidx");
&process_cmd("$cdbfasta_exec $genome_sequence_file");

&run_exonerate();


if ($exonerate_file) {

  &correct_exonerate_gff($exonerate_file) unless -s $exonerate_file.'.corrected.gff3' && $exonerate_file.'.corrected.passed';
  die unless -s $exonerate_file.'.corrected.gff3';
  &process_for_gene_prediction($exonerate_file.'.corrected.gff3');
  print "Done!\n";
}


############################
sub correct_exonerate_gff() {
  my $exonerate_file = shift;
  print "Processing $exonerate_file using these cut-offs: identical:$identical_fraction similar:$similar_fraction;  No more than $mismatch_cutoff mismatches in total; Methionine and * in protein sequence; Start codon and stop codon in genome; No overlapping genes; Splice sites (GT..AG or GC..AG)\n";

  #pod2usage unless $exonerate_file && -s $exonerate_file;
  # IF WE NEED TO GET WHOLE SCAFFOLD SEQUENCE:
  pod2usage "Cannot find a file $exonerate_file or $genome_sequence_file\n"
    unless $exonerate_file
      && $genome_sequence_file
      && -s $exonerate_file
      && -s $genome_sequence_file;

  #  my $genome_dir = $genome_sequence_file.'_dir';
  #  &splitfasta( $genome_sequence_file, $genome_dir, 1 );
## --ryo "RYOAP_START\nRYOAP_STATS_D\tAlignment length\tidentical\tsimilar\tmismatch\tidentical%%\tsimilar%%\nRYOAP_STATS\t%s\t%et\t%ei\t%es\t%em\t%pi\t%ps\nRYOAP_CODING_QUERY\t%qi\t%qab\t%qae\n>%qi\n%qas\nRYOAP_CODING_GENOME\t%ti\t%tcb\t%tce\n>%ti\n%tcs\nRYOAP_END\n"
  open( EXONERATE, $exonerate_file );
  open( CORRECTED_EXONERATE_GFF,    ">$exonerate_file.corrected.gff3" );
  open( GENEIDGFF, ">$exonerate_file.corrected.gff3.geneid.gff3" );
  my ( %details, $flag, $scounter, $ecounter, %vulgar_data, $gene_counter, %already_printed );
###########################################################################################
# EXONERATE protein has no UTR, the gene is actually the coding part of the mRNA
# get gene data
###########################################################################################
  my $header = <EXONERATE>;
    die "It seems that the exonerate was run using protein2genome but -cdna was requested. Aborting\n" if ($header =~/protein2genome/ && $is_cdna);
    die "It seems that the exonerate was run using cdna2genome but -cdna was not requested. Aborting\n" if (($header =~/cdna2genome/ || $header =~/coding2genome/)  && !$is_cdna);
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
      # making other corrections. and fixing exonerate, geneid, snap etc bugs
###########################################################################################
      my ($gene_start, $gene_end );
      my $exon_counter   = int(0);
      my $cds_counter    = int(0);
      my $intron_counter = int(0);
      my $splice_counter = int(0);
      my ( $mRNA, $mRNA_id );
      $exonerate_data[0]=~/sequence (\S+)/;
      my $gene_id = $1 if $1;
      
 	  # if it has already been printed, skip it. someitmes there is a bug and this creates non-unique IDs
	  if ($already_printed{$gene_id}){
		warn "Warning: $gene_id was already printed before. skipping...\n";
		next;
	  }
	  $already_printed{$gene_id}=1;

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
          } else {
            $mRNA_id = $gene_id.".mRNA";
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
        } elsif ( $data[2] =~ /utr/ ) {
          $data[2] = 'UTR';
          $data[8] = "ID=$gene_id." . $data[2] . ";Parent=$mRNA_id\n";
        } elsif ( $data[2] =~ /splice/ ) {
          $splice_counter++;
          $data[8] =
            "ID=$gene_id." . $data[2] . ":$splice_counter;Parent=$mRNA_id\n";
          $data[2] = 'splice_junction';
        } elsif ( $data[2] eq 'cds' ) {
          if ( $cds_counter > 0 ) {
            ### check if the next line is an exon. if it is not, check if there is an overlapping previous CDS (exonerate cdna2genome bug)
            # if there is we have to delete $exonerate_data[$i]; and move on
            my $previous = $overlap_check[$cds_counter-1];
            my @next_data = split( "\t", $exonerate_data[ $i + 1 ] ) if $exonerate_data[ $i + 1 ];
            if ( $previous && $next_data[2] && $next_data[2] ne 'exon' ) {
                if ( $data[6] eq '+' ) {
                  # start must be higher than previous end
                  if ( $data[3] < $previous->{'end'} ) {
                    #warn $data[8].":Start is smaller than previous end:\n".$data[3]." vs ".$previous->{'end'}."\n";
                    delete( $exonerate_data[$i] );
                    next;
                  }
                } else {

                  # end must be lower than previous start
                  if ( $data[4] > $previous->{'start'} ) {
                    delete( $exonerate_data[$i] );
                    next;
                  }
                }
            }
          } elsif ( scalar(@overlap_check) == 1 ) {

# if this is the first and only CDS, then it is a single exon gene. exonerate cdna2genome bug includes UTR co-ordidates. fix it
# by getting previous UTR line and adjust co-ordidates to by utr[end]+1 if + or utr[start]-1 if -.
            my @previous_data = split( "\t", $exonerate_data[ $i - 1 ] );
            if ( @previous_data && $previous_data[2] =~ /utr/i ) {
              if ( $data[6] eq '+' ) {
                if ( $data[3] < $previous_data[4] ) {
                  $data[3] = $previous_data[4] + 1;
                }
              } else {
                if ( $data[4] > $previous_data[3] ) {
                  $data[4] = $previous_data[3] - 1;
                }
              }
            }
          }
          $cds_counter++;
          $data[2] = 'CDS';
          $data[8] = "ID=cds.$gene_id;Parent=$mRNA_id\n";
        } elsif ( $data[2] eq 'exon' ) {
          $exon_counter++;
          $data[8] = "ID=$gene_id.exon.$exon_counter;Parent=$mRNA_id\n";
        } elsif ( $data[2] eq 'intron' ) {
          $intron_counter++;
          $data[8] = "ID=$gene_id.intron.$intron_counter;Parent=$mRNA_id\n";
        }
        $data[8] =~ s/\s+\;\s+/\;/g;
        $exonerate_data[$i] = join( "\t", @data );
      }
#################################################################################
      # now printing data
################################################################################
      my ( @geneid_data, $geneid_strand );
      my $max_exon_counter   = $exon_counter;
      my $max_cds_counter    = $cds_counter;
      my $max_intron_counter = $intron_counter;
      my $max_splice_counter = $splice_counter;
      $exon_counter   = 0;
      $cds_counter    = 0;
      $intron_counter = 0;
      $splice_counter = 0;

      

      for ( my $i = 0 ; $i < scalar(@exonerate_data) ; $i++ ) {
        my $gff_line = $exonerate_data[$i] || next;
        my @data = split( "\t", $gff_line );
        if ( $data[2] eq 'gene' ) {
          print CORRECTED_EXONERATE_GFF join( "\t", @data );
          print CORRECTED_EXONERATE_GFF $mRNA if $mRNA;
          my $geneid_mRNA_id = $mRNA_id ? $mRNA_id : $gene_id;

          # save geneid strand
          $geneid_strand = $data[6];
          $gene_counter++;
        } else {
          print CORRECTED_EXONERATE_GFF join( "\t", @data );

          #geneID. glimmer
          # geneid asks for exons, not necessarely coding. right?
          # but ZFF asks for exons and checks they are all coding?! argh....
          if ( $data[2] eq 'exon' ) {
            my $exon_type;
            $exon_counter++;
            if ( $max_exon_counter == 1
                 && ( $data[3] == $gene_start && $data[4] == $gene_end )
                 || ( $data[4] == $gene_start && $data[3] == $gene_end ) )
            {
              $exon_type = "Single";
            } elsif ( $exon_counter == 1 ) {
              $exon_type = "First";
            } elsif ( $data[6] eq '-' && $data[3] == $gene_start ) {
              $exon_type = 'Terminal';
            } elsif ( $data[6] eq '+' && $data[4] == $gene_end ) {
              $exon_type = 'Terminal';
            } else {
              $exon_type = 'Internal';
            }
            my $geneid_mRNA_id = $mRNA_id ? $mRNA_id : $gene_id;
            my $geneid_gffstr =
                $data[0] . "\t"
              . "FORGENEID" . "\t"
              . $exon_type . "\t"
              . $data[3] . "\t"
              . $data[4] . "\t.\t"
              . $data[6]
              . "\t.\t$geneid_mRNA_id\n";
            $geneid_data[ ( $exon_counter - 1 ) ] = $geneid_gffstr
              if $geneid_gffstr;
          }
        }
      }

      # delimit end of gene model
      print CORRECTED_EXONERATE_GFF "##\n";

      # print geneid
      if (@geneid_data) {
        if ( $geneid_strand eq '+' ) {
          print GENEIDGFF join( '', @geneid_data );
        } else {
          print GENEIDGFF join( '', reverse(@geneid_data) );
        }
        print GENEIDGFF '#$' . "\n";
        undef(@geneid_data);
        undef($geneid_strand);
      }
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
        $query_data[2] <= $query_data[3] ? $query_data[2] : $query_data[3];
      my $query_end =
        $query_data[3] >= $query_data[2] ? $query_data[3] : $query_data[2];
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
        $target_data[2] <= $target_data[3] ? $target_data[2] : $target_data[3];
      my $target_end =
        $target_data[3] >= $target_data[2] ? $target_data[3] : $target_data[2];
      my $alignment_strand =
        $target_data[2] <= $target_data[3] ? int(1) : int(-1);
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
      $details{$query_id} = {
                              'counter'          => $scounter,
                              'qlength'          => $qlength,
                              'tlength'          => $tlength + $offset,
                              'score'            => $stats_data[1],
                              'mismatch_number'  => $stats_data[5],
                              'identical_frac'   => $stats_data[6],
                              'similar_frac'     => $stats_data[7],
                              'query_id'         => $query_id,
                              'query_seq'        => uc($query_seq),
                              'target_seq'       => uc($target_seq),
                              'query_start'      => $query_start,
                              'query_end'        => $query_end,
                              'target_id'        => $target_id,
                              'target_start'     => $target_start + $offset,
                              'target_end'       => $target_end + $offset,
                              'alignment_strand' => $alignment_strand,
      };
    }
  }
  print "Processed $gene_counter gene models.\n"; 
  close(EXONERATE);
  close(GENEIDGFF);
  close CORRECTED_EXONERATE_GFF;
  die "Report seems incomplete $scounter != $ecounter!\n"
    unless $scounter && $ecounter && ( $scounter == $ecounter );
  &sort_gff3( $exonerate_file.'.corrected.gff3', "##\n" );
  &order_fasta($genome_sequence_file,$exonerate_file.'.corrected.gff3');
  &gff3_fix_phase_GTF($exonerate_file.'.corrected.gff3') unless -s $exonerate_file.'.corrected.gff3.cds.gtf';

  if ($stop_after_correction){
    print "User asked to stop after correction\n";
    exit();
  }
  &create_golden_gffs($exonerate_file.'.corrected',\%details,\%vulgar_data) unless -s "$exonerate_file.corrected.passed";
}

sub create_golden_gffs(){
  my $exonerate_file = shift;
  my $details_hashref = shift;
  my $vulgar_data_hashref = shift;
###########################################################################################
  print "Finding golden subset...\n";
  my $query_outfile =
    $is_cdna ? "$exonerate_file.passed.cDNA" : "$exonerate_file.passed.protein";
  open( OUT,    ">$exonerate_file.data" ) unless $nodataprint;
  open( VULGAR, ">$exonerate_file.passed.vulgar" );
  open( LOG,    ">$exonerate_file.golden.log" );
  open( PASSP,  ">$query_outfile" );
#  open( PASSG,  ">$exonerate_file.passed.genome" )  unless $nodataprint;
  open( PASSGC, ">$exonerate_file.passed.genome.cds" )  unless $nodataprint;
  print OUT
    "#Hit number\tQuery length\tTarget length\tScore\tNumber of mismatches\tFraction of identical\tFraction of similar\tquery name\tquery start\tquery end\tquery sequence\tAlignment direction\ttarget name\ttarget start\ttarget end\ttarget sequence\n" 
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
    ){
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
      . $query_id . "\t"
      . $details_hashref->{$query_id}{'query_start'} . "\t"
      . $details_hashref->{$query_id}{'query_end'} . "\t"
      . $details_hashref->{$query_id}{'query_seq'} . "\t"
      . $details_hashref->{$query_id}{'alignment_strand'} . "\t"
      . $details_hashref->{$query_id}{'target_id'} . "\t"
      . $details_hashref->{$query_id}{'target_start'} . "\t"
      . $details_hashref->{$query_id}{'target_end'} . "\t"
      . $details_hashref->{$query_id}{'target_seq'} . "\n"  unless $nodataprint;
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
    unless ( $details_hashref->{$query_id}{'identical_frac'} >= $identical_fraction ) {
      print LOG
        "$query_id rejected: identical_frac less than $identical_fraction: "
        . $details_hashref->{$query_id}{'identical_frac'} . "\n";
      $failed_cutoff++;
      next;
    }
    unless ( $details_hashref->{$query_id}{'similar_frac'} >= $similar_fraction ) {
      print LOG "$query_id rejected: similar_frac less than $similar_fraction: "
        . $details_hashref->{$query_id}{'similar_frac'} . "\n";
      $failed_cutoff++;
      next;
    }
    unless ( $details_hashref->{$query_id}{'mismatch_number'} <= $mismatch_cutoff ) {
      print LOG "$query_id rejected: mismatches more than $mismatch_cutoff: "
        . $details_hashref->{$query_id}{'mismatch_number'} . "\n";
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
    unless ( substr( $details_hashref->{$query_id}{'target_seq'}, 0, 3 ) eq 'ATG' ) {
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
    my $scaffold_seq = &get_id_seq_from_fasta($details_hashref->{$query_id}{'target_id'},$genome_sequence_file);
    die "Cannot find genome entry " . $details_hashref->{$query_id}{'target_id'} . "\n" if ( !$scaffold_seq );

    #bug check
    unless (    int( $details_hashref->{$query_id}{'target_start'} )
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
                                 substr($scaffold_seq,$details_hashref->{$query_id}{'target_start'} - 2 - 1,$details_hashref->{$query_id}{'target_end'} - 1) # starts from 0
      );
    } else {
      $scaffold_subseq =
                        substr($scaffold_seq,$details_hashref->{$query_id}{'target_start'} ,$details_hashref->{$query_id}{'target_end'} +2)
    }
    if ( $no_single_exon
       && length($scaffold_subseq) <= length( $details_hashref->{$query_id}{'target_seq'} ) )
    {
      $failed_cutoff++;
      print LOG "$query_id rejected: single coding exon\n";
      next;
    }
    
    my $number_of_ns = $scaffold_subseq =~ tr/N/N/;
    if ( $number_of_ns > ( length($scaffold_subseq) * 0.40 ) ) {
      $failed_cutoff++;
      print LOG "$query_id rejected: genome sequence is more than 40% Ns\n";
      next;
    }
    print LOG "$query_id accepted\n";
    print VULGAR $vulgar_data_hashref->{$query_id} if $vulgar_data_hashref->{$query_id};
    $pass_check{ 's'
        . $details_hashref->{$query_id}{'target_start'} . 'e'
        . $details_hashref->{$query_id}{'target_end'} } = $query_id;
    print PASSP ">$query_id\n" . $details_hashref->{$query_id}{'query_seq'} . "\n";
#    print PASSG ">$query_id "     . $details_hashref->{$query_id}{'target_id'} . " genome sequence\n$scaffold_subseq\n" unless $nodataprint;
    print PASSGC ">$query_id "
      . $details_hashref->{$query_id}{'target_id'}
      . " genome CDS\n"
      . $details_hashref->{$query_id}{'target_seq'} . "\n" unless $nodataprint;
    $accepted{$query_id} = sha1_hex($scaffold_subseq);
    $pass_counter++;
  }
  print "\r$counter/$query_number    $pass_counter passed           \n";
  close OUT  unless $nodataprint;
  close PASSP unless $nodataprint;
#  close PASSG unless $nodataprint;
  close VULGAR;

  my $number_of_passing_genes = scalar( keys %accepted );
  open( PASS, ">$exonerate_file.corrected.gff3.passed" );
  my @shuffled_genes = shuffle( keys %accepted );
  foreach my $gene (@shuffled_genes) {
    print PASS $gene . "\t" . $accepted{$gene} . "\n";
  }
  close PASS;

  print LOG "Processed $counter sequences.\nFound $number_of_passing_genes sequences passing criteria.\n\n";
  close LOG;
  print "Processed $counter sequences.\nFound $number_of_passing_genes sequences passing criteria.\n\n";


  if ($stop_after_golden) {
    print "User stop requested\n";
    exit(0);
  }
}

sub process_for_gene_prediction(){
  print "Processing for gene prediction software...\n";
  my $gff_file = shift;
  my (%accepted, %training_genes, %aug_optimization_genes );
  my $accepted_hashref= \%accepted;
  my $training_genes_hashref = \%training_genes;
  my $aug_optimization_genes_hashref = \%aug_optimization_genes;

  my $passed_genes_in_training = int(0);
  my $passed_genes_in_optimization = int(0);
  my $number_of_passing_genes = int(0);

  my $pass_file = $gff_file;
  $pass_file =~s/.gff3$/.passed/;
  open( PASS, $pass_file ) || die "Cannot find pass file $pass_file";
  while (my $ln=<PASS>) {
    chomp($ln);
    my @data = split("\t",$ln);
    next unless $data[1];
    $accepted{$data[0]} = $data[1];
    $number_of_passing_genes++;
  }
  close PASS;

  my $training_set_size = int( $number_of_passing_genes * 0.40 ) if $training_set_size > int( $number_of_passing_genes * 0.40 );
  foreach my $gene (shuffle (keys %accepted)){
    $training_genes{$gene} = 1;
    $passed_genes_in_training++;
    last if $passed_genes_in_training >= $training_set_size;
  }
  foreach my $gene (shuffle (keys %accepted)) {
    next if $training_genes{$gene};
    $aug_optimization_genes{$gene} = 1;
    $passed_genes_in_optimization++;
    last if $passed_genes_in_optimization >= $aug_optimization_geneset;
  }
  $aug_optimization_geneset = scalar( keys %{$aug_optimization_genes_hashref} );
  $training_set_size        = scalar( keys %{$training_genes_hashref} );
  die "Augustus optimization gene set is zero!\n" if $aug_optimization_geneset < 1;
  die "Augustus Training gene set is zero!\n" if $training_set_size < 1;



  print "Checking for overlaps using main GFF3 file\n";
  $failed_cutoff += &filter_gff( "$gff_file",
                 $accepted_hashref, "$gff_file.golden.gff3",
                 1 );
  &gff2zff("$gff_file.golden.gff3","$gff_file.golden.zff");
  
  my @failed = `$fathom_exec $gff_file.golden.zff $gff_file.golden.gff3.fasta -validate 2>&1|grep 'error\b'`;
  if (@failed){
    print "Fathom tells us that ".scalar(@failed)." annotations failed... Taking them into account...\n";
    $failed_cutoff+=scalar(@failed);
    foreach my $fail (@failed){
      $fail=~/^(\S+)/;
      my $id = $1;
      delete($accepted_hashref->{$id});
      delete($training_genes_hashref->{$id});
      delete($aug_optimization_genes_hashref->{$id});
    }
    &filter_gff( "$gff_file",
                 $accepted_hashref, "$gff_file.golden.gff3"
                  );
    &gff2zff("$gff_file.golden.gff3","$gff_file.golden.zff");
  } 
  my @warnings = `$fathom_exec $gff_file.golden.zff $gff_file.golden.gff3.fasta -validate 2>&1|grep warnings|egrep -v 'GC\.\.AG|short'`;
  open (FATHOM,">$gff_file.golden.zff.warnings");
  print FATHOM @warnings;
  close FATHOM;
  if (scalar(@warnings)>1){
    pop(@warnings);
   print "Fathom gave us ".scalar(@warnings)." warnings (possible non-canonical splice sites, see $gff_file.golden.zff.warnings). Fixing...\n"; 
  $failed_cutoff+=scalar(@warnings);
  foreach my $fail (@warnings){
      $fail=~/^\S+:\s+(\S+)/;
      my $id = $1;
      delete($accepted_hashref->{$id});
      delete($training_genes_hashref->{$id});
      delete($aug_optimization_genes_hashref->{$id});
    }
    &filter_gff( "$gff_file",
                 $accepted_hashref, "$gff_file.golden.gff3"
                 );
    &gff2zff("$gff_file.golden.gff3","$gff_file.golden.zff");
  }

  print "Checking for overlaps using geneid file\n";
  &sort_gff3( "$gff_file.geneid.gff3", '#$' . "\n" );
  &filter_gff2_geneid(
                         "$gff_file.geneid.gff3",
                         $accepted_hashref,
                         "$gff_file.geneid.golden.gff3",
                         1
  );
  

  print "Processing for gene predictors...\n";

  print "'Golden': $number_of_passing_genes. Filtering GFFs...\n"
#    . "Training set: $training_set_size. Test: "
#    . ( $number_of_passing_genes - $training_set_size ) . ".\n"
#    . "Augustus optimization set: "
#    . $aug_optimization_geneset . ".\n"
#    . " (as the optimization set is taken from the test set,\n"
#    . " the Augustus test set is "
#    . ( $number_of_passing_genes - $training_set_size - $aug_optimization_geneset )
#    . ". more or less \n"
    . " - Filters may reject some during checking/conversion\n";
  print "Processing...\n";


  #augustus
  print "\tAugustus\n";

  #for training
  if ( -s "$gff_file" ) {
    
    &gff2hints( "$gff_file.golden.gff3", 1 )
      if -s "$gff_file.golden.gff3";
    &filter_gff( "$gff_file.golden.gff3",
                 $training_genes_hashref,
                 "$gff_file.golden.train.gff3" );
    &filter_gff(
                 "$gff_file.golden.train.gff3.rest",
                 $aug_optimization_genes_hashref,
                 "$gff_file.golden.optimization.gff3"
    );
    rename( "$gff_file.golden.optimization.gff3.rest","$gff_file.golden.test.gff3" );
    system("$gff2gb_exec $gff_file.golden.train.gff3 $genome_sequence_file $augustus_flank_region $gff_file.golden.train.gb >/dev/null"    );
    system("$augustus_train_exec --species=generic $gff_file.golden.train.gb 2>&1 | grep 'n sequence' | perl -pe 's/.*n sequence (\\S+):.*/\$1/' | sort -u > $gff_file.golden.train.gb.bad.lst"    );
    system("$augustus_filterGenes_exec $gff_file.golden.train.gb.bad.lst $gff_file.golden.train.gb > $gff_file.golden.train.good.gb"    );

    system("$gff2gb_exec $gff_file.golden.test.gff3 $genome_sequence_file $augustus_flank_region $gff_file.golden.test.gb  >/dev/null"    );
    system("$augustus_train_exec --species=generic $gff_file.golden.test.gb 2>&1 | grep 'n sequence' | perl -pe 's/.*n sequence (\\S+):.*/\$1/' | sort -u > $gff_file.golden.test.gb.bad.lst"    );
    system("$augustus_filterGenes_exec $gff_file.golden.test.gb.bad.lst $gff_file.golden.test.gb > $gff_file.golden.test.good.gb"    );

    system("$gff2gb_exec $gff_file.golden.optimization.gff3 $genome_sequence_file $augustus_flank_region $gff_file.golden.optimization.gb  >/dev/null"    );
    system("$augustus_train_exec --species=generic $gff_file.golden.optimization.gb 2>&1 | grep 'n sequence' | perl -pe 's/.*n sequence (\\S+):.*/\$1/' | sort -u > $gff_file.golden.optimization.gb.bad.lst"    );
    system("$augustus_filterGenes_exec $gff_file.golden.optimization.gb.bad.lst $gff_file.golden.optimization.gb > $gff_file.golden.optimization.good.gb"    );
  }

  #parse GB
  my $train_ref = &parse_gb("$gff_file.golden.train.good.gb");
  my $test_ref  = &parse_gb("$gff_file.golden.test.good.gb");
  my $opt_ref   = &parse_gb("$gff_file.golden.optimization.good.gb");

  #hints
  &gff2hints("$gff_file")
    if -s "$gff_file";
  &gff2hints("$gff_file.golden.gff3.rest")
    if -s "$gff_file.golden.gff3.rest";

  #snap already have made zff. $gff_file.golden.zff and $gff_file.golden.gff3.fasta
  print "\tsnap\n";
#  &gff2zff("$gff_file.golden.train.gff3","$gff_file.golden.train.zff");
#  &gff2zff("$gff_file.golden.test.gff3","$gff_file.golden.test.zff");
#  &gff2zff("$gff_file.golden.optimization.gff3","$gff_file.golden.optimization.zff");
#    &order_fasta( $genome_sequence_file, "$gff_file.golden.optimization.gff3" );
#    &order_fasta( $genome_sequence_file, "$gff_file.golden.train.gff3" );
#    &order_fasta( $genome_sequence_file, "$gff_file.golden.test.gff3" );


  #geneid
  print "\tgeneid\n";
  if ( -s "$gff_file.geneid.gff3" ) {
    &filter_gff2_geneid(
                         "$gff_file.geneid.golden.gff3",
                         $training_genes_hashref,
                         "$gff_file.geneid.golden.train.gff3",
                         1
    );
    rename( "$gff_file.geneid.golden.train.gff3.rest",
            "$gff_file.geneid.golden.test.gff3" );
    &order_fasta($genome_sequence_file,
                  "$gff_file.geneid.golden.gff3" );
    &order_fasta($genome_sequence_file,
                  "$gff_file.geneid.golden.train.gff3" );
    &order_fasta($genome_sequence_file,
                  "$gff_file.geneid.golden.test.gff3" );

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
    &gb2geneid( $train_ref,
                "$gff_file.golden.train.good.gb.geneid" )
      if $train_ref;
    &gb2geneid( $test_ref,
                "$gff_file.golden.test.good.gb.geneid" )
      if $test_ref;
    &gb2geneid( $opt_ref,
                "$gff_file.golden.optimization.good.gb.geneid" )
      if $opt_ref;
    &order_fasta("$gff_file.golden.optimization.good.gb.fasta","$gff_file.golden.optimization.good.gb.geneid" );
    &order_fasta( "$gff_file.golden.train.good.gb.fasta", "$gff_file.golden.train.good.gb.geneid" );
    &order_fasta( "$gff_file.golden.test.good.gb.fasta", "$gff_file.golden.test.good.gb.geneid" );
  }

  #glimmer
  print "\tglimmer\n";
  &gb2glimmer( $train_ref,
               "$gff_file.golden.train.good.gb.glimmer" )
    if $train_ref;
  &gb2glimmer( $test_ref,
               "$gff_file.golden.test.good.gb.glimmer" )
    if $test_ref;
  &gb2glimmer( $opt_ref,
               "$gff_file.golden.optimization.good.gb.glimmer" )
    if $opt_ref;

  #Different thresholds that can be chosen for the splice sites
  #can be consulted in: - false.acc  false.don  false.atg


}

sub revcomp {
  my $dna     = shift;
  my $revcomp = reverse($dna);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}

sub filter_gff() {
  my $gff           = shift;
  my $hash_ref      = shift;
  my $out           = shift;
  my $overlap_check = shift;
  my $out2          = $out . ".rest";
  my $filter_out    = "$gff.filtered";
  my $filter_out2   = "$gff.filtered.rest";
  open( GFF, $gff ) || die;
  open( OUT, ">$filter_out" );
  open( OUT2, ">$filter_out2" );
  my $rejects;
  while ( my $ln = <GFF> ) {
    my $id;
    next unless $ln =~ /\sgene\s/;
    $ln =~ /ID=([^;]+);?/;
    $id = $1 || die $ln;
    chomp($id);

    #pasa
    my $mRNA_id = $id;
    if ( $id =~ /^g\.(\d+)/ ) {
      $mRNA_id = 'm.' . $1;
    }
    if ( $hash_ref->{$id} || $hash_ref->{$mRNA_id} ) {
      print OUT $ln;
      while ( my $ln2 = <GFF> ) {
        if ( $ln2 =~ /^\s*$/ || $ln2 =~ /^##\n/ ) {
          print OUT "\n";
          last;
        }
        print OUT $ln2;
      }
    } else {
      print OUT2 $ln;
      while ( my $ln2 = <GFF> ) {
        if ( $ln2 =~ /^\s*$/  || $ln2 =~ /^##\n/ ) {
          print OUT2 "\n";
          last;
        }
        print OUT2 $ln2;
      }
    }
  }
  close GFF;
  close OUT;
  close OUT2;
  if ($overlap_check) {
    $rejects = &remove_overlapping_sorted_gff( $filter_out,  $out );
    &remove_overlapping_sorted_gff( $filter_out2, $out2 );
    unlink($filter_out);
    unlink($filter_out2);
  } else {
    rename( $filter_out,  $out );
    rename( $filter_out2, $out2 );
  }
  return $rejects;
}

sub filter_gff2_geneid() {
  my $gff           = shift;
  my $hash_ref      = shift;
  my $out           = shift;
  my $overlap_check = shift;
  $out = "$gff.filtered" if !$out;
  my $out2        = $out . ".rest";
  my $filter_out  = "$gff.filtered";
  my $filter_out2 = "$gff.filtered.rest";

  my $previous_ref='';
  open( GFF, $gff ) || die;
  open( OUT, ">$filter_out" );
  open( OUT2, ">$filter_out2" );
  while ( my $ln = <GFF> ) {
    next if $ln =~ /^\s*$/;
    $ln =~ /(\S+)$/;
    my $id = $1 || die;
    $id=~s/\.mRNA$//;
    $ln =~ /^(\S+)/;
    my $ref = $1;
    $previous_ref = $ref if !$previous_ref;
    next if $ref =~/#$/ && $previous_ref =~/#$/;
    if ( $hash_ref->{$id} ) {
      if ( $ref ne $previous_ref ) {
        $previous_ref = $ref;
        print OUT '#$' . "\n";
      }
      print OUT $ln;
    } else {
      if ( $ref ne $previous_ref ) {
        $previous_ref = $ref;
        print OUT2 '#$' . "\n";
      }
      print OUT2 $ln;
    }
  }
  close GFF;
  close OUT;
  close OUT2;
  if ($overlap_check) {
    &remove_overlapping_gff_geneid( $filter_out,  $out );
    &remove_overlapping_gff_geneid( $filter_out2, $out2 );
    unlink($filter_out);
    unlink($filter_out2);
  } else {
    rename( $filter_out,  $out );
    rename( $filter_out2, $out2 );
  }
}

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

sub gff2zff() {
  ## NB snap uses the word exon to mean CDS... argh
  # don't use hardmasked.
  my $gff = shift;
  my $out = shift;
  my $out_hints = $out.".xdef";
  open( GFF, $gff ) || die;
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
      $zff_line .= ">$ref\n";
      $xdef_line.= ">$ref\n";
      $previous_ref = $ref;
    }
    $mRNA_data[8] =~ /ID=([^;]+)/;
    my $mRNA_id = $1 || die Dumper $record;
    my $cds_counter = 0;
    my @cds_data;
    foreach my $line (@lines) {
      next if $line =~ /^\s*$/ || $line=~/^##/;
      my @data = split( "\t", $line );
      die "GFF $gff had a line that didn't begin with the refence $ref".$line unless $data[0] eq $ref;
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

      } elsif ( $i == 0 && $data[6] eq '+' ) {
        $type = 'Einit';
      } elsif ( $i == ( $cds_counter - 1 ) && $data[6] eq '+' ) {
        $type = 'Eterm';
      } else {
        $type = 'Exon';
      }

      $zff_line .= join( "\t", ( $type, $start, $end, $mRNA_id ) ) . "\n";#unless $type eq 'Esngl'; 
      $xdef_line .=join( "\t", ( "Coding", $data[3], $data[4], $data[6], '+100','.','.','.','ADJ' ) ) . "\n";#unless $type eq 'Esngl';
    }
    print OUT $zff_line if $zff_line;
    print XDEF $xdef_line if $xdef_line;
  }
  $/ = $orig_delim;
  close OUT;
  close XDEF;
  close GFF;
}

sub sort_gff3() {
  my $gff       = shift;
  my $delimiter = shift;
  my $orig_sep  = $/;
  open( GFF, $gff ) ||die;
  open( OUT, ">$gff.sorted" );
  $/ = $delimiter;
  my @records = <GFF>;
  chomp(@records);

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
  close GFF;
  close OUT;
  $/ = $orig_sep;
  rename( "$gff.sorted", $gff );
}

sub order_fasta() {
  my $fasta = shift || die;
  my $gff   = shift || die;
  &process_cmd("$cdbfasta_exec $fasta 2>/dev/null") unless -s $fasta.'.cidx';
  my @order = `cut -f 1 $gff|uniq `;
  open( OUT, ">$gff.fasta" );
  my %done;

  foreach my $ln (@order) {
    next if $ln =~ /^#/ || $ln =~ /^\s*$/;
    $ln =~ /^(\S+)/;
    my $id = $1;
    next unless $id;
    if ( !$done{$id} ) {
      my $seq = &get_id_seq_from_fasta($id,$fasta);
      die "Cannot get $id from $fasta for $gff\n" unless $seq;
      print OUT ">$id\n$seq\n";
      $done{$id}++;
    }
  }
  close OUT;
}

sub remove_overlapping_gff_geneid() {
  my $file = shift;
  my $out  = shift;
  open( IN,   $file );
  open( OUT,  ">$out" );
  open( LOG2, ">$out.log" );
  my $previous_ref;
  my $previous_gene_id;
  my $last_position = 1;
  my $skipped       = int(0);
  my %hash;
  my %ignore;
  my @order;

  while ( my $ln = <IN> ) {
    if ( $ln =~ /^#/ ) {
      if (@order) {
        my $last_gene = $order[-1] || next;
        next if !$order[-2];
        if (    !$hash{$last_gene}{'checkF'}
             && !$hash{$last_gene}{'checkT'}
             && $hash{$last_gene}{'ref'} ne $hash{ $order[-2] }{'ref'} )
        {
          print LOG2 "Excluding singleton $last_gene\n";
          pop(@order);
          delete( $hash{$last_gene} );
        }
      }
      next;
    }
    my @data = split( "\t", $ln );
    next if !$data[8];
    chomp( $data[8] );
    next if $ignore{ $data[8] };
    push( @order, $data[8] ) if !$hash{ $data[8] };
    $hash{ $data[8] }{'checkF'}++ if $data[2] eq 'First';
    $hash{ $data[8] }{'checkT'}++ if $data[2] eq 'Terminal';

    if ( $hash{ $data[8] }{'checkF'} && $hash{ $data[8] }{'checkF'} > 1 ) {
      print LOG2 "Gene $data[8] has more than 2 first exons! Excluding\n";
      $ignore{ $data[8] } = 1;
      next;
    } elsif ( $hash{ $data[8] }{'checkT'} && $hash{ $data[8] }{'checkT'} > 1 ) {
      print LOG2 "Gene $data[8] has more than 2 terminal exons! Excluding\n";
      $ignore{ $data[8] } = 1;
      next;
    }
    $hash{ $data[8] }{'ref'} = $data[0];
    push( @{ $hash{ $data[8] }{'data'} }, $ln );
  }
  foreach my $gene_id (@order) {
    next if !$gene_id || !$hash{$gene_id};
    my @gene_data_array = @{ $hash{$gene_id}{'data'} };
    my ( $ref_id, $stop, $start );
    foreach my $ln (@gene_data_array) {
      my @data = split( "\t", $ln );
      $stop  = $data[3] if ( !$stop  || $stop < $data[3] );
      $stop  = $data[4] if ( !$stop  || $stop < $data[4] );
      $start = $data[3] if ( !$start || $start > $data[3] );
      $start = $data[4] if ( !$start || $start > $data[4] );
      $ref_id = $data[0];
    }
    $last_position    = $stop    if !$last_position;
    $previous_gene_id = $gene_id if !$previous_gene_id;
    $previous_ref     = $ref_id  if !$previous_ref;
    if ( $previous_ref eq $ref_id && $start <= $last_position ) {
      print LOG2 "Overlapping gene found $gene_id. Skipping it\n";
      $skipped++;
      next;
    }
    if ( $previous_ref ne $ref_id ) {
      print OUT '#$' . "\n";
    }
    print OUT join( '', @gene_data_array );
    $previous_gene_id = $gene_id;
    $previous_ref     = $ref_id;
    $last_position    = $stop;
  }
  close IN;
  close OUT;
  close LOG2;
  print "\tOverlaps checked. Skipped $skipped genes. See $out.log for details\n";
  return $skipped;
}

sub remove_overlapping_sorted_gff() {
 # ideally we could just harvest the gene and check for start/stop but in practice it didn't work 
# as i had to loop through every gene (nested as well)
  my $file = shift;
  my $out  = shift;
  my %hash;
  open( IN,   $file );
  open( OUT,  ">$out" );
  open( LOG2, ">$out.log" );
  my $previous_ref;
  my $previous_gene_id;
  my $last_position = 1;
  my $toprint       = '';
  my $skipped       = int(0);

  while ( my $ln = <IN> ) {
    if ( $ln =~ /^#/ || $ln =~ /^\s*$/ ) { $toprint .= $ln; next; }
    my @data = split( "\t", $ln );
    if ( $data[2] eq 'gene' ) {
      my $gene_data = '';
      while ( my $ln2 = <IN> ) {
        $gene_data .= $ln2;
        last if ( $ln2 =~ /^\s*$/ );
      }
      my $gene_id = $data[8] || die;
      my $ref_id  = $data[0] || die;
      chomp($gene_id);
      $previous_gene_id = $gene_id if !$previous_gene_id;
      $previous_ref     = $ref_id  if !$previous_ref;
      if ( $previous_ref eq $ref_id
           && ( $data[3] <= $last_position || $data[4] <= $last_position ) )
      {
        print LOG2 "Overlapping gene found $gene_id. Skipping it and keeping $previous_gene_id\n";
        $skipped++;
        $toprint = '';
        next;
      }
      $toprint .= $ln . $gene_data;
      print OUT $toprint;
      $toprint          = '';
      $previous_gene_id = $gene_id;
      $previous_ref     = $ref_id;
      $last_position    = $data[4] > $data[3] ? $data[4] : $data[3];
    }
  }
  print OUT $toprint;
  close IN;
  close OUT;
  close LOG2;
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
          push(
                @{ $hash{$id}{'exons'}{$exon_count} },
                ( $loc->start, $loc->end )
          );
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
    my @sorted_exons = sort { $a <=> $b } keys %{ $hash_ref->{$id}{'exons'}};
    if ( $hash_ref->{$id}{'strand'} == 1 ) {
      my $strand = '+';
      my ($gstart,$e) = @{ $hash_ref->{$id}{'exons'}->{$sorted_exons[0]} };
      my ($s,$gend) = @{ $hash_ref->{$id}{'exons'}->{$sorted_exons[-1]} };
      print  GFF $id
          . "\tGB\texon\t"
          . $gstart . "\t"
          . $gend
          . "\t.\t$strand\t.\t$id.1" . "\n";
      
      foreach my $exon_counter (@sorted_exons)      {
        my ( $start, $stop ) = @{ $hash_ref->{$id}{'exons'}->{$exon_counter} };
        print GFF $id
          . "\tGB\texon\t"
          . $start . "\t"
          . $stop
          . "\t.\t$strand\t.\t$id.1" . "\n";
      }
    } elsif ( $hash_ref->{$id}{'strand'} == -1 ) {
      my $strand = '-';
      my ($gstart,$e) = @{ $hash_ref->{$id}{'exons'}->{$sorted_exons[-1]} };
      my ($s,$gend) = @{ $hash_ref->{$id}{'exons'}->{$sorted_exons[0]} };
      print  GFF $id
          . "\tGB\texon\t"
          . $gstart . "\t"
          . $gend
          . "\t.\t$strand\t.\t$id.1" . "\n";
      foreach my $exon_counter (@sorted_exons){
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
  &sort_gff3( $out,"\n\n" );
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
        } elsif ( $exon_counter == 1 ) {
          $exon_type = "First";
        } elsif (
            $exon_counter == ( scalar( keys %{ $hash_ref->{$id}{'exons'} } ) ) )
        {
          $exon_type = 'Terminal';
        } else {
          $exon_type = 'Internal';
        }
        my ( $start, $stop ) = @{ $hash_ref->{$id}{'exons'}->{$exon_counter} };
        print GENEID $id
          . "\tFORGENEID\t$exon_type\t"
          . $start . "\t"
          . $stop
          . "\t.\t$strand\t.\t$id.1" . "\n";
      }
    } elsif ( $hash_ref->{$id}{'strand'} == -1 ) {
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
        } elsif ( $exon_counter == 1 ) {
          $exon_type = "Terminal";
        } elsif (
            $exon_counter == ( scalar( keys %{ $hash_ref->{$id}{'exons'} } ) ) )
        {
          $exon_type = 'First';
        } else {
          $exon_type = 'Internal';
        }
        my ( $start, $stop ) = @{ $hash_ref->{$id}{'exons'}->{$exon_counter} };
        print GENEID $id
          . "\tFORGENEID\t$exon_type\t"
          . $start . "\t"
          . $stop
          . "\t.\t$strand\t.\t$id.1" . "\n";
      }
    }
    print GENEID '#$' . "\n";
  }
  close GENEID;
  &sort_gff3( $out, '#$' . "\n" );
}

sub gb2glimmer() {
  my $hash_ref = shift;
  my $out      = shift;
  open( GLIMMER, ">$out" );
  foreach my $id ( sort keys %{$hash_ref} ) {
    next unless $hash_ref->{$id}{'strand'};
    if ( $hash_ref->{$id}{'strand'} == 1 ) {
      foreach
        my $exon ( sort { $a <=> $b } keys %{ $hash_ref->{$id}{'exons'} } )
      {
        my ( $start, $stop ) = @{ $hash_ref->{$id}{'exons'}->{$exon} };
        print GLIMMER $id . " " . $start . " " . $stop . "\n";
      }
    } elsif ( $hash_ref->{$id}{'strand'} == -1 ) {
      foreach
        my $exon ( sort { $b <=> $a } keys %{ $hash_ref->{$id}{'exons'} } )
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
  pod2usage "Can't find the Augustus directory, easiest way is to create a symlink of the augustus executable somewhere in your PATH (e.g. \$HOME/bin) or add the Augustus bin directory in your PATH\n" unless $augustus_dir && -d $augustus_dir;
  $gff2gb_exec = $augustus_dir . '/scripts/gff2gbSmallDNA.pl' if !$gff2gb_exec;
  $augustus_train_exec = $augustus_dir . '/bin/etraining' if !$augustus_train_exec;
  $augustus_filterGenes_exec = $augustus_dir.'/scripts/filterGenes.pl' if !$augustus_filterGenes_exec;
  pod2usage "Can't find etraining from Augustus\n" if !-s $augustus_train_exec;
  pod2usage "Can't find gff2gbSmallDNA.pl from Augustus\n" if !-s $gff2gb_exec;
  pod2usage "Can't find filterGenes.pl from Augustus\n" if !-s $augustus_filterGenes_exec;
}

sub parse_transdecoder_gene_gff3() {

  # add utr co-ordinates to GFF
  my $transdecoder_gff = shift;
  my $fasta            = shift;
  my $hash_ref         = shift;
  my $orig_sep = $/;
  my %seq_hash;
  $/ = ">";
 open (IN,$fasta);
 <IN>;
 while (my $record=<IN>){
	my @lines = split("\n",$record);
	pop(@lines) if $lines[-1] eq '>';
        my $id = shift (@lines);
	$seq_hash{$id} = join('',@lines);
 }
 close IN;
  my $transdecoder_contig = $transdecoder_gff . '.contigs';
  open( ANNOT, ">$transdecoder_contig.annotations" );
  my %contigs_used;
  $/ = "\n\n";
  open( TRANSDECODER, $transdecoder_gff ) || die;
  open(TROUT,"> $transdecoder_contig");
  while ( my $record = <TRANSDECODER> ) {
    chomp($record);
    my @lines     = split( "\n", $record );
    my $mRNA_line = $lines[1];
    my @mRNA_data = split( "\t", $mRNA_line );
    $mRNA_data[8] =~ /ID=([^;]+)/;
    my $mRNA_id = $1 || next;
    next unless !$only_complete || $mRNA_data[8] =~ /type%3Acomplete/;
    next if ( $hash_ref && !$hash_ref->{$mRNA_id} );
    my $ref = $mRNA_data[0];
    $contigs_used{$ref}++;
    print TROUT ">$mRNA_id\n".$seq_hash{$ref}."\n";
    for ( my $i = 2 ; $i < scalar(@lines) ; $i++ ) {
      my @data = split( "\t", $lines[$i] );
      print ANNOT join(
                        ' ',
                        (
                          $mRNA_id, $data[6],
                          $data[3], ( abs( $data[4] - $data[3] ) + 1 )
                        )
        )
        . "\n"
        if $data[2] eq 'CDS';
    }
  }
  close TRANSDECODER;
  close ANNOT;
  close TROUT;
  $/ = $orig_sep;
  print "Prepared $transdecoder_contig with "
    . scalar( keys %contigs_used )
    . " sequences\n";
  return $transdecoder_contig;
}

sub splitfasta() {
  my @files;
  my $file2split         = shift;
  my $outdir             = shift;
  my $how_many_in_a_file = shift;
  return unless $file2split && -s $file2split && $outdir && $how_many_in_a_file;
  mkdir($outdir) unless -d $outdir;
  my ($flag);
  my $filecount = 0;
  my $seqcount  = 0;
  open( FILE, $file2split );

  while ( my $line = <FILE> ) {
    if ( $line =~ /^\s*$/ ) { next; }    #empty line
    elsif ( $line =~ /^>(\S+)/ ) {
      my $id = $1;
      $seqcount++;
      if ( !$flag ) {
        $filecount++;
        my $outfile = $file2split;
        $outfile .= "_" . $filecount;
        $outfile = $id if $how_many_in_a_file == 1;
        open( OUT, ">$outdir/$outfile" ) || die("Cannot open $outdir/$outfile");
        push( @files, "$outdir/$outfile" );
        $flag = 1;
      } elsif ( $seqcount >= $how_many_in_a_file ) {
        $seqcount = 0;
        $filecount++;
        my $outfile = $file2split;
        $outfile .= "_" . $filecount;
        $outfile = $id if $how_many_in_a_file == 1;
        close(OUT);
        open( OUT, ">$outdir/$outfile" );
        push( @files, "$outdir/$outfile" );
      }
      print OUT ">$id\n";
    } else {
      print OUT $line;
    }
  }
  close(FILE);
  close(OUT);
  return \@files;
}

sub gff2hints() {
  my $gff       = shift;
  my $golden    = shift;
  my $src       = $golden ? 'GLD' : 'XNT';
  my $type_suffix = $golden ? '' : 'part';
  my $priority = $golden ? 7 : 5;
  my $delimiter = "\n\n";
  my $orig_sep  = $/;
  open( GFF, $gff )          || die;
  open( OUT, ">$gff.hints" ) || die;
  $/ = $delimiter;

  while ( my $record = <GFF> ) {
    my @lines = split( "\n", $record );

    #my $gene_line = $lines[0];
    my $mRNA_line = $lines[1];
    $mRNA_line =~ /ID=([^;]+)/;
    my $mRNA_id = $1;
    $mRNA_id=~s/[^\w\.]+/_/g;
    # CDS intron
    for ( my $i = 2 ; $i < scalar(@lines) ; $i++ ) {
      my @data = split( "\t", $lines[$i] );
      next unless $data[8];
      if ( $data[2] eq 'mRNA' ) {
        $data[2] = 'genicpart';
        print OUT join( "\t", @data[ 0 .. 7 ] )
          . "\tsrc=$src;pri=$priority;grp=$mRNA_id\n";
      } elsif (    $data[2] eq 'exon'
                || $data[2] eq 'CDS'
                || $data[2] eq 'intron' )
      {
        $data[2].=$type_suffix;
        print OUT join( "\t", @data[ 0 .. 7 ] )
          . "\tsrc=$src;pri=$priority;grp=$mRNA_id\n";
      } elsif ( $data[2] eq 'splice_junction' ) {
        $data[2] = $data[8] =~ /splice3/ ? 'ass' : 'dss';
        print OUT join( "\t", @data[ 0 .. 7 ] )
          . "\tsrc=$src;pri=$priority;grp=$mRNA_id\n";
      } elsif ( $data[2] eq 'UTR' ) {
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
  my ( $cmd) = @_;
  my $ret = system($cmd);
  if ( $ret && $ret != 256 ) {
    die "Error, cmd died with ret $ret\n";
  }
  return $ret;
}


sub gff3_fix_phase_GTF(){
    my $gff3_file = shift;
    open (IN,$gff3_file);
    open(OUT,">$gff3_file.cds");
    while (my $ln=<IN>){
	print OUT $ln unless $ln=~/\tintron\t|\tsplice_junction\t/;
    }
    close IN;
    close OUT;
    $gff3_file .= '.cds';
    my $index_file = "$gff3_file.inx";
    my $gene_obj_indexer = new Gene_obj_indexer( { "create" => $index_file } );
    my $asmbl_id_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer);
    open (OUT,">$gff3_file.t");
    open (GTF,">$gff3_file.gtf");
    foreach my $asmbl_id (keys %$asmbl_id_to_gene_list_href) {
       my $genome_seq = &get_id_seq_from_fasta($asmbl_id, $genome_sequence_file);
       my @gene_ids = @{$asmbl_id_to_gene_list_href->{$asmbl_id}};
       foreach my $gene_id (@gene_ids) {
           my $gene_obj_ref = $gene_obj_indexer->get_gene($gene_id);
           foreach my $isoform ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
               $isoform->delete_isoforms();
               my $isoform_id = $isoform->{Model_feat_name};
               eval {
                  $isoform->set_CDS_phases (\$genome_seq);
               };
               print OUT $isoform->to_GFF3_format() . "\n";
                my $gtf_text = "";
                eval {
                    $gtf_text = $isoform->to_GTF_format(\$genome_seq);
                };
                if ($@) {
                    # do it in pseudogene mode
                    $isoform->{is_pseudogene} = 1;
                    $gtf_text = $isoform->to_GTF_format(\$genome_seq);
                }
                print GTF "$gtf_text\n";
           }
       }
   }
   close OUT;
   close GTF;
   unlink("$gff3_file");
   rename("$gff3_file.t","$gff3_file");
   unlink $index_file;
}


sub convert_to_gtf(){

    my $gff3_file = shift;
    open (OUT,">$gff3_file.gtf");
    my $inx_file = "$gff3_file.inx";
    my $gene_obj_indexer = new Gene_obj_indexer( { "create" => $inx_file } );
    
    my $asmbl_id_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer);

    foreach my $asmbl_id (sort keys %$asmbl_id_to_gene_list_href) {
        
        ## get the genome sequence
        my $genome_seq = &get_id_seq_from_fasta($asmbl_id, $genome_sequence_file);
        my @gene_ids = @{$asmbl_id_to_gene_list_href->{$asmbl_id}};
        foreach my $gene_id (@gene_ids) {
            ## note models of isoforms are bundled into the same gene object.
            my $gene_obj_ref = $gene_obj_indexer->get_gene($gene_id);
            foreach my $gene_obj ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
                $gene_obj->delete_isoforms(); # unbundle the model object!
                my $gtf_text = "";
                eval {
                    $gtf_text = $gene_obj->to_GTF_format(\$genome_seq);
                };
                if ($@) {
                    # do it in pseudogene mode
                    $gene_obj->{is_pseudogene} = 1;
                    $gtf_text = $gene_obj->to_GTF_format(\$genome_seq);
                }
                print OUT "$gtf_text\n";
            }
        }
    }
    close OUT;
    unlink $inx_file;
}


sub check_program(){
	my @paths;
	foreach my $prog (@_){
	        my $path = `which $prog`;
        	die "Error, path to required $prog cannot be found" unless $path =~ /^\//;
		chomp($path);
		$path = readlink($path) if -l $path;
		push(@paths,$path);
	}
	return @paths;
}

sub check_program_optional(){
        my @paths;
        foreach my $prog (@_){
                my $path = `which $prog`;
                warn "Warning: path to optional $prog cannot be found" unless $path =~ /^\//;
                chomp($path);
		$path = readlink($path) if -l $path;
                push(@paths,$path);
        }
        return @paths;
}


sub get_id_seq_from_fasta {
    my ($acc, $fasta_db) = @_;
    my $seq;
    if (!$get_id_seq_from_fasta_hash{$fasta_db}{$acc}){
       $seq = `$cdbyank_exec $fasta_db.cidx -a '$acc'`;
       chomp($seq);
       $get_id_seq_from_fasta_hash{$fasta_db}{$acc}=$seq;
    }else{
       $seq = $get_id_seq_from_fasta_hash{$fasta_db}{$acc};
    }
    die "FATAL: couldn't retrieve seq for $acc from $fasta_db\n" unless ($seq);
    my @x = split (/\n/, $seq);
    my $id = shift @x;
    $seq = join ("", @x);
    $seq = uc ($seq);
    $seq =~ s/\s//g;
    return $seq;
}


sub run_exonerate(){
if ( ( $peptide_file || ( $transdecoder_gff && $transdecoder_assembly_file )  )
     && $genome_file )
{
    my ($aatpackage,$parafly_exec) = &check_program('AAT.pl','ParaFly');
    $aatpackage = dirname($aatpackage) if $aatpackage;
  my $genome_dir = basename($genome_file) . '_dir';
  my $files_ref;
  if ( -d $genome_dir ) {
    my @genome_files = glob( $genome_dir . "/*" );
    $files_ref = \@genome_files;
  } else {
    print "Splitting genome into $genome_dir...\n";
    $files_ref = &splitfasta( $genome_file, $genome_dir, 1 );
  }
  $exonerate_file =
      $transdecoder_gff
    ? $transdecoder_gff . '.exonerate.results'
    : $peptide_file . '.exonerate.results';
  my @commands;
  my $aat_command_file = "./".basename($genome_dir).".commands";
  # check if previous exonerate run exists and if so run it?
    unless ($no_rerun_exonerate){
        if (-s 'run_exonerate_commands.cmd' &&(!-s 'run_exonerate_commands.cmd.completed' || (-s 'run_exonerate_commands.cmd' != -s 'run_exonerate_commands.cmd.completed'))){
            print "Processing with AAT\n";
  	    &process_cmd("$parafly_exec -shuffle -CPU $threads -c run_exonerate_commands.cmd -failed_cmds run_exonerate_commands.cmd.failed -v 2>&1 |grep -v 'successfully completed'");
	}
  }
  if ($peptide_file) {
    unless ($no_rerun_exonerate){

    print "Preparing/running AAT...\n";
    my $matrix_file = "$aatpackage/../matrices/BS";
    die "Cannot find AAT's matrices/BS as $matrix_file\n" unless -s $matrix_file;
    foreach my $genome_file (@$files_ref) {
      next if $genome_file =~ /\.aat\./ || -d $genome_file;
      push( @commands,"$aatpackage/dps $genome_file $peptide_file $aatpackage/../matrices/BS  -f 400 -i 30 -a $intron_size > $genome_file.aat.d ;"
          . "$aatpackage/ext $genome_file.aat.d -f 400 > $genome_file.aat.ext ;"
          . "$aatpackage/extCollapse.pl $genome_file.aat.ext > $genome_file.aat.extCol ;"
          . "$aatpackage/filter $genome_file.aat.extCol -c 1 > $genome_file.aat.filter ; rm -f $genome_file.aat.d $genome_file.aat.ext $genome_file.aat.extCol \n"
      ) unless -s "$genome_file.aat.filter";
    }
    open( CMD, ">$aat_command_file" );
    print CMD shuffle(@commands);
    close CMD;
    &process_cmd("$parafly_exec -shuffle -CPU $threads -c $aat_command_file -v -failed_cmds $aat_command_file.failed"
    );
    system("find $genome_dir -empty -delete");
    print "Running exonerate...\n";
    unlink($exonerate_file);
    my $exonerate_options = " -minorf $minorf -protein -in $peptide_file -separate -aat $genome_dir/*filter -threads $threads -intron_max $intron_size  $same_species ";
    $exonerate_options .=" -softmask -ref $softmasked_genome " if ($softmasked_genome);
    $exonerate_options .=" -ref $genome_file " if (!$softmasked_genome);
    $exonerate_options .=" -norefine " if $norefine;
    &process_cmd('run_exonerate.pl'.$exonerate_options);
    die "Exonerate run failed for some reason....\n"
      if ( !-d $peptide_file . "_queries" );
    system(   "cat "
            . $peptide_file
            . "_queries/*exonerate_results >> $exonerate_file" );
    }
  } else {
    my ($aatpackage,$parafly_exec) = &check_program('AAT.pl','ParaFly');
    $aatpackage = dirname($aatpackage) if $aatpackage;
    print "Processing transdecoder output via AAT and exonerate to produce output $exonerate_file\n";
    print "Finding contigs from transdecoder...\n";
    my $fasta_contigs = "$transdecoder_gff.contigs";
    if ($transdecoder_peptides) {
      $is_cdna = 1;
      my %hash;
      open( IN, $transdecoder_peptides );
      while ( my $ln = <IN> ) {
        next unless $ln =~ /^>(\S+)/;
        $hash{$1} = 1;
      }
      close IN;
      $fasta_contigs =
        &parse_transdecoder_gene_gff3( $transdecoder_gff,
                                       $transdecoder_assembly_file, \%hash )
        unless -s $fasta_contigs;
    } else {
      $fasta_contigs =
        &parse_transdecoder_gene_gff3( $transdecoder_gff,
                                       $transdecoder_assembly_file )
        unless -s $fasta_contigs;
    }
    die "Could not create $fasta_contigs\n"
      unless -s $fasta_contigs && -s $fasta_contigs . '.annotations';

    unless ($no_rerun_exonerate){

    print "Preparing/running AAT...\n";
    foreach my $genome_file (@$files_ref) {
      next if $genome_file =~ /\.aat\./ || -d $genome_file;
      push(
        @commands,"$aatpackage/dds $genome_file $fasta_contigs -o 75 -p 75 -c 30000 -f 200 -i 30 -a $intron_size > $genome_file.aat.d ;"

#   "$aatpackage/dps $genome_file $transdecoder_peptides $aatpackage/../matrices/BS  -f 400 -i 30 -a $intron_size > $genome_file.aat.d ;"
          . "$aatpackage/ext $genome_file.aat.d -f 400 > $genome_file.aat.ext ;"
          . "$aatpackage/extCollapse.pl $genome_file.aat.ext > $genome_file.aat.extCol ;"
          . "$aatpackage/filter $genome_file.aat.extCol -c 1 > $genome_file.aat.filter ; rm -f $genome_file.aat.d $genome_file.aat.ext $genome_file.aat.extCol \n"
      ) unless -s "$genome_file.aat.filter";
    }
      open( CMD, ">$aat_command_file" );
      print CMD shuffle(@commands);
      close CMD;
      &process_cmd("$parafly_exec -shuffle -CPU $threads -c $aat_command_file -v -failed_cmds $aat_command_file.failed"
      );
    system("find $genome_dir -empty -delete");
      print "Running exonerate...\n";
      my $exonerate_options = " -minorf $minorf -annotation $fasta_contigs.annotations -in $fasta_contigs -separate -aat $genome_dir/*filter -threads $threads -intron_max $intron_size $same_species ";
      $exonerate_options .=" -softmask -ref $softmasked_genome " if ($softmasked_genome);
      $exonerate_options .=" -ref $genome_file " if (!$softmasked_genome);
      $exonerate_options .=" -norefine " if $norefine;
      &process_cmd('run_exonerate.pl'.$exonerate_options);;
      die "Exonerate run failed for some reason....\n"
        if ( !-d $fasta_contigs . "_queries" );
      unlink ($exonerate_file);
      system(   "cat "
            . $fasta_contigs
            . "_queries/*exonerate_results >> $exonerate_file" );
    }

  print "Completed processing exonerate file as $exonerate_file\n";
  }
}
}
