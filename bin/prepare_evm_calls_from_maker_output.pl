#!/usr/bin/env perl

=pod

=head1 USAGE

Use EVM to call new Maker models

Mandatory:

  -maker_output_dir   :s Directory that hosts maker output
  -genome             :s FASTA file with genome sequences

Optional:

  -weights            :s A weight file for EVM. Defaults created
  -verbose               Be slighly verbose
  
 =head1 AUTHORS

 Brian Haas, Alexie Papanicolaou
  
=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use POSIX qw(ceil);
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/:$RealBin:$RealBin/../3rd_party/evidencemodeler/";
use Gene_obj;
use Gene_validator;

my ( $maker_dir, $genome_sequence, $evm_weights_file, $debug );
my %allowed_types =
  ( 'augustus_masked' => 1, 'genemark' => 1, 'snap_masked' => 1 );
my $outdir = 'maker_output.' . $$;
my ($evm_exec) = &check_program("evidence_modeler.pl");

pod2usage $! unless GetOptions(
            'maker_output_dir:s' => \$maker_dir,
            'genome:s'           => \$genome_sequence,
            'weights'            => \$evm_weights_file,
            'debug|verbose'      => \$debug
);

pod2usage "Cannot find $maker_dir\n" unless $maker_dir && -d $maker_dir;
pod2usage "No genome provided...\n"
  unless $genome_sequence && -s $genome_sequence;
$evm_weights_file = &create_evm_weights()
  unless -s $evm_weights_file && $evm_weights_file;
mkdir $outdir unless -d $outdir;

my $maker_gff = $maker_dir . '.gff';

&process_cmd("find . -name '*.gff' -exec cat '{}' \+ | cut -s -f 1-9 | grep -vP '\tcontig\t1'  > $maker_gff");

open( IN, $maker_gff ) || die;
my %hash;
while ( my $ln = <IN> ) {
 next if $ln =~ /^#/;
 my @data = split( "\t", $ln );
 next unless $data[8];
 $ln = "\n" . $ln if $data[2] eq 'gene';
 push( @{ $hash{ $data[1] } }, $ln );
}
close IN;

my (@abinitio_files);
foreach my $type ( keys %hash ) {
 open( OUT, ">$outdir/$type" );
 foreach my $ln ( @{ $hash{$type} } ) {
  print OUT $ln;
 }
 close OUT;
 next unless $allowed_types{$type};
 my $gene_gff = &maker_match_gff_to_gene_gff3("$outdir/$type");
 push( @abinitio_files, $gene_gff );
}

%hash = ();
open( OUTGP, ">$outdir/gene_predictions.gff3" );
foreach my $file (@abinitio_files) {
 open( IN, $file );
 while ( my $ln = <IN> ) {
  print OUT $ln;
 }
 print "\n";
 close IN;
}
close OUTGP;

my $evm_output = $outdir . '/evm.out';
&process_cmd(   $evm_exec
              . " --gene_predictions $outdir/gene_predictions.gff3 "
              . " --transcript_alignments $outdir/est2genome.gff "
              . " --protein_alignments $outdir/protein2genome.gff "
              . " --weights $evm_weights_file "
              . " --genome $genome_sequence "
              . " --trellis_search_limit 200 "
              . " > $evm_output " );

print "Done.\n See $evm_output\n\n";

##########################################################
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

sub evm_to_gff_gtf() {

 my $evm_file                 = shift;
 my $genome_sequences_hashref = &read_fasta($genome_sequence);
 open( GFFOUT, ">$evm_file.gff" );
 open( GTF,    ">$evm_file.gtf" );

 foreach my $contig_id ( keys @{$genome_sequences_hashref} ) {
  my ( %model_num_to_coords, %model_id_to_ev_type, %end5_to_phase );
  my $model_id = 1;

  open( my $fh, $evm_file ) or die "Error, cannot open $evm_file\n";
  while (<$fh>) {
   if (/^\!/) { next; }    ## comment line

   if (/^\#/) {
    my $ev_type;
    if (/ELIMINATED/) {
     $ev_type = "EVM_elm";
    }
    else {
     $ev_type = "EVM";
    }
    $model_id_to_ev_type{$model_id} = $ev_type;
    next;
   }

   chomp;

   unless (/\w/) {
    $model_id++;
    next;
   }

   my @x = split(/\t/);
   if ( scalar(@x) == 6 && $x[0] =~ /^\d+$/ && $x[1] =~ /^\d+$/ ) {

    if ( $x[2] eq 'INTRON' ) { next; }

    my $coords_ref = $model_num_to_coords{$model_id};
    unless ( ref $coords_ref ) {
     $coords_ref = $model_num_to_coords{$model_id} = {};
    }

    $coords_ref->{ $x[0] } = $x[1];

    $end5_to_phase{ $x[0] } = $x[3];

   }
  }
  close $fh;

  my %phase_conversion = (
                           1 => 0,
                           2 => 1,
                           3 => 2,
                           4 => 0,
                           5 => 1,
                           6 => 2
  );

  my @gene_objs;

  foreach my $model_id ( sort { $a <=> $b } keys %model_num_to_coords ) {

   my $ev_type = $model_id_to_ev_type{$model_id};

   # if ($ev_type ne "EVM") { next; } #ignoring the EVM_elm for now.

   my $coords_ref = $model_num_to_coords{$model_id};

   my $gene_obj = new Gene_obj();
   $gene_obj->populate_gene_obj( $coords_ref, $coords_ref );

   ## set phase:
   foreach my $exon ( $gene_obj->get_exons() ) {
    my $cds_obj = $exon->get_CDS_obj();
    my ( $end5, $end3 ) = $cds_obj->get_coords();
    $cds_obj->set_phase( $phase_conversion{ $end5_to_phase{$end5} } );
   }

   $gene_obj->{Model_feat_name} = "evm.model.$contig_id.$model_id";
   $gene_obj->{TU_feat_name}    = "evm.TU.$contig_id.$model_id";
   $gene_obj->{com_name}        = "EVM prediction $contig_id.$model_id";

   $gene_obj->{asmbl_id} = "$contig_id";

   print GFFOUT $gene_obj->to_GFF3_format( source => $ev_type ) . "\n";

   my $gtf_text = "";
   eval {
    $gtf_text =
      $gene_obj->to_GTF_format( \$genome_sequences_hashref->$contig_id );
   };
   if ($@) {

    # do it in pseudogene mode - if not then UTR is printed as CDS... not good!
    $gene_obj->{is_pseudogene} = 1;
    $gtf_text =
      $gene_obj->to_GTF_format( \$genome_sequences_hashref->$contig_id );
   }
   print GTF "$gtf_text\n";
  }
 }
 close GFFOUT;
 close GTF;
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

sub maker_match_gff_to_gene_gff3 {
 my $match_gff = shift;
 my $out       = $match_gff . '.gene.gff3';

 #from brian;
 my %genes;

 open( my $fh, $match_gff ) or die "Error, cannot open file $match_gff";
 open( OUT, ">$out" );
 while (<$fh>) {
  unless (/\w/)  { next; }
  if     (/^\#/) { next; }
  chomp;

  my @x         = split(/\t/);
  my $scaff     = $x[0];
  my $source    = $x[1];
  my $feat_type = $x[2];
  my $lend      = $x[3];
  my $rend      = $x[4];
  my $orient    = $x[6];
  my $info      = $x[8];

  ( $lend, $rend ) = sort { $a <=> $b } ( $lend, $rend );

  my ( $end5, $end3 ) =
    ( $orient eq '+' ) ? ( $lend, $rend ) : ( $rend, $lend );

  my %tags = &get_tags($info);

  if ( my $parent = $tags{Parent} ) {

   $genes{$parent}->{coords}->{$end5} = $end3;
  }
  elsif ( my $id = $tags{ID} ) {
   my $name = $tags{Name} || "No name";
   $genes{$id}->{name}     = $name;
   $genes{$id}->{scaffold} = $scaff;
   $genes{$id}->{source}   = $source;
  }
 }
 close $fh;

 foreach my $gene_id ( keys %genes ) {
  my $data_href = $genes{$gene_id};

  my $coords_href = $data_href->{coords};

  my $gene_obj = new Gene_obj();
  $gene_obj->populate_gene_object( $coords_href, $coords_href );
  $gene_obj->{asmbl_id}        = $data_href->{scaffold};
  $gene_obj->{com_name}        = $data_href->{name};
  $gene_obj->{TU_feat_name}    = $gene_id;
  $gene_obj->{Model_feat_name} = "m.$gene_id";

  print OUT $gene_obj->to_GFF3_format( source => $data_href->{source} ) . "\n";
 }
 return $out;
}

####
sub get_tags {
 my ($info) = @_;

 my %tags;

 my @fields = split( /;/, $info );
 foreach my $field (@fields) {
  my ( $key, $val ) = split( /=/, $field );

  $tags{$key} = $val;
 }

 return (%tags);
}

sub create_evm_weights() {
 my $out = "$outdir/evm_weights";
 open( OUT, ">$out" );
 print OUT
"ABINITIO_PREDICTION\taugustus_masked\t1\nABINITIO_PREDICTION\tsnap_masked\t1\nABINITIO_PREDICTION\tgenemark\t1\nPROTEIN\tprotein2genome\t1\nTRANSCRIPT\test2genome\t1\n";
 close OUT;
 return $out;
}
