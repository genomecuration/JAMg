#!/usr/bin/env perl

=pod

=head1 NAME

 create_features_from_gff3.pl

=head1 USAGE

Provide a GFF3 file and a genome FASTA file to phase and create sequence feature (CDS,PEP,MRNA,GENE) information for
 
 Mandatory
  -gff|infile   :s GFF file to process
  -genome|fasta :s FASTA file of genome
 
 Optional
  -name            Use Transcript common name as the main ID in the output.
  -lettername      Use -R[two letters] notation for alternative transcript instead of .[digits]  
  -verbose         Print progress and debug info
 
NB: -name means that the common name has no spaces and is unique (will be checked). Useful for WebApollo
 
=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Carp;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use GFF3_utils;
use GTF_utils;
use URI::Escape;

$|=1;
our $SEE;
my $minorf = 3;    #minimum orf size in bp
my ( $gfffile, $genome, $change_name,$lettername,$verbose );
GetOptions(
            'gff|infile:s'     => \$gfffile,
            'genome|fasta:s'   => \$genome,
            'name|change_name' => \$change_name,
            'debug'          => \$SEE,
            'verbose'          => \$verbose,
            'lettername'       => \$lettername
);
$gfffile = shift if !$gfffile;
$genome  = shift if !$genome;
pod2usage unless $gfffile && $genome;

die "GFF file $gfffile does not exist\n"  unless -s $gfffile;
die "Fasta file $genome does not exist\n" unless -s $genome;

my $scaffold_seq_hashref = &read_fasta($genome);

my %unique_names_check;

&gff3_process($gfffile);

#########################################################################
sub gff3_process() {
 my $gff3_file = shift;
 open( IN, $gff3_file ) || confess( "Cannot find $gff3_file " . $! );
 my $index_file = "$gff3_file.inx";
 my $gene_obj_indexer = new Gene_obj_indexer( { "create" => $index_file } );
 my $genome_id_to_gene_list_href =
   &GFF3_utils::index_GFF3_gene_objs( $gff3_file, $gene_obj_indexer );
 open( GFF3, ">$gff3_file.gff3" );
 open( PEP,  ">$gff3_file.pep" );
 open( CDS,  ">$gff3_file.cds" );
 open( GENE, ">$gff3_file.gene" );
 open( MRNA, ">$gff3_file.mRNA" );
 
 print "Indexing complete\n";
 foreach my $genome_id (sort  keys %$genome_id_to_gene_list_href ) {

  my $genome_seq = $scaffold_seq_hashref->{$genome_id};
  if ( !$genome_seq ) {
   warn "Cannot find sequence $genome_id from genome\n";
   next;
  }
  my @gene_ids = sort @{ $genome_id_to_gene_list_href->{$genome_id} };
  print "\nProcessing scaffold $genome_id\n" if $verbose;
  foreach my $gene_id (@gene_ids) {
   print "\rprocessing gene $gene_id"
     . "                                                  "
     if $verbose;
   my (%params,%preferences);
   $preferences{'sequence_ref'} = \$genome_seq;
   $params{unspliced_transcript} = 1;    # highlights introns

   my $gene_obj_ref = $gene_obj_indexer->get_gene($gene_id);

   $gene_obj_ref->create_all_sequence_types( \$genome_seq, %params );
   my $gene_seq = $gene_obj_ref->get_gene_sequence();
   $gene_seq =~ s/(\S{80})/$1\n/g;
   chomp $gene_seq;
   print GENE ">$gene_id type:gene\n$gene_seq\n";

   foreach my $isoform ( $gene_obj_ref, $gene_obj_ref->get_additional_isoforms() )   {
    my $isoform_id  = $isoform->{Model_feat_name};
    print "\rprocessing gene $gene_id isoform $isoform_id                                              " if $verbose;
    next unless $isoform->has_CDS() || !$isoform->get_CDS_span();
    my @model_span  = $isoform->get_CDS_span();
    next if ( abs( $model_span[0] - $model_span[1] ) < 3 );

    my $common_name = $isoform->{transcript_name} || $isoform->{com_name};
    my $description = '';
    my $alt_name    = '';
    my $main_id     = $isoform_id;
    
    if ( $common_name && $change_name ) {
     if ( $common_name =~ /\s/ ) {
      $description = $common_name;
      $description =~ s/^\s*(\S+)\s*//;
      $description = uri_escape($description);
      $common_name = $1 || die;
     }

     if ( $common_name =~ /\.\d+$/ && !$lettername){ 
      $main_id = $common_name;
     }elsif ($lettername && $common_name =~ /-R[A-Z]+$/ ) {
      $main_id = $common_name;
     }elsif ($lettername){
       $main_id = $common_name . '-RA';
     }else{
       $main_id = $common_name . '.1';
     }
     $alt_name = "($isoform_id) ";
     if ( $unique_names_check{$common_name} ) {
      if ($lettername && $common_name =~ /-R[A-Z]+$/ ) {
       die "Common name $common_name ends in transcript notation but it is not unique!\n";
      }
      elsif (!$lettername && $common_name =~ /\.\d+$/) {
       die "Common name $common_name ends in transcript notation but it is not unique!\n";
      }
      if ($lettername){
       my $letter = 'B';
       for (my $i=1;$i<$unique_names_check{$common_name};$i++){
        $letter++;
       }
       $main_id  = $common_name . '-R'.$letter;
      }else{
       $main_id  = $common_name . '.'.($unique_names_check{$common_name}+1);
      }
      $alt_name = "($isoform_id) ";
     }
     $unique_names_check{$common_name}++;
     # set description as note and update name
     $isoform->{transcript_name} =$main_id;
     $isoform->{com_name} = $main_id;
     $isoform->{pub_comment} = $description if $description;
    }

    # get sequences
    # CDS
    my $seq = $isoform->get_CDS_sequence();
    $seq =~ s/(\S{80})/$1\n/g if $seq;
    chomp $seq if $seq;
    if ( $seq && length($seq) >= $minorf ) {
     print CDS
       ">$main_id ".$alt_name."type:CDS  gene:$gene_id$description\n$seq\n";
    }
    else {
     warn "Transcript $main_id has no coding sequence. Will not process and will skip it\n";
     next;
    }

    # proteins
    $seq = $isoform->get_protein_sequence();
    $seq =~ s/(\S{80})/$1\n/g;
    chomp $seq;
    print PEP
      ">$main_id ".$alt_name."type:polypeptide  gene:$gene_id$description\n$seq\n";

    # mRNA (all exons)
    $seq = $isoform->get_cDNA_sequence();
    $seq =~ s/(\S{80})/$1\n/g;
    chomp $seq;
    print MRNA
      ">$main_id ".$alt_name."type:mRNA gene:$gene_id$description\n$seq\n";


    eval { $isoform->set_CDS_phases( \$genome_seq ); };

   }

   # GFF3
   print GFF3 $gene_obj_ref->to_GFF3_format_extended(%preferences) . "\n";

  }
 }
 print "\nDone!\n";
 close GFF3;
 close PEP;
 close CDS;
 close MRNA;
 close GENE;

 # rename( $gff3_file, "$gff3_file.original" ) unless -s "$gff3_file.original";
 # unlink $index_file;
 # rename( "$gff3_file.gff3", $gff3_file );
 &sort_gff3("$gff3_file.gff3");
 return;
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

sub sort_gff3() {

 # find out why two empty lines still cause sort to crash
 my $gff       = shift;
 my $delimiter = shift;
 $delimiter = &get_gff_delimiter($gff) if !$delimiter;
 my $orig_sep = $/;
 open( GFF, $gff ) || confess "Cannot open $gff $!";
 open( GFF3, ">$gff.sorted" );
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
  print GFF3 $rec . $/;
 }

 close GFF3;
 $/ = $orig_sep;
 rename( "$gff.sorted", $gff );
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
