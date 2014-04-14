#!/usr/bin/env perl

=pod

=head1 NAME

 create_features_from_gff3.pl

=head1 USAGE

 Provide a GFF3 file and a genome FASTA file to phase and create sequence feature (CDS,PEP,MRNA,GENE) information for
 
 create_features_from_gff3.pl genes.gff3 genome.fasta

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

my $minorf                   = 3;     #minimum orf size in bp

my $gfffile = shift || die ("Provide a GFF3 file and a genome FASTA file to phase and create sequence feature (CDS,PEP,MRNAS) information for\n");
my $genome = shift  || die ("Provide a GFF3 file and a genome FASTA file to phase and create sequence feature (CDS,PEP,MRNAS) information for\n");

die "GFF file $gfffile does not exist\n" unless -s $gfffile;
die "Fasta file $genome does not exist\n" unless -s $genome;

my $scaffold_seq_hashref = &read_fasta($genome);

my %unique_names_check;

&gff3_fix_phase($gfffile);

sub gff3_fix_phase() {
 my $gff3_file = shift;
 open( IN, $gff3_file ) || confess( "Cannot find $gff3_file " . $! );
 my $index_file = "$gff3_file.inx";
 my $gene_obj_indexer = new Gene_obj_indexer( { "create" => $index_file } );
 my $asmbl_id_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs( $gff3_file, $gene_obj_indexer );
 open( OUT,  ">$gff3_file.gff3" );
 open( PEP,  ">$gff3_file.pep" );
 open( CDS,  ">$gff3_file.cds" );
 open( GENE, ">$gff3_file.gene" );
 open( MRNA, ">$gff3_file.mRNA" );

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

   foreach my $isoform ( $gene_obj_ref, $gene_obj_ref->get_additional_isoforms() )   {
    $isoform->delete_isoforms();
    my $isoform_id = $isoform->{Model_feat_name};
    my $common_name = $isoform->{com_name};
    my @model_span = $isoform->get_CDS_span();
    next
      if ( !$isoform->get_CDS_span()
           || abs( $model_span[0] - $model_span[1] ) < 3 );

    eval { $isoform->set_CDS_phases( \$genome_seq ); };

    # check name
    my ($original_id,$alt_name,$main_id,$gene_name) = ($isoform_id,$common_name,$isoform_id,$gene_id);
    my $description='';
    if ($common_name){
	if ($common_name=~/\s/){
		$description = $common_name;
		$description =~s/^\s*(\S+)//;
		$common_name=$1||die;
	}
	if (!$unique_names_check{$common_name}){
		if ($common_name=~/\.\d{1,2}$/){
			$main_id = $common_name;
		}else{
			$main_id = $common_name.'.1';
		}
		$alt_name = $isoform_id;
	}else{
		if ($common_name=~/\.\d{1,2}$/){
			die "Common name $common_name ends in transcript notation (.[digit]) but it is not unique!\n";
		}
		my $counter = $unique_names_check{$common_name} + 1;
		$main_id = $common_name.".$counter";
		$alt_name = $isoform_id;
 	}
	$unique_names_check{$common_name}++;
    }


    # get sequences
    # CDS
    my $seq = $isoform->get_CDS_sequence();
    $seq =~ s/(\S{80})/$1\n/g if $seq;
    chomp $seq if $seq;
    if ( $seq && length($seq) >= $minorf ) {
     print CDS ">$main_id ($alt_name) gene:$gene_name$description\n$seq\n";
    }
    else {
     warn "Gene $gene_name has no coding sequence. Will not process and will skip it\n";
     next;
    }

    # proteins
    $seq = $isoform->get_protein_sequence();
    $seq =~ s/(\S{80})/$1\n/g;
    chomp $seq;
    print PEP ">$main_id ($alt_name) gene:$gene_name$description\n$seq\n";

    # gene (+introns)
    $seq = $isoform->get_gene_sequence();
    $seq =~ s/(\S{80})/$1\n/g;
    chomp $seq;
    print GENE ">$main_id ($alt_name) gene:$gene_name$description\n$seq\n";

    # mRNA (all exons)
    $seq = $isoform->get_cDNA_sequence();
    $seq =~ s/(\S{80})/$1\n/g;
    chomp $seq;
    print MRNA ">$main_id ($alt_name) gene:$gene_name$description\n$seq\n";

    # GFF3
    print OUT $isoform->to_GFF3_format_extended(%preferences) . "\n";

   }
  }
 }
 close OUT;
 close PEP;
 close CDS;
 close MRNA;
 close GENE;
 rename( $gff3_file, "$gff3_file.original" ) unless -s "$gff3_file.original";
 unlink $index_file;
 rename( "$gff3_file.gff3", $gff3_file );
 &sort_gff3("$gff3_file");
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
