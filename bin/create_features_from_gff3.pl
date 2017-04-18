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
  -name               Use Transcript common name as the main ID in the output.
  -lettername         Use -R[two letters] notation for alternative transcript instead of .[digits]  
  -verbose            Print progress and debug info
  -skip_delete        Skip Status=Delete mRNAs
  -one_isoform        Process only one isoform per gene
  -rename             Rename the IDs with a new JAMg IDs, up to 32 characters (needed for hhblits). Don't use with -name
  -strip_name         Remove Name tag
  -change_source  :s  Change GFF Source to this value
  -simple             Don't add introns and splice sites in GFF
 
NB: -name means that the common name has no spaces and is unique (will be checked). Useful for WebApollo
 
=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Carp;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use GFF3_utils;
use GTF_utils;

$|=1;
our $SEE;
my $minorf = 3;    #minimum orf size in bp
my ($simple_gff3, $gfffile, $genome, $change_name,$lettername,$verbose, $one_iso, $do_rename, $change_source, $strip_name, $split_single );
pod2usage $! unless &GetOptions(
	          'one_isoform'      => \$one_iso,
            'gff|infile:s'     => \$gfffile,
            'genome|fasta:s'   => \$genome,
            'name|change_name' => \$change_name,
            'debug'            => \$SEE,
            'verbose'          => \$verbose,
            'lettername'       => \$lettername,
            'change_source'    => \$change_source,
      	    'rename'	    	   => \$do_rename,
            'strip_name'       => \$strip_name,
	    'simple' => \$simple_gff3,
	    'split_single' => \$split_single,
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
 confess "Cannot find $gff3_file\n" unless -s $gff3_file;
 my $index_file = "$gff3_file.inx";
 my $gene_obj_indexer;
 if (-s $index_file){
	$gene_obj_indexer = new Gene_obj_indexer({"use" => $index_file});
 }else{
	$gene_obj_indexer =  new Gene_obj_indexer( { "create" => $index_file } );
 }
 confess "Cannot index with $index_file\n" unless $gene_obj_indexer && -s $index_file;
 my $genome_id_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs( $gff3_file, $gene_obj_indexer );

 open( IN, $gff3_file ) || confess( "Cannot find $gff3_file " . $! );
 open( GFF3, ">$gff3_file.gff3" );
 open( GFF3_SINGLE, ">$gff3_file.gff3.single" ) if $split_single;
 open( PEP,  ">$gff3_file.pep" );
 open( CDS,  ">$gff3_file.cds" );
 open( GENE, ">$gff3_file.gene" );
 open( MRNA, ">$gff3_file.mRNA" );
 open( TRACK, ">$gff3_file.track"); 
 print "Indexing complete\n";
 my $gene_count = 1;
 foreach my $genome_id (sort {
	($a =~ /(\d*)$/)[0] cmp ($b =~ /(\d*)$/)[0] || $a cmp $b
  } keys %$genome_id_to_gene_list_href ) {

  my $genome_seq = $scaffold_seq_hashref->{$genome_id};
  if ( !$genome_seq ) {
   warn "Cannot find sequence $genome_id from genome\n";
   next;
  }

  my @gene_ids = sort {
    ($a =~ /(\d*)$/)[0] cmp ($b =~ /(\d*)$/)[0] || $a cmp $b
    } @{ $genome_id_to_gene_list_href->{$genome_id} };

  print "\nProcessing scaffold $genome_id\n" if $verbose;

 foreach my $gene_id (@gene_ids) {
    
   my (%params,%preferences,$is_single_exon);
   $preferences{'sequence_ref'} = \$genome_seq;
   $preferences{'source'}  = $change_source if $change_source;
   $params{unspliced_transcript} = 1;    # highlights introns

   my $gene_obj_ref = $gene_obj_indexer->get_gene($gene_id);
   next if ($gene_id=~/temp_model/  
	|| ($gene_obj_ref->{gene_name} && $gene_obj_ref->{gene_name}=~/temp_model/) 
	|| ($gene_obj_ref->{TU_feat_name} && $gene_obj_ref->{TU_feat_name}=~/temp_model/)
	);

   $gene_obj_ref->{TU_feat_name} = $gene_id if !$gene_obj_ref->{TU_feat_name};

   my $jamg_id_gene  = "JAMg_model_".$gene_count;


   $gene_obj_ref->create_all_sequence_types( \$genome_seq, %params );
   my $gene_seq = $gene_obj_ref->get_gene_sequence();
   next if !$gene_seq;
   $gene_seq =~ s/(\S{80})/$1\n/g;
   chomp $gene_seq;

   # sort longest CDS first.
   my (@isoforms,%cdslength_hash);

   foreach my $iso ( $gene_obj_ref, $gene_obj_ref->get_additional_isoforms() ){
	   if ($iso->{Model_feat_name} =~/temp_model/ || $iso->{transcript_name}=~/temp_model/){
      
      next;
     }

	$cdslength_hash{$iso} = length($iso->get_CDS_sequence());
	if ($one_iso){
		$isoforms[0] = $iso if !$isoforms[0] || $cdslength_hash{$iso} > $cdslength_hash{$isoforms[0]};
	}else{
		push(@isoforms,$iso); # all data
	}
   }

   next unless $isoforms[0];

   if ($do_rename){
     print TRACK "GENE\t$gene_id\t$jamg_id_gene\n";
     $gene_id = $jamg_id_gene;
     $gene_obj_ref->{gene_name} = $gene_id;
     $gene_obj_ref->{TU_feat_name} = $gene_id;
   }
   print GENE ">$gene_id type:gene\n$gene_seq\n";
   print "\rprocessing gene $gene_id"
     . "                                                  "
     if $verbose;

   my $isoform_count;
   foreach my $isoform ( @isoforms )   {
    my $isoform_id  = $isoform->{Model_feat_name};

    next unless $isoform->has_CDS() || !$isoform->get_CDS_span();
    my @model_span  = $isoform->get_CDS_span();
    next if ( abs( $model_span[0] - $model_span[1] ) < 3 );
    
   # there is a real issue here. Most times EVM and PASA are ok with phasing
   # but sometimes there are issues.
   # and these issues happen both ways (with and without rephasing). 
   # so i decided to do rephasing always unless it is a single CDS gene on the negative side.
   # unfortunately this has to happen across the entire gene
# this is now solved within the library via set_CDS (basically first CDS is reset to have phase = 0 (frame = 1)
   if ($isoform->{'num_exons'} == 1){
	$is_single_exon++;
   }
    print "\rprocessing gene $gene_id isoform $isoform_id                                       " if $verbose;
    $isoform_count++;
    my $common_name = $isoform->{transcript_name} || $isoform->{com_name};
    my $description = '';
    my $alt_name    = '';
    my $transcript_main_id     = $isoform_id;
    
    if ( $common_name && $change_name ) {
     if ( $common_name =~ /\s/ ) {
      $description = $common_name;
      $description =~ s/^\s*(\S+)\s*//;
      $description = $description;
      $common_name = $1 || die;
     }

     if ( $common_name =~ /\.\d+$/ && !$lettername){ 
      $transcript_main_id = $common_name;
     }elsif ($lettername && $common_name =~ /-R[A-Z]+$/ ) {
      $transcript_main_id = $common_name;
     }elsif ($lettername){
       $transcript_main_id = $common_name . '-RA';
     }else{
       $transcript_main_id = $common_name . '.1';
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
       $transcript_main_id  = $common_name . '-R'.$letter;
      }else{
       $transcript_main_id  = $common_name . '.'.($unique_names_check{$common_name}+1);
      }
      $alt_name = "($isoform_id) ";
     }
     $unique_names_check{$common_name}++;
     print TRACK "TRANSCRIPT\t$isoform_id\t$transcript_main_id\n";
     # set description as note and update name
     $isoform->{transcript_name} = $transcript_main_id;
     $isoform->{com_name} = $transcript_main_id;
     $isoform->{pub_comment} = $description if $description;
    } # end change name
    elsif($do_rename){
      print TRACK "TRANSCRIPT\t$isoform_id\t$jamg_id_gene.$isoform_count\n";
      $transcript_main_id = $jamg_id_gene.'.'.$isoform_count;
      $isoform->{Model_feat_name} = $transcript_main_id;
    }

    if ($strip_name){
      $isoform->{transcript_name} = $transcript_main_id;
      $isoform->{com_name} = $transcript_main_id;
    }

    # get sequences
    # CDS
    my $seq = $isoform->get_CDS_sequence();
    $seq =~ s/(\S{80})/$1\n/g if $seq;
    chomp $seq if $seq;
    if ( $seq && length($seq) >= $minorf ) {
     print CDS
       ">$transcript_main_id ".$alt_name."type:CDS  gene:$gene_id$description\n".uc($seq)."\n";
    }
    else {
     warn "Transcript $transcript_main_id has no coding sequence. Will not process and will skip it\n";
     next;
    }

    # proteins
    $seq = $isoform->get_protein_sequence();
    $seq =~ s/(\S{80})/$1\n/g;
    chomp $seq;


    print PEP
      ">$transcript_main_id ".$alt_name."type:polypeptide  gene:$gene_id$description\n".uc($seq)."\n";

    # mRNA (all exons)
    $seq = $isoform->get_cDNA_sequence();
    $seq =~ s/(\S{80})/$1\n/g;
    chomp $seq;
    print MRNA ">$transcript_main_id ".$alt_name."type:mRNA gene:$gene_id$description\n".uc($seq)."\n";
    die "OK, this is unexpected: there is a stop codon inside the ORF for transcript $transcript_main_id!\n" if $seq=~/\*\S/;
    eval { $isoform->set_CDS_phases( \$genome_seq ); } unless $preferences{'norephase'};
   }
   # GFF3
   if ($simple_gff3){
	if ($is_single_exon && $split_single){
	        print GFF3_SINGLE $gene_obj_ref->to_GFF3_format(%preferences) . "\n";
	}else{
	        print GFF3 $gene_obj_ref->to_GFF3_format(%preferences) . "\n";
	}
   }else{
	if ($is_single_exon && $split_single){
	   	print GFF3_SINGLE $gene_obj_ref->to_GFF3_format_extended(%preferences) . "\n";
	}else{
	   	print GFF3 $gene_obj_ref->to_GFF3_format_extended(%preferences) . "\n";
	}
   }
   $gene_count++;
  }
 }
 close GFF3;
 close GFF3_SINGLE;
 close PEP;
 close CDS;
 close MRNA;
 close GENE;
 close TRACK;

 &sort_gff3("$gff3_file.gff3");
 unlink("$gff3_file.track") if !-s "$gff3_file.track";

 unlink("$gff3_file.gff3.single") unless -s "$gff3_file.gff3.single";
 &sort_gff3("$gff3_file.gff3.single") if -s "$gff3_file.gff3.single";
 # rename( "$gff3_file.gff3", $gff3_file );
 print "\nDone!\n";

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


sub fix_single_exon_phase(){
	my $obj = shift;
	my $strand = $obj->get_orientation();
	 foreach my $exon ( $obj->get_exons() ) {
	  if ( my $cds = $exon->get_CDS_exon_obj() ) {
		my $phase = $cds->get_phase();
		next if !$phase || $phase == 0;
		#delete this many from beginning
		if ($strand eq '+'){
			$exon->{end5}+= $phase;
			$cds->{end5}+= $phase;
			$cds->{phase}= 0;			
		}else{
			# i think
			$exon->{end5}-= $phase;
			$cds->{end5}-= $phase;
			$cds->{phase}= 0;	
		}

	  }
	 }

	return $obj;

}
