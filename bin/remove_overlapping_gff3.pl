#!/usr/bin/env perl

#NB: This was coded in a rush, it's crap but had to get something done

use strict;
use warnings;
use Carp;
use Data::Dumper;

my $gff_file = shift || die "Provide GFF and Priority Source ".$!;
my $priority_source = shift || die $!;

&remove_overlapping_gff($gff_file,$gff_file.'.purged',$priority_source);


#########################################

sub remove_overlapping_gff() {

 # we allow genes on opposite strands
 my $file               = shift;
 my $out                = shift;
 my $source_keep = shift;
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
  my $source = $gene_data[1];

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
    if($master_gene_list{$mrna_id}{'source'} eq $source_keep ){
	   print LOG2 "Gene $mrna_id found more than once. Skipping new data and keeping what I found first\n";
	   $skipped++;
	   next GENE;
	# we assume source_keep has no overlapping genes we want to keep!
    }
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
       && $smallest_coord <= $last_position_check{$ref_id}{$strand}
       && $strand eq $master_gene_list{$previous_mrna_id}{'strand'}
	 )
  {
	   print LOG2 "Gene overlapping $previous_mrna_id found: $mrna_id\n";
	   $overlap_found{$previous_mrna_id}{$mrna_id}{'data'} = $record . $delimiter;
	   $overlap_found{$previous_mrna_id}{$mrna_id}{'source'} = $source;
	   
  }

  $master_gene_list{$mrna_id}{'data'}    = $record . $delimiter;
  $master_gene_list{$mrna_id}{'source'}   = $source;
  $master_gene_list{$mrna_id}{'strand'}   = $strand;
  $previous_mrna_id                      = $mrna_id;
  $previous_ref                          = $ref_id;
  $last_position_check{$ref_id}{$strand} = $largest_coord;
 }
 $/ = $orig_sep;
 close IN;
 print LOG2 "\nDeciding on overlaps:\n";

 my %kept;
 my %skipped_hash;

 open( OUT, ">$out" );
 foreach my $gene ( sort keys %overlap_found){
  if ( $overlap_found{$gene} ) {
   if ($master_gene_list{$gene}{'source'} eq $source_keep){
    foreach my $id (keys %{$overlap_found{$gene}}){
	    $skipped_hash{$id}=1 unless $id eq $gene;
    }
   }
  }
 }

 foreach my $gene ( sort keys %master_gene_list ) {
  next if $skipped_hash{$gene};

  if ( $overlap_found{$gene} ) {
   my $master_size = $cds_sizes{$gene};
   my $gene_longest;
   my $right_source;
   if ($master_gene_list{$gene}{'source'} eq $source_keep){
    print LOG2 "1 $gene from $source_keep found: Keeping only $gene\n";
    print OUT $master_gene_list{$gene}{'data'};
    my @lines = split( "\n",  $master_gene_list{$gene}{'data'} );
    my @data  = split( "\t", $lines[0] );
    $kept{$gene} = $data[1];
    $skipped += scalar( keys %{ $overlap_found{$gene} } );
    foreach my $id (keys %{$overlap_found{$gene}}){
	    $skipped_hash{$id}=1 unless $id eq $gene;
    }
    next;	
   }

   foreach my $overlap_gene ( keys %{ $overlap_found{$gene} } ) {
    my $size = $cds_sizes{$overlap_gene};
    $gene_longest = $overlap_gene if $size > $master_size;
    my $this_source = $overlap_found{$gene}{$overlap_gene}{'source'};
    $right_source = $overlap_gene if $this_source eq $source_keep;
   }

   if (!$right_source){
	print LOG2 "No genes with source $source_keep found. Keeping all\n";
        print OUT $master_gene_list{$gene}{'data'};
        my @lines = split( "\n", $master_gene_list{$gene}{'data'} );
        my @data  = split( "\t", $lines[0] );
        $kept{$gene} = $data[1];
   }
   elsif($right_source){
    print LOG2 "2 $gene from $source_keep found: Keeping only $right_source\n";
#    print OUT $overlap_found{$gene}{$right_source}{'data'};
    
    my @lines = split( "\n", $overlap_found{$gene}{$right_source} );
    my @data  = split( "\t", $lines[0] );
    $kept{$right_source} = $data[1];
    $skipped += scalar( keys %{ $overlap_found{$gene} } );
    $skipped_hash{$gene} = 1 unless $gene eq $right_source;
    foreach my $id (keys %{$overlap_found{$gene}}){
	    $skipped_hash{$id}=1 unless $id eq $right_source;
    }
    
   }
  }
  else {
   # no overlaps
   print OUT $master_gene_list{$gene}{'data'};
   my @lines = split( "\n", $master_gene_list{$gene}{'data'} );
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
 
 print "\tOverlaps checked. Kept ".scalar(keys %kept)." genes. See $out.log for details\n";
# return $skipped;
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
