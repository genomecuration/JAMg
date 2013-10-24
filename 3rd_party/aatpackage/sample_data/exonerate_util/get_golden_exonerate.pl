#!/usr/bin/env perl

=pod

=head1 NAME

 get_golden_exonerate.pl

=head1 USAGE

	'exonerate:s' => Exonerate results file with special format (see perldoc)
	'genome:s'   => Genome sequence FASTA file
	'cdna'     => Exonerate was run with cDNA otherwise expect protein (no start/stop searched in query, just genome)
	'nosingle' => Do not allow single coding exon genes (Exonerate has bugs relating to these genes)
	'training' => number of genes to go into a random training set (def. 66% of total or 5000, whichever is higher)

=head1 DESCRIPTION

 use the following option in exonerate:

 --ryo "RYOAP_START\nRYOAP_STATS_D\tAlignment length\tidentical\tsimilar\tmismatch\tidentical%%\tsimilar%%\nRYOAP_STATS\t%s\t%et\t%ei\t%es\t%em\t%pi\t%ps\nRYOAP_CODING_QUERY\t%qi\t%qab\t%qae\n>%qi\n%qas\nRYOAP_CODING_GENOME\t%ti\t%tcb\t%tce\n>%ti\n%tcs\nRYOAP_END\n"

=head1 Criteria for pass

 Methionine and * in protein sequence
 Start codon and stop codon in genome
 Fraction of identical residues > 95%
 Fraction of similar residues > 90%
 No more than 20 mismatches allowed
 TODO:
 For a particular genomic start/stop co-ordinate & length, only one gene is printed (one with best score) but isoforms (alt. splicings) are allowed.

=cut

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use Bio::Index::Fasta;
use Digest::SHA qw/sha1_hex/;
use Data::Dumper;
my ( $exonerate_file, $genome_file, $is_cdna, $no_single_exon );
my $training_set_size = 5000;
GetOptions(
            'exonerate:s' => \$exonerate_file,
            'genome:s'    => \$genome_file,
            'cdna'        => \$is_cdna,
            'nosingle'    => \$no_single_exon,
            'training:i'       => \$training_set_size,
);
$exonerate_file = shift if !$exonerate_file;


my $shuf_exec = `which shuf 2>/dev/null`;chomp($shuf_exec); #coreutils

#pod2usage unless $exonerate_file && -s $exonerate_file;
# IF WE NEED TO GET WHOLE SCAFFOLD SEQUENCE:
pod2usage
  unless $exonerate_file
   && $genome_file
   && -s $exonerate_file
   && -s $genome_file;
unless ( -f "$genome_file.index" ) {
 print "Indexing $genome_file...\n";
 my $inx = Bio::Index::Fasta->new( -filename   => "$genome_file.index",
                                   -write_flag => 1 );
 $inx->make_index($genome_file);
}
my $genome_inx = Bio::Index::Fasta->new( -filename => "$genome_file.index" )
  || die("Could not get index for $genome_file\n");
## --ryo "RYOAP_START\nRYOAP_STATS_D\tAlignment length\tidentical\tsimilar\tmismatch\tidentical%%\tsimilar%%\nRYOAP_STATS\t%s\t%et\t%ei\t%es\t%em\t%pi\t%ps\nRYOAP_CODING_QUERY\t%qi\t%qab\t%qae\n>%qi\n%qas\nRYOAP_CODING_GENOME\t%ti\t%tcb\t%tce\n>%ti\n%tcs\nRYOAP_END\n"
print "Processing $exonerate_file\n";
open( EXONERATE, $exonerate_file );
open( NEWGFF,    ">$exonerate_file.corrected.gff3" );
open( GENEIDGFF, ">$exonerate_file.corrected.geneid.gff3" );
open( GLIMMER, ">$exonerate_file.corrected.glimmer.txt" );
my ( %details, $flag, $scounter, $ecounter );
my $gene_flag;
# EXONERATE has no UTR, the gene is actually the coding part of the mRNA
while ( my $ln = <EXONERATE> ) {
 next if $ln =~ /^\s*$/;
 if ( $ln =~ /^# --- START OF GFF DUMP/ ) {
  my ($gene_id,$gene_start,$gene_end);
  my $exon_counter   = 0;
  my $intron_counter = 1;
  my ($mRNA,$mRNA_id);
  while ( my $ln2 = <EXONERATE> ) {
   last if $ln2 =~ /^# --- END OF GFF DUMP/;
   next if $ln2 =~ /^#/;
   my @data = split( "\t", $ln2 );
   if ( $data[8] ) {
    next if $data[2] eq 'similarity';
    my $offset = int(0);
    if ( $data[0] =~ s/:(\d+)-(\d+)$// ) {
     $offset = $1 - 1;
     $data[3] += $offset;
     $data[4] += $offset;
    }
    if ( $data[2] eq 'gene' ) {
     $data[8] =~ /sequence\s+(\S+)/;
     $gene_id = $1;
     $gene_start = $data[3];
     $gene_end = $data[4];
     undef($mRNA);#reset
     undef($mRNA_id);
     #pasa
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
     else {
      $data[8] = "ID=$gene_id\n";
     }
    }
    elsif ( $data[2] eq 'cds' ) {
     $data[2] = 'CDS';
     $data[8] = "ID=cds.$gene_id;Parent=$gene_id\n";
    }
    elsif ( $data[2] eq 'exon' ) {
     $exon_counter++;
     $data[8] = "ID=$gene_id.exon.$exon_counter;Parent=$gene_id\n";
    }
    elsif ( $data[2] eq 'intron' ) {
     $data[8] = "ID=$gene_id.intron.$intron_counter;Parent=$gene_id\n";
    }
    else {

     # should we checking for splice site here as part of the golden checks?
     next;
    }
    $data[8] =~ s/\s+\;\s+/\;/g;
    if ( $data[2] eq 'gene'){
     print NEWGFF "\n" if $gene_flag;
     print NEWGFF join( "\t", @data );
     print NEWGFF $mRNA if ($mRNA );
     my $geneid_mRNA_id = $mRNA_id ? $mRNA_id : $gene_id;
     print GENEIDGFF "##\n"; # this is used to delimit records so that they can be sorted later. use it for all GFFs? (it is not kept after sorting)
     print GENEIDGFF $data[0] . "\t" . "FORGENEID" . "\t" . 'Sequence' . "\t".$data[3]. "\t".$data[4]. "\t.\t.\t.\t$geneid_mRNA_id\n";
     print GLIMMER "\n" if $gene_flag;
     $gene_flag++;
    }
    else {
     print NEWGFF join( "\t", @data );

     #geneID. glimmer
     if ( $data[2] eq 'exon' ) {
      my $exon_type;
      if ($exon_counter == 1 && ($data[3] == $gene_start && $data[4] == $gene_end) || ($data[4] == $gene_start && $data[3] == $gene_end)){
       $exon_type = "Single";
      }elsif ($exon_counter == 1){
       $exon_type = "First";
      }elsif ($data[6] eq '-' && $data[3] == $gene_start){
       $exon_type = 'Terminal';
      }elsif ($data[6] eq '+' && $data[4] == $gene_end){
       $exon_type = 'Terminal';
      }else{
        $exon_type = 'Internal';
      }
      
      my $geneid_mRNA_id = $mRNA_id ? $mRNA_id : $gene_id;
      print GENEIDGFF $data[0] . "\t" . "FORGENEID" . "\t" . $exon_type . "\t".$data[3]. "\t".$data[4]. "\t.\t".$data[6]."\t.\t$geneid_mRNA_id\n";
      
      
      print GLIMMER $data[0]."\t$data[3]\t$data[4]\n" if $data[6] eq '+';
      print GLIMMER $data[0]."\t$data[4]\t$data[3]\n" if $data[6] eq '-';
      
      
     }
    }
   }
  }

 }
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

  #		$query_fasta .= $query_seq;
  #TARGET
  my $target_data_str = $ln;
  chomp($target_data_str);
  my @target_data = split( "\t", $target_data_str );
  my $target_id = $target_data[1];
  my $target_start =
    $target_data[2] <= $target_data[3] ? $target_data[2] : $target_data[3];
  my $target_end =
    $target_data[3] >= $target_data[2] ? $target_data[3] : $target_data[2];
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

  #		$target_fasta .= $target_seq;
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
close(EXONERATE);
close (GENEIDGFF);
close GLIMMER;
close NEWGFF;
die "Report seems incomplete $scounter != $ecounter!\n"
  unless $scounter && $ecounter && ( $scounter == $ecounter );
my $query_outfile =
  $is_cdna ? "$exonerate_file.passed.cDNA" : "$exonerate_file.passed.protein";
open( OUT,  ">$exonerate_file.data" );
open( LOG,  ">$exonerate_file.golden.log" );
open( PASS, ">$exonerate_file.passed" );
#print PASS "#Query_ID\tSHAHEX of genome-based ORF\n";
open( PASSP,  ">$query_outfile" );
open( PASSG,  ">$exonerate_file.passed.genome" );
open( PASSGC, ">$exonerate_file.passed.genome.cds" );
print OUT
"#Hit number\tQuery length\tTarget length\tScore\tNumber of mismatches\tFraction of identical\tFraction of similar\tquery name\tquery start\tquery end\tquery sequence\tAlignment direction\ttarget name\ttarget start\ttarget end\ttarget sequence\n";
my $pass_counter = int(0);
my %pass_check;
my %accepted;
my $region_exists = int(0);
my $failed_cutoff = int(0);

my $counter = int(0);
$| = 1;
foreach my $query_id (
                       sort { $details{$b}{'score'} <=> $details{$a}{'score'} }
                       keys %details
  )
{
 $counter++;
 print OUT "Hit "
   . $details{$query_id}{'counter'} . "\t"
   . $details{$query_id}{'qlength'} . "\t"
   . $details{$query_id}{'tlength'} . "\t"
   . $details{$query_id}{'score'} . "\t"
   . $details{$query_id}{'mismatch_number'} . "\t"
   . $details{$query_id}{'identical_frac'} . "\t"
   . $details{$query_id}{'similar_frac'} . "\t"
   . $query_id . "\t"
   . $details{$query_id}{'query_start'} . "\t"
   . $details{$query_id}{'query_end'} . "\t"
   . $details{$query_id}{'query_seq'} . "\t"
   . $details{$query_id}{'alignment_strand'} . "\t"
   . $details{$query_id}{'target_id'} . "\t"
   . $details{$query_id}{'target_start'} . "\t"
   . $details{$query_id}{'target_end'} . "\t"
   . $details{$query_id}{'target_seq'} . "\n";
 if (
      $pass_check{
           's'
         . $details{$query_id}{'target_start'} . 'e'
         . $details{$query_id}{'target_end'}
      }
   )
 {
  print LOG
    "$query_id rejected: higher scoring alignment with same co-ordinates ("
    . 's'
    . $details{$query_id}{'target_start'} . 'e'
    . $details{$query_id}{'target_end'} . "):"
    . $pass_check{ 's'
     . $details{$query_id}{'target_start'} . 'e'
     . $details{$query_id}{'target_end'} }
    . "\n";
  $region_exists++;
  next;
 }
 unless ( $details{$query_id}{'identical_frac'} >= 90 ) {
  print LOG "$query_id rejected: identical_frac less than 90: "
    . $details{$query_id}{'identical_frac'} . "\n";
  $failed_cutoff++;
  next;
 }
 unless ( $details{$query_id}{'similar_frac'} >= 95 ) {
  print LOG "$query_id rejected: similar_frac less than 95: "
    . $details{$query_id}{'similar_frac'} . "\n";
  $failed_cutoff++;
  next;
 }
 unless ( $details{$query_id}{'mismatch_number'} <= 20 ) {
  print LOG "$query_id rejected: mismatches more than 20: "
    . $details{$query_id}{'mismatch_number'} . "\n";
  $failed_cutoff++;
  next;
 }
 unless ( $is_cdna || substr( $details{$query_id}{'query_seq'}, 0, 1 ) eq 'M' )
 {
  print LOG "$query_id rejected: query does not start with M: "
    . substr( $details{$query_id}{'query_seq'}, 0, 1 ) . "\n";
  $failed_cutoff++;
  next;
 }
 unless ( $is_cdna
          || substr( $details{$query_id}{'query_seq'}, -1, 1 ) eq '*' )
 {
  print LOG "$query_id rejected: query does not end with stop codon (*): "
    . substr( $details{$query_id}{'query_seq'}, -1, 1 ) . "\n";
  $failed_cutoff++;
  next;
 }
 unless ( substr( $details{$query_id}{'target_seq'}, 0, 3 ) eq 'ATG' ) {
  print LOG "$query_id rejected: target does not start with ATG: "
    . substr( $details{$query_id}{'target_seq'}, 0, 3 ) . "\n";
  $failed_cutoff++;
  next;
 }
 unless (
          (
               substr( $details{$query_id}{'target_seq'}, -3, 3 ) eq 'TAG'
            || substr( $details{$query_id}{'target_seq'}, -3, 3 ) eq 'TAA'
            || substr( $details{$query_id}{'target_seq'}, -3, 3 ) eq 'TGA'
          )
   )
 {
  print LOG
"$query_id rejected: target does not end with stop codon (TAG, TAA or TGA): "
    . substr( $details{$query_id}{'target_seq'}, -3, 3 ) . "\n";
  $failed_cutoff++;
  next;
 }
 my $scaffold_obj = $genome_inx->fetch( $details{$query_id}{'target_id'} );
 my $scaffold_seq;
 my $scaffold_length = $scaffold_obj->length();

 #bug check
 unless (    int( $details{$query_id}{'target_start'} )
          && int( $details{$query_id}{'target_end'} )
          && $details{$query_id}{'target_start'} > 0
          && $details{$query_id}{'target_end'} > 0
          && $details{$query_id}{'target_end'} <= $scaffold_length )
 {
  print LOG "$query_id had to be skipped: has issues with target end ("
    . $details{$query_id}{'target_id'}
    . ") being negative or larger than scaffold size - possible corruption of exonerate file! ("
    . $details{$query_id}{'target_end'}
    . " is <0 or > "
    . $scaffold_length . ")\n";
  next;
 }

 #consider: adding stop codon as it is not reported at end??
 if ( $details{$query_id}{'alignment_strand'} == -1 ) {
  $scaffold_seq = &revcomp(
                            $scaffold_obj->subseq(
                                        $details{$query_id}{'target_start'} - 2,
                                        $details{$query_id}{'target_end'}
                            )
  );
 }
 else {
  $scaffold_seq =
    $scaffold_obj->subseq( $details{$query_id}{'target_start'} + 1,
                           $details{$query_id}{'target_end'} + 3 );
 }
 if ( $no_single_exon
      && length($scaffold_seq) <= length( $details{$query_id}{'target_seq'} ) )
 {
  $failed_cutoff++;
  print LOG "$query_id rejected: single coding exon\n";
  next;
 }
 print LOG "$query_id accepted\n";
 $accepted{$query_id} = 1;
 $pass_check{ 's'
    . $details{$query_id}{'target_start'} . 'e'
    . $details{$query_id}{'target_end'} } = $query_id;
 print PASSP ">$query_id\n" . $details{$query_id}{'query_seq'} . "\n";
 print PASSG ">$query_id "
   . $details{$query_id}{'target_id'}
   . " genome sequence\n$scaffold_seq\n";
 print PASSGC ">$query_id "
   . $details{$query_id}{'target_id'}
   . " genome CDS\n"
   . $details{$query_id}{'target_seq'} . "\n";
 my $scaffold_seq_sha = sha1_hex($scaffold_seq);
 print PASS "$query_id\t$scaffold_seq_sha\n";
 $pass_counter++;
}
close OUT;
close PASS;
close PASSP;
close PASSG;
my $number_of_passing_genes = scalar(keys %accepted);
$training_set_size = int($number_of_passing_genes * 0.66) if $training_set_size > int($number_of_passing_genes * 0.66);   
my %training_genes;
if ($shuf_exec){
 system("$shuf_exec -n $training_set_size $exonerate_file.passed > $exonerate_file.passed.shuffled.trainset");
 open(IN,"$exonerate_file.passed.shuffled.trainset");
 while (my $ln=<IN>){
  $ln=~/^(\S+)/;
  $training_genes{$1}=1;
 }
 close IN;
}else{
 system("head -n $training_set_size $exonerate_file.passed > $exonerate_file.passed.trainset");
 open(IN,"$exonerate_file.passed.trainset");
 while (my $ln=<IN>){
  $ln=~/^(\S+)/;
  $training_genes{$1}=1;
 }
 close IN;
}

if (-s "$exonerate_file.corrected.gff3"){
 &filter_gff( "$exonerate_file.corrected.gff3", \%accepted,"$exonerate_file.corrected.golden.gff3" ) ;
 &filter_gff( "$exonerate_file.corrected.gff3", \%training_genes,"$exonerate_file.corrected.train.gff3" ) ;
}
if (-s "$exonerate_file.corrected.geneid.gff3"){
 &sort_gff2_geneid("$exonerate_file.corrected.geneid.gff3");
 &filter_gff2_geneid( "$exonerate_file.corrected.geneid.gff3", \%accepted,"$exonerate_file.corrected.geneid.golden.gff3" );
 &filter_gff2_geneid( "$exonerate_file.corrected.geneid.gff3", \%training_genes,"$exonerate_file.corrected.geneid.train.gff3" );
}
print LOG
"Done! Processed $counter sequences.\nFound $pass_counter sequences passing criteria.\nAnother $region_exists because another sequence had higher score with same starting and ending coding exons\nAnother $failed_cutoff failed cutoffs and requirements\n\n";
close LOG;
print
"Done! Processed $counter sequences.\nFound $pass_counter sequences passing criteria.\nAnother $region_exists because another sequence had higher score with same starting and ending coding exons\nAnother $failed_cutoff failed cutoffs and requirements\n\n";
############################
sub revcomp {
 my $dna     = shift;
 my $revcomp = reverse($dna);
 $revcomp =~ tr/ACGTacgt/TGCAtgca/;
 return $revcomp;
}

sub filter_gff() {
 my $gff      = shift;
 my $hash_ref = shift;
 my $out = shift;
 $out      = "$gff.filtered" if !$out;
 
 open( GFF, $gff ) || die;
 open( OUT, ">$out" );
 while ( my $ln = <GFF> ) {
  my $id;
  next unless $ln =~ /\sgene\s/;
  $ln =~ /ID=([^;]+);?/;
  $id = $1 || die $ln;

  #pasa
  my $id2 = $id;
  if ( $id =~ /^g\.(\d+)$/ ) {
   $id2 = 'm.' . $1;
  }
  if ( $hash_ref->{$id} || $hash_ref->{$id2} ) {
   print OUT $ln;
   while ( my $ln2 = <GFF> ) {
    last if $ln2 =~ /^\s*$/;
    print OUT $ln2;
   }
  }
 }
}

sub filter_gff2_geneid() {
 my $gff      = shift;
 my $hash_ref = shift;
 my $out = shift;
 my $genes_written;
 $out      = "$gff.filtered" if !$out;
 open( GFF, $gff ) || die;
 open( OUT, ">$out" );
 while ( my $ln = <GFF> ) {
  next if $ln=~/^\s*$/;
  $ln=~/(\S+)$/;
  my $id=$1||die;
  if(!$genes_written && $ln=~/\sSequence\s/ && $hash_ref->{$id}){
   print OUT $ln;
   $genes_written++;
  }elsif($genes_written && $ln=~/\sSequence\s/ && $hash_ref->{$id}){
   print OUT '#$'."\n".$ln; 
  }else{
   print OUT $ln if ( $hash_ref->{$id} );
  }
 }
}

sub sort_gff2_geneid(){
 my $gff      = shift;
 my $orig_sep = $/;
 open( GFF, $gff ) || die;
 open( OUT, ">$gff.sorted" ) || die;
 $/ = <GFF>; #this swas put here so that it automatically determines the delimiter
 my @records = <GFF>;
 chomp(@records);
 print OUT sort {
        my @first = split("\n",$a);
        my @second = split("\n",$b);
        my @split1 = split("\t",$first[0]);
        my @split2 = split("\t",$second[0]);
        return $split1[0] cmp $split2[0] || $split1[3] <=> $split2[3];
 }@records;
 close GFF;
 close OUT;
 rename("$gff.sorted",$gff);
 $/ = $orig_sep;

}


