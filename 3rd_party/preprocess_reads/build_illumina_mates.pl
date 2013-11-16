#!/usr/bin/env perl

=pod

=head1 NAME

 build_illumina_mates.pl

=head1 USAGE

    'i:s'          => Input FASTQ files (expects 2)
    'length_min:i' => Dump in a rejected seqs file those sequences which are shorter than this length, regardless of mate presence
    'fasta'        => Input is a fasta file. Otherwise it is expecting a fastq
    'outfile'      => optionally give name for combined outfile
    'comb'         => Produce a combined pair file (interleaved)
    'sorted'       => Files are sorted by read ID
    'tosort'       => Sort files before processing
    'maxmemory'    => Do not use more memory than this many kb (for sorting). Defaults to 20971520
    cpu:i           => How many CPUs to use (def 1)

=head1 AUTHORS

 Alexie Papanicolaou 1  

        Ecosystem Sciences, CSIRO, Black Mountain Labs, Clunies Ross Str, Canberra, Australia
        alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

This software is released under the GNU General Public License version 3 (GPLv3).
It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html.
Please note that incorporating the whole software or parts of its code in proprietary software
is prohibited under the current license.

=head1 BUGS & LIMITATIONS

None known so far.

=cut

use strict;
use warnings;
use Time::Progress;
use Getopt::Long;
use Pod::Usage;

my $cwd = `pwd`;
chomp($cwd);
my $debug =0;
my $cpu =1;
$| = 1;
my (
     $infile1,  $infile2,          $min_length,
     $is_fasta, $outfile_combined, $comb,
     $sorted,   $tosort,           @infiles
);
my $maxmemory = 8971520;    # kb.8gb
GetOptions(
            'i:s{2}'       => \@infiles,
            'length_min:i' => \$min_length,
            'fasta'        => \$is_fasta,
            'out:s'        => \$outfile_combined,
            'comb'         => \$comb,
            'sorted'       => \$sorted,
            'tosort'       => \$tosort,
            'maxmemory:i'  => \$maxmemory,
	    'cpu:i' => \$cpu
);
pod2usage "Files must be sorted or -tosort given prior processing\n"
  unless $sorted || $tosort;
@infiles = @ARGV if !@infiles;
$infile1 = $infiles[0];
$infile2 = $infiles[1];

unless (    $infile1
         && $infile2
         && -s $infile1
         && -s $infile2
         && $infile1 ne $infile2 )
{
 warn "Need two fastq files!\n";
 pod2usage;
}

my $sort_exec = `which sort` || die "Cannot find sort";
chomp($sort_exec);
die "Cannot find sort command\n" unless $sort_exec && -s $sort_exec;
my $sort_version = `$sort_exec --version|head -n 1`;
$sort_version =~/(\d+)\.\d+$/;
$sort_version = $1;
$sort_exec .= " --parallel=$cpu " if $sort_version && $sort_version >= 8 && $cpu > 1;

my $fasta_check = `head -n1 $infile1|grep "^>"` if !$is_fasta;
if ( $fasta_check || $is_fasta ) {
 print "File $infile1 is probably a FASTA\n";
 $is_fasta = 1;
}
else {
 print "File $infile1 is a FASTQ\n";
}
$fasta_check = `head -n1 $infile2|grep "^>"` if !$is_fasta;
if ( $is_fasta || $fasta_check ) {
 print "File $infile2 is probably a FASTA\n";
 $is_fasta++;
}
else {
 print "File $infile2 is a FASTQ\n";
}
die "Both files must be FASTA or FASTQ\n"
  unless ( !$is_fasta || $is_fasta == 2 );
my $head_sign = $is_fasta ? '>' : '@';
my $max_seq;
if ( -s $infile1 >= -s $infile2 ) {
 $max_seq = -s $infile1;
}
else {
 $max_seq = -s $infile2;
 my $tmpfile = $infile2;
 $infile2 = $infile1;
 $infile1 = $tmpfile;
}
my $timer = new Time::Progress;
$outfile_combined = $infile1 . ".mated.combined.fastq" if !$outfile_combined;
my ( $total, $length, $bytes1b, $bytes2b ) = ( 0, 0, 0, 0 );
my ( $counter, %hash );
if ($sorted) {
 ### cat $1 |awk '{OFS="\t"; getline seqs; getline sep; getline quals; print 0$0,seqs,sep,quals}'|'. $sort_exec.' -k1,1|awk '{OFS="\n";sub(/^0/,"");print $1,$2,"+",$4}'
 &process_sorted( $infile1, $infile2 );
}
elsif ( $tosort && $infile1 !~ /.sorted$/ && $infile2 !~ /.sorted$/ ) {
 print "Sorting sequences...\n";
 print "File $infile1\n";
 my $buffer = int( ( ( -s $infile1 ) + ( ( -s $infile1 ) * 0.5 ) ) / 1024 );
 $buffer = ( $buffer < $maxmemory ) ? $buffer : $maxmemory;
 system( "cat $infile1 "
  . '|awk \'{OFS="\t"; getline seqs; getline sep; getline quals; print 0$0,seqs,sep,quals}\'| '.$sort_exec.' -S'
  . $buffer
  . ' -k1,1|awk \'{OFS="\n";sub(/^0/,"");print $1,$2,"+",$4}\'  >'
  . " $infile1.sorted" )
   if !-s "$infile1.sorted" && !$is_fasta;
 system(   "cat $infile1 "
         . '|awk \'{OFS="\t"; getline seqs; print 0$0,seqs}\'| '.$sort_exec.' -S'
         . $buffer
         . ' -k1,1|awk \'{OFS="\n";sub(/^0/,"");print $1,$2}\'  >'
         . " $infile1.sorted" )
   if !-s "$infile1.sorted" && $is_fasta;
 #sleep(3);
 print "File $infile2\n";
 $buffer = int( ( ( -s $infile2 ) + ( ( -s $infile2 ) * 0.5 ) ) / 1024 );
 $buffer = ( $buffer < $maxmemory ) ? $buffer : $maxmemory;
 system( "cat $infile2 "
  . '|awk \'{OFS="\t"; getline seqs; getline sep; getline quals; print 0$0,seqs,sep,quals}\'| '.$sort_exec.' -S'
  . $buffer
  . ' -k1,1|awk \'{OFS="\n";sub(/^0/,"");print $1,$2,"+",$4}\'  >'
  . " $infile2.sorted" )
   if !-s "$infile2.sorted" && !$is_fasta;
 system(   "cat $infile2 "
         . '|awk \'{OFS="\t"; getline seqs; print 0$0,seqs}\'| '.$sort_exec.' -S'
         . $buffer
         . ' -k1,1|awk \'{OFS="\n";sub(/^0/,"");print $1,$2}\'  >'
         . " $infile2.sorted" )
   if !-s "$infile2.sorted" && $is_fasta;
 my $sinfile1 = $infile1 . ".sorted";
 my $sinfile2 = $infile2 . ".sorted";
 &process_sorted( $sinfile1, $sinfile2 );
 unlink($sinfile1);
 unlink($sinfile2);
}
elsif ( $tosort && $infile1 =~ /.sorted$/ && $infile2 =~ /.sorted$/ ) {
 print "Assuming files are already sorted\n";
 &process_sorted( $infile1, $infile2 );
}
else {
 %hash = &build_hash($infile2);
 $timer->restart;
 &print_hashed_data($infile1);
}

print "Finished with $total mated sequences of a total length of $length b.p.\n";

#########################
sub process_sorted($$) {
 print "Starting sorted comparison...\n";
 $timer->attr( min => 0, max => $max_seq );
 $timer->restart;
 my $counter = 0;
 open( OUT, ">$infile1.mated.combined" ) if $comb;
 open( OUTS1, ">$infile1.mated.singletons.fastq" );
 open( OUTS2, ">$infile2.mated.singletons.fastq" );
 open( OUT1,  ">$infile1.mated.fastq" );
 open( OUT2,  ">$infile2.mated.fastq" );
 open( IN1,   $infile1 ) || die;
 open( IN2,   $infile2 ) || die;

 while ( my $sid1 = <IN1> ) {
  chomp($sid1);
  $sid1 =~ s/^[>@](\S+)\s*.*$/$1/;
  die "1:Id line does not look like fastq/a\n$sid1\n" . <IN1> . <IN1> if !$1;
  "a"=~/a/;
  $sid1 =~ s/(\/\d+)$//;
  my $suffix1 = $1 ? $1 : '';
  my $seq1 = <IN1>;
  my ( $qid1, $qlt1 );
  $qid1 = <IN1> if !$is_fasta;
  $qlt1 = <IN1> if !$is_fasta;

  my $sid2 = <IN2>;
  my ( $seq2, $qid2, $qlt2,$suffix2 );
  if ($sid2) {
   chomp($sid2);
   $sid2 =~ s/^[>@](\S+)\s*.*$/$1/;
   die "2:Id line does not look like fastq/a\n$sid2\n" . <IN2> . <IN2> if !$1;
   "a"=~/a/;
   $sid2 =~ s/(\/\d+)$//;
   $suffix2 = $1 ? $1 : '';
   $seq2 = <IN2>;
   $qid2 = <IN2> if !$is_fasta;
   $qlt2 = <IN2> if !$is_fasta;
  }

  if ($is_fasta) {
   &decide_sort_print( $sid1, $seq1, $suffix1, $sid2, $seq2, $suffix2 );
   $counter += length($sid1) + length($seq1);
  }
  else {
   &decide_sort_print( $sid1, $seq1, $suffix1, $qid1, $qlt1, $sid2, $seq2, $suffix2, $qid2, $qlt2 );
   $counter += length($sid1) + length($seq1) + length($qid1) + length($qlt1);
  }
  if ( $counter =~ /00000$/ ) {
   print $timer->report( "eta: %E min, %40b %p\r", $counter );
  }
 }
 while ( my $sid2 = <IN2> ) {
  my $seq2 = <IN2>;
  my ( $qid2, $qlt2 );
  $qid2 = <IN2> if !$is_fasta;
  $qlt2 = <IN2> if !$is_fasta;
  print OUTS2 "$sid2\n" . $seq2 . "+\n$qlt2" if !$is_fasta;
  print OUTS2 "$sid2\n" . $seq2 if $is_fasta;
 }
 close OUTS1;
 close OUTS2;
 close OUT1;
 close OUT2;
 close OUT if $comb;
 close IN1;
 close IN2;
}

sub decide_sort_print() {
 no warnings 'recursion';
 my $sid1 = shift;
 my $seq1 = shift;
 my $suffix1 = shift;
 $suffix1 = '' if !$suffix1;
 my $qid1 = shift if !$is_fasta;
 my $qlt1 = shift if !$is_fasta;
 my $sid2 = shift;
 my $seq2 = shift;
 my $suffix2 = shift;
 $suffix2 = '' if !$suffix2;

 my $qid2 = shift if !$is_fasta;
 my $qlt2 = shift if !$is_fasta;
 if ( !$sid1 && $sid2 ) {
  print OUTS2 $head_sign . $sid2 . "$suffix2\n" . $seq2 . "+\n" . $qlt2 if !$is_fasta;
  print OUTS2 $head_sign . $sid2 . "$suffix2\n" . $seq2 if $is_fasta;
  return;
 }
 elsif ( !$sid2 && $sid1 ) {
  print OUTS1 $head_sign . $sid1 . "$suffix1\n" . $seq1 . "+\n" . $qlt1 if !$is_fasta;
  print OUTS1 $head_sign . $sid1 . "$suffix1\n" . $seq1 if $is_fasta;
  return;
 }
 elsif ( !$sid1 && !$sid2 ) {
  return;
 }
 if ( $sid1 eq $sid2 ) {
  $total++;
  $length += length($seq1) + length($seq2);
  print OUT1 $head_sign . $sid1 . "$suffix1\n" . $seq1 . "+\n" . $qlt1 if !$is_fasta;
  print OUT1 $head_sign . $sid1 . "$suffix1\n" . $seq1                 if $is_fasta;
  print OUT2 $head_sign . $sid2 . "$suffix2\n" . $seq2 . "+\n" . $qlt2 if !$is_fasta;
  print OUT2 $head_sign . $sid2 . "$suffix2\n" . $seq2                 if $is_fasta;
  if ($comb) {
   print OUT $head_sign . $sid1 . "$suffix1\n" . $seq1 . "+\n" . $qlt1 if !$is_fasta;
   print OUT $head_sign . $sid1 . "$suffix1\n" . $seq1                 if $is_fasta;
   print OUT $head_sign . $sid2 . "$suffix2\n" . $seq2 . "+\n" . $qlt2 if !$is_fasta;
   print OUT $head_sign . $sid2 . "$suffix2\n" . $seq2                 if $is_fasta;
  }
 }
 elsif ( $sid1 ge $sid2 ) {
  print OUTS2 $head_sign . $sid2 . "$suffix2\n" . $seq2 . "+\n" . $qlt2 if !$is_fasta;
  print OUTS2 $head_sign . $sid2 . "$suffix2\n" . $seq2 if $is_fasta;
  $sid2 = <IN2>;
  chomp($sid2);
  $sid2 =~ s/^[>@](\S+)\s*.*$/$1/;
  "a"=~/a/;
  $sid2 =~ s/(\/\d+)$//;
  $suffix2 = $1;
  $seq2 = <IN2>;
  $qid2 = <IN2> if !$is_fasta;
  $qlt2 = <IN2> if !$is_fasta;
  &decide_sort_print( $sid1, $seq1, $suffix1, $sid2, $seq2,$suffix2 ) if $is_fasta;
  &decide_sort_print( $sid1, $seq1, $suffix1, $qid1, $qlt1, $sid2, $seq2, $suffix2, $qid2, $qlt2 )
    if !$is_fasta;
 }
 else {
  print OUTS1 $head_sign . $sid1 . "$suffix1\n" . $seq1 . "+\n" . $qlt1 if !$is_fasta;
  print OUTS1 $head_sign . $sid1 . "$suffix1\n" . $seq1 if $is_fasta;
  $sid1 = <IN1>;
  if ($sid1) {
   chomp($sid1);
   $sid1 =~ s/^[>@](\S+)\s*.*$/$1/;
   "a"=~/a/;
   $sid1 =~ s/(\/\d+)$//;
   $suffix1 = $1;
  }
  $seq1 = <IN1>;
  $qid1 = <IN1> if !$is_fasta;
  $qlt1 = <IN1> if !$is_fasta;
  &decide_sort_print( $sid1, $seq1, $suffix1, $sid2, $seq2, $suffix2 ) if $is_fasta;
  &decide_sort_print( $sid1, $seq1, $suffix1, $qid1, $qlt1, $sid2, $seq2, $suffix2, $qid2, $qlt2 )
    if !$is_fasta;
 }
}

sub build_hash($) {
 my $file_to_hash = shift;
 open( IN2, $file_to_hash ) || die;
 print "Building HASH\n";
 $timer->attr( min => 0, max => -s $file_to_hash );
 $timer->restart;
 while ( my $seqdesc = <IN2> ) {
  my $id = $seqdesc;
  $id =~ s/^(\S+)\s*.*$/$1/;
  "a"=~/a/;
  $id =~ s/\d+$//;
  my $seq = <IN2>;
  $seq .= "\n" if ( $is_fasta && $seq =~ /.+\z$/ );
  my $ql = '';
  if ($is_fasta) {
   $hash{$id}{'data'} = $seqdesc . $seq if $is_fasta;
  }
  else {
   $ql = <IN2> . <IN2>;
   $ql .= "\n" if ( $ql =~ /.+\z$/ );
   $hash{$id}{'data'} = $seqdesc . $seq . $ql;
  }
  $hash{$id}{'length'} = length($seq);
  $counter += length($seqdesc) + length($seq) + length($ql);
  if ( $counter =~ /00000$/ ) {
   print $timer->report( "eta: %E min, %40b %p\r", $counter );
  }
 }
 close IN2;
}

sub print_hashed_data($) {
 my $file_to_print = shift;

 print "Will also produce a combined file $infile1.mated.combined\n" if $comb;
 open( OUT, ">$infile1.mated.combined" ) if $comb;
 open( OUT1B, ">$infile1.mated.singletons" );
 if ($min_length) {
  open( OUT1C, ">$infile1.short.$min_length" );
  open( OUT2C, ">$infile2.short.$min_length" );
  open( OUT1,  ">$infile1.long.$min_length" );
  open( OUT2,  ">$infile2.long.$min_length" );
 }
 else {
  open( OUT1, ">$infile1.mated" );
  open( OUT2, ">$infile2.mated" );
 }
 open( IN1, $file_to_print ) || die;
 $counter = 0;
 print "\nPrinting data...\n";
 while ( my $seqdesc = <IN1> ) {
  $counter += 2;
  if ( $counter =~ /00000$/ ) {
   print $timer->report( "eta: %E min, %40b %p\r", $counter );
  }
  my $id = $seqdesc;
  $id =~ s/^(\S+)\s*.*$/$1/;
  "a"=~/a/;
  $id =~ s/\d+$//;
  my $seq = <IN1>;
  $seq .= "\n" if ( $is_fasta && $seq =~ /.+\z$/ );
  my $ql = '';
  unless ($is_fasta) {
   $ql = <IN1> . <IN1>;
   $ql .= "\n" if ( $ql =~ /.+\z$/ );
  }
  if ( $hash{$id}{'data'} ) {

   # too small
   if ( $min_length
       && ( $hash{$id}{'length'} < $min_length || length($seq) < $min_length ) )
   {
    print OUT2C $hash{$id}{'data'};
    print OUT1C $seqdesc . $seq . $ql if !$is_fasta;
    print OUT1C $seqdesc . $seq if $is_fasta;
   }

   # ok
   else {
    $total++;
    $length += $hash{$id}{'length'} + length($seq);
    print OUT $seqdesc . $seq . $ql . $hash{$id}{'data'} if !$is_fasta && $comb;
    print OUT $seqdesc . $seq . $hash{$id}{'data'} if $is_fasta && $comb;
    print OUT2 $hash{$id}{'data'};
    print OUT1 $seqdesc . $seq . $ql if !$is_fasta;
    print OUT1 $seqdesc . $seq if $is_fasta;
   }
   delete( $hash{$id} );
  }

  # singleton
  else {
   print OUT1B $seqdesc . $seq . $ql if !$is_fasta;
   print OUT1B $seqdesc . $seq if $is_fasta;
  }
 }
 print "\nPrinting more singletons...\n";
 $timer->attr( min => 0, max => scalar( keys %hash ) );
 $timer->restart;
 $counter = 0;
 foreach my $id ( keys %hash ) {
  $counter++;
  next if !$id;
  if ( $counter =~ /00000$/ ) {
   print $timer->report( "eta: %E min, %40b %p\r", $counter );
  }

  # singleton
  print OUT1B $hash{$id}{'data'} if $hash{$id}{'data'};
 }
 close IN1;
 close OUT1;
 close OUT2;
 close OUT1B;
 close OUT1C;
 close OUT2C;
 close OUT if $comb;
}

