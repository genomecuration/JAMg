#!/usr/bin/env perl

=pod 

=head1 USAGE

 -in1|1   :s  => File 1 of pair
 -in2|2   :s  => File 2 of pair
 -samples :i  => Number of sequences to sample (give a float of <= 1.0 to get a proportion). Defaults to 500000

 Giving the same number of -sample as sequences available (or -sample 1 ) will just result in same sequences but shuffled

=cut

use strict;
use warnings;
use Pod::Usage;
use Time::Progress;
use Tie::File;
use Getopt::Long;

$|=1;
my ($file1,$file2,$samples);

GetOptions(
	'in1|1:s' => \$file1,
	'in2|2:s' => \$file2,
	'samples:i' => \$samples
);

$samples = 500000 if !$samples;
my (%jackknive);

pod2usage unless $file1 && -s $file1 && $file2 && -s $file2;
my $sequence_counter=int(0);
my $timer   = new Time::Progress;
$timer->attr( min => 0, max => $samples  );

tie my @file1, 'Tie::File', $file1 or die $!;
tie my @file2, 'Tie::File', $file2 or die $!;

print "Calculating...\n";
my $max_lines=`wc -l < $file1`; chomp($max_lines);
die "Number of lines is not divisible by 4... weird ($max_lines)\n" unless $max_lines % 4 == 0;
my $max_sequences = $max_lines / 4;

if ($samples <= 1){
	$samples = int($max_sequences * $samples);
}

warn "You're asking for more samples ($samples) than available sequences ($max_sequences). This will result in a shuffle with same number of sequences\n" if $samples > $max_sequences;
warn "You're asking for same number of samples ($samples) as available sequences ($max_sequences). This will result in a shuffle with same number of sequences\n" if $samples == $max_sequences;

open (OUT1,">$file1.sample.$samples");
open (OUT2,">$file2.sample.$samples");
print "Printing...\n";
OUTER: for (my $i=0;$i<$samples;$i++){
	if ($sequence_counter == $max_sequences){
		warn "\nMaximum sequences ($max_sequences) reached. Stopping\n";
		last;
	}
        my $random_line = rand($max_lines);
	next if ($random_line > ($max_lines+3)) ||  $jackknive{$random_line} || (($random_line % 4) != 0);
        my $ln=$file1[$random_line]."\n";
	my $seq=$file1[$random_line+1]."\n";
	my $id=$file1[$random_line+2]."\n";
	my $qual=$file1[$random_line+3]."\n";
        print OUT1 $ln.$seq.$id.$qual;
        my $ln2=$file2[$random_line]."\n";
	my $seq2=$file2[$random_line+1]."\n";
	my $id2=$file2[$random_line+2]."\n";
	my $qual2=$file2[$random_line+3]."\n";
        print OUT2 $ln2.$seq2.$id2.$qual2;
	$jackknive{$random_line} = 1;
	$sequence_counter++;
	print $timer->report( "eta: %E min, %40b %p\r", $sequence_counter ) if $sequence_counter =~ /000$/;
}
print "\nDone, see $file1.sample.$samples & $file2.sample.$samples\n";

close OUT1;
close OUT2;


