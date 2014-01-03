#!/usr/bin/env perl

use strict;
use warnings;

my $usage = "usage: $0 CGskew_output\n\n";

my $infile = $ARGV[0] or die $usage;

my @pos_to_sums;
my $num_seqs = 0;

open (my $fh, $infile) or die "Error, cannot open file $infile";
while (<$fh>) {
	chomp;
	my ($acc, @pts) = split (/\s+/);
	for (my $i = 0; $i <= $#pts; $i++) {
		$pos_to_sums[$i] += $pts[$i];
	}
	$num_seqs++;
}
close $fh;


print "Avg";
foreach my $sum (@pos_to_sums) {
	printf("\t%.4f", $sum/$num_seqs);
}
print "\n";

exit(0);

