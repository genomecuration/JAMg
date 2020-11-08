#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::RealBin/../../PerlLib");
use WigParser;

my $usage = "usage: $0 strand_coverage.wig strand[+-]\n\n";

# need strand value to include in the GFF file.

my $wig_file = $ARGV[0] or die $usage;
my $strand = $ARGV[1] or die $usage;

main: {

	my $scaffold = "";
	my $lend = undef;
	my $rend = undef;

	open (my $fh, $wig_file) or die "Error, cannot open file $wig_file";
	while (<$fh>) {
		chomp;
		if (/variableStep chrom=(\S+)/) {
			if ($lend && $rend) {
				print "$scaffold\tpartition\tregion\t$lend\t$rend\t.\t$strand\t.\t.\n";
				$lend = undef;
				$rend = undef;
			}
			
			$scaffold = $1;
			next;
		}
		
		my @x = split(/\t/);
		my ($pos, $cov) = @x;
		
		if ($cov) {

			if (! defined $lend) {
				$lend = $pos;
				$rend = $pos;
			}
			else {
				$rend = $pos;
			}
		}
		else {
			# no coverage... may have ended a block.
			if ($lend) { 
				# report coverage block:
				my $len = $rend - $lend + 1;
				print "$scaffold\tpartition\tregion\t$lend\t$rend\t.\t$strand\t.\t.\t$len\n";
				
				# reinit
				$lend = undef;
				$rend = undef;
			}
			
		}
	}
	

	# get last block
	if ($lend && $rend) {
		print "$scaffold\tpartition\tregion\t$lend\t$rend\t.\t$strand\t.\t.\n";
	}
	
	exit(0);

}
		
