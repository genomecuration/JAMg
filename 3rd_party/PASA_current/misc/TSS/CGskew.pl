#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;

my $usage = "usage: $0 fastaFileSameLengths winLength\n\n";

my $fasta = $ARGV[0] or die $usage;
my $winLength = $ARGV[1] or die $usage;



main: {

	my $seqLength;

	my $fasta_reader = new Fasta_reader($fasta);
	while (my $seq_obj = $fasta_reader->next()) {
		my $sequence = $seq_obj->get_sequence();
		my $acc = $seq_obj->get_accession();
		
		## ensure constant seqLength:
		if ($seqLength && length($sequence) != $seqLength) {
			die "Error, inconsistent sequence lengths";
		}
		elsif (! $seqLength) {
			$seqLength = length($sequence);
		}
		
		my @CGskew;
		
		my @chars = split (//, $sequence);
		for (my $i = 0; $i <= length($sequence) - $winLength; $i++) {
			my $left_win_pos = $i;
			my $right_win_pos = $i + $winLength - 1;
			my @win_chars = @chars[$left_win_pos..$right_win_pos];
			my @G = grep { /G/i } @win_chars;
			my @C = grep { /C/i } @win_chars;

			my $numG = scalar(@G);
			my $numC = scalar(@C);
			
			my $total = $numG + $numC;
			
			my $skew = ($total > 0) ? ( ($numC-$numG)/$total) : 0;
			
			$CGskew[$i] = sprintf("%.4f", $skew);
		}

		print "$acc\t" . join (" ", @CGskew) . "\n";
	}
	
	
	exit(0);
}


