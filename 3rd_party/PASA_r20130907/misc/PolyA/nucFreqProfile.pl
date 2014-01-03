#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;

my $usage = "usage: $0 fastaFile\n\n";

my $fastaFile = $ARGV[0] or die $usage;

main: {
	my $seqLength;

	my @char_array;

	my $num_seqs = 0;
	
	my $fasta_reader = new Fasta_reader($fastaFile);
	while (my $seq_obj = $fasta_reader->next()) {
		my $sequence = uc $seq_obj->get_sequence();
		
		$num_seqs++;
		
		if ($seqLength && length($sequence) != $seqLength) {
			die "Error, inconsistent sequence lengths: $seqLength";
		}
		elsif (! $seqLength) {
			$seqLength = length($sequence);
		}
		
		my @chars = split (//, $sequence);
		
		for (my $i = 0; $i <= $#chars; $i++) {
			my $char = $chars[$i];
			$char_array[$i]->{$char}++;
		}
	}
	
	## Report findings:
	my @nucs = qw (G A T C);
	
	foreach my $nuc (@nucs) {
		print "$nuc";
		for (my $i = 0; $i <= $#char_array; $i++) {
			my $num_nuc = $char_array[$i]->{$nuc} || 0;
			my $ratio = sprintf("%.3f", $num_nuc / $num_seqs);
			print "\t$ratio";
		}
		print "\n";
	}
	


	exit(0);
}



