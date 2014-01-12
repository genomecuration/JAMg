#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use CdbTools;
use Nuc_translator;

my $usage = "usage: $0 TSS_info genome.db flank_5prime flank_3prime\n\n";

my $TSS_info = $ARGV[0] or die $usage;
my $genomeDB = $ARGV[1] or die $usage;
my $flank_5prime = $ARGV[2] or die $usage;
my $flank_3prime = $ARGV[3] or die $usage;

main: {
	my %contig_to_TSS;
	
	open (my $fh, $TSS_info) or die "Error, cannot open file $TSS_info";
	while (<$fh>) {
		chomp;
		my ($contig, $gene, $len, $orient, $lend, $rend) = split (/\t/);
		
		my $struct = { gene => $gene,
					   orient => $orient,
					   lend => $lend,
					   rend => $rend,
		};
		push (@{$contig_to_TSS{$contig}}, $struct);
	}
	close $fh;


	foreach my $contig (sort keys %contig_to_TSS) {
		
		my $genome = &cdbyank_linear($contig, $genomeDB);

		my $genome_length =length($genome);
		
		foreach my $struct (@{$contig_to_TSS{$contig}}) {
			my ($gene, $orient, $lend, $rend) = ($struct->{gene},
												 $struct->{orient},
												 $struct->{lend},
												 $struct->{rend});
			
			my $region_seq = "";
			my $TSS;

			if ($orient eq '+') {
				$TSS = $lend;
				my $left_pos = $TSS - $flank_5prime;  ## include TSS base on the 3' side.
				my $right_pos = $TSS + $flank_3prime -1;
			
				## TSS 5prime side coordinates: (TSS - flank_5prime, TSS -1)
				## TSS 3prime side coordinates: (TSS, TSS + flank_3prime -1)
	
				if ($left_pos < 1 || $right_pos > $genome_length) { next; } # ignore
				
				$region_seq = substr($genome, $left_pos-1, $right_pos - $left_pos + 1);
			}
			else {
				$TSS = $rend;
				my $left_pos = $TSS - $flank_3prime + 1;
				my $right_pos = $TSS + $flank_5prime;
				
				## TSS 5prime side coordinates: (TSS + flank_5prime, TSS + 1)
				## TSS 3prime sice coordinates: (TSS - flank_3prime + 1, TSS)



				if ($left_pos < 1 || $right_pos > $genome_length) { next; } # ignore
				
				$region_seq = substr($genome, $left_pos-1, $right_pos - $left_pos + 1);
				$region_seq = &reverse_complement($region_seq);
			}

			$region_seq = uc $region_seq;
			my @chars = split (//, $region_seq);
			for my $i (0..$flank_5prime-1) {
				$chars[$i] = lc $chars[$i];
			}
			$region_seq = join ("", @chars);
			
			print ">$gene $contig $orient $TSS\n$region_seq\n";
		}
	}


	exit(0);
	
}
				
