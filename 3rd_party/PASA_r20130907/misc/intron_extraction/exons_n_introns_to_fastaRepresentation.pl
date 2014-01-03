#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use CdbTools;
use Fasta_reader;

my $usage = "usage: $0 parsed_introns collapsed_exons genome_fasta\n\n";

################################
## Coding output here:
#       U: Unknown
#       E: Exon
#       I: Intron
################################


my $parsed_introns = $ARGV[0] or die $usage;
my $collapsed_exons = $ARGV[1] or die $usage;
my $genome_fasta_file = $ARGV[2] or die $usage;

my %contig_intron_coords;
my %contig_exon_regions;

main: {
	
	my %all_contigs;

	{ # parse the introns
		open (my $fh, $parsed_introns) or die "Error, cannot open file $parsed_introns";
		while (<$fh>) {
			chomp;
			my ($contig_id, $orient, $seq_range, $len, $seq) = split (/\t/);
			my ($lend, $rend) = split (/-/, $seq_range);

			push (@{$contig_intron_coords{$contig_id}}, [$lend, $rend]);
		}
		close $fh;
		
	}

	{ # parse the exons
		open (my $fh, $collapsed_exons) or die "Error, cannot open file $collapsed_exons";
		while (<$fh>) {
			chomp;
			my ($contig, $coord_range) = split (/\t/);
			my ($lend, $rend) = split (/-/, $coord_range);
			push (@{$contig_exon_regions{$contig}}, [$lend, $rend]);
		}
		close $fh;
	}
	
	my $fasta_reader = new Fasta_reader($genome_fasta_file);
	
	while (my $seq_obj = $fasta_reader->next()) {
		
		my $accession = $seq_obj->get_accession();
		my $sequence = $seq_obj->get_sequence();
		my $seq_len = length($sequence);
		
		my @chars;
		for (1..$seq_len) {
			push (@chars, 'U');
		}
		
		#$sequence =~ s/(\S{60})/$1\n/g;
		#print ">$accession\n$sequence\n";
		

		## add in the exons:
		my $exon_regions_aref = $contig_exon_regions{$accession};
		if (ref $exon_regions_aref) {
			foreach my $exon_coordset (@$exon_regions_aref) {
				my ($lend, $rend) = sort {$a<=>$b} @$exon_coordset;
				for my $i ($lend .. $rend) {
					$chars[$i-1] = 'E';
				}
			}
		}
		
		my $introns_aref = $contig_intron_coords{$accession};
		if (ref $introns_aref) {
			my @introns = @$introns_aref;
			foreach my $intron (@introns) {
				my ($lend, $rend) = sort {$a<=>$b} @$intron;
				for my $i ($lend .. $rend) {
					$chars[$i-1] = 'I';
				}
			}
		}

		my $coding_string = join ("", @chars);
		$coding_string =~ s/(\S{60})/$1\n/g;
		chomp $coding_string;
		
		print ">$accession\n$coding_string\n";
		
	}
	
	exit(0);


}




