#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

use lib ("$FindBin::RealBin/../../PerlLib");
use Gene_obj;

my $usage = "usage: $0 alignments.gff3\n\n";

my $trans_gff3_file = $ARGV[0] or die $usage;


main: {
	
	my %genome_trans_to_coords;
	
	open (my $fh, $trans_gff3_file) or die "Error, cannot open file $trans_gff3_file";
	while (<$fh>) {
		chomp;
		
		unless (/\w/) { next; }
		
		my @x = split(/\t/);

		unless (scalar (@x) >= 8 && $x[8] =~ /ID=/) {
			print STDERR "ignoring line: $_\n";
			next;
		}
		
		my $scaff = $x[0];
		my $type = $x[2];
		my $lend = $x[3];
		my $rend = $x[4];

		my $orient = $x[6];
		
		my $info = $x[8];
		
		my @parts = split(/;/, $info);
		my %atts;
		foreach my $part (@parts) {
			$part =~ s/^\s+|\s+$//;
			$part =~ s/\"//g;
			my ($att, $val) = split(/=/, $part);
			
			if (exists $atts{$att}) {
				die "Error, already defined attribute $att in $_";
			}
			
			$atts{$att} = $val;
		}

		my $gene_id = $atts{ID} or die "Error, no gene_id at $_";
		my $trans_id = $atts{Target} or die "Error, no trans_id at $_";
		{
			my @pieces = split(/\s+/, $trans_id);
			$trans_id = shift @pieces;
		}
		
		my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

		$genome_trans_to_coords{$scaff}->{$gene_id}->{$trans_id}->{$end5} = $end3;

	}


	## Output genes in gff3 format:

	foreach my $scaff (sort keys %genome_trans_to_coords) {

		my $genes_href = $genome_trans_to_coords{$scaff};

		foreach my $gene_id (keys %$genes_href) {

			my $trans_href = $genes_href->{$gene_id};

			foreach my $trans_id (keys %$trans_href) {

				my $coords_href = $trans_href->{$trans_id};

				my $gene_obj = new Gene_obj();

				$gene_obj->{TU_feat_name} = $gene_id;
				$gene_obj->{Model_feat_name} = $trans_id;
				$gene_obj->{com_name} = "$gene_id $trans_id";
				
				$gene_obj->{asmbl_id} = $scaff;
				
				$gene_obj->populate_gene_object($coords_href, $coords_href);
			
				print $gene_obj->to_BED_format();
								
			}
		}
	}


	exit(0);
}

