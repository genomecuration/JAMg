#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use GFF3_utils;
use Gene_obj;

my $usage = "usage: $0 annotation.gff3\n\n";

my $annot_gff3 = $ARGV[0] or die $usage;



main: {
	
	my $gene_obj_indexer_href = {};
	
	my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($annot_gff3, $gene_obj_indexer_href);

	foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
		
		my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};

		foreach my $gene_id (@gene_ids) {

			my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};

			my $orient = $gene_obj_ref->get_orientation();

			my @introns = $gene_obj_ref->get_intron_coordinates();
			
			my $gene_id = $gene_obj_ref->{TU_feat_name};
			my $com_name = $gene_obj_ref->{com_name};
			my $pub_locus = $gene_obj_ref->{pub_locus};
			
			my ($gene_lend, $gene_rend) = sort {$a<=>$b} $gene_obj_ref->get_gene_span();
			my ($model_lend, $model_rend) = sort {$a<=>$b} $gene_obj_ref->get_model_span();

			my @block_starts;
			my @block_lengths;

			foreach my $exon (sort {$a->{end5}<=>$b->{end5}} $gene_obj_ref->get_exons()) {
				if (my $cds = $exon->get_CDS_exon_obj()) {
					my ($lend, $rend) = sort {$a<=>$b} $cds->get_coords();
					
					my $block_len = $rend - $lend + 1;
					push (@block_lengths, $block_len);
					
					my $block_start = $lend - $gene_lend;
					push (@block_starts, $block_start);
				}
			}
			
			my $name_info = $pub_locus . " " . $com_name;
			$name_info =~ s/\s+/_/g;

			print join("\t", $asmbl_id,
					   $gene_lend, $model_rend,
					   $name_info,
					   1,
					   $orient,
					   $model_lend, 
					   $model_rend,
					   "",
					   scalar(@block_lengths),
					   join(",", @block_lengths),
					   join(",", @block_starts)) . "\n";
		}
	}


	exit(0);
	
}



