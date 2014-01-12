#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ($ENV{EUK_MODULES});
use Gene_obj;
use CdbTools;
use GFF3_utils;
use Carp;
use Nuc_translator;

my $usage = "\n\nusage: $0 gff3_file\n\n";

my $gff3_file = $ARGV[0] or die $usage;

## associate gene identifiers with contig id's.
my $gene_obj_indexer_href = {};
my $contig_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs($gff3_file, $gene_obj_indexer_href);


open (my $prime3_utr_fh, ">3primeUTR_lengths.$$.dat");
open(my $prime5_utr_fh, ">5primeUTR_lengths.$$.dat");

foreach my $asmbl_id (sort keys %$contig_to_gene_list_href) {
    
    my @gene_ids = @{$contig_to_gene_list_href->{$asmbl_id}};
    
    foreach my $gene_id (@gene_ids) {
        print "// processing $gene_id\n";
		my $gene_obj_ref = $gene_obj_indexer_href->{$gene_id};
		&analyze_5prime_UTR($gene_obj_ref);
		&analyze_3prime_UTR($gene_obj_ref);
	}
}


exit(0);


####
sub analyze_5prime_UTR {
	my ($gene_obj_ref) = @_;
	
	my $gene_id = $gene_obj_ref->{TU_feat_name};
	my $orient = $gene_obj_ref->get_orientation();
	my $contig = $gene_obj_ref->{asmbl_id};
	
	my @gene_UTRs;
	
	foreach my $model ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
		
		my @prime5_utr = $model->get_5prime_UTR_coords();
		
		if (@prime5_utr) {
			my $utr_len = &sum_UTR(@prime5_utr);
			push (@gene_UTRs, [$model->{Model_feat_name}, \@prime5_utr, $utr_len]);
		}
	}
	
	if (@gene_UTRs) {
		@gene_UTRs = sort {$a->[2]<=>$b->[2]} @gene_UTRs;
		my $longest_UTR = pop @gene_UTRs;
		my $model_id = $longest_UTR->[0];
		my $utr_list_aref = $longest_UTR->[1];
		my $utr_len = $longest_UTR->[2];
		my @utr_span = &get_utr_span(@$utr_list_aref);
		print $prime5_utr_fh "$contig\t$model_id\t$utr_len\t$orient\t" . join ("\t", @utr_span) . "\n";
	}
	
	return;
}


####
sub analyze_3prime_UTR {
	my ($gene_obj_ref) = @_;
	
	my $gene_id = $gene_obj_ref->{TU_feat_name};
	my $orient = $gene_obj_ref->get_orientation();
	my $contig = $gene_obj_ref->{asmbl_id};
	
	my @gene_UTRs;
	
	foreach my $model ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
		
		my @prime3_utr = $model->get_3prime_UTR_coords();
		
		if (@prime3_utr) {
			my $utr_len = &sum_UTR(@prime3_utr);
			push (@gene_UTRs, [$model->{Model_feat_name}, \@prime3_utr, $utr_len]);
		}
	}
	
	if (@gene_UTRs) {
		@gene_UTRs = sort {$a->[2]<=>$b->[2]} @gene_UTRs;
		my $longest_UTR = pop @gene_UTRs;
		my $model_id = $longest_UTR->[0];
		my $utr_list_aref = $longest_UTR->[1];
		my $utr_len = $longest_UTR->[2];
		my @utr_span = &get_utr_span(@$utr_list_aref);
		print $prime3_utr_fh "$contig\t$model_id\t$utr_len\t$orient\t" . join ("\t", @utr_span) . "\n";
	}
	
	return;
}


####
sub sum_UTR {
	my (@utr_coords) = @_;

	my $sum_len = 0;
	foreach my $utr_region (@utr_coords) {
		my ($end5, $end3) = @$utr_region;
		my $len = abs ($end5 - $end3) + 1;
		$sum_len += $len;
	}

	return($sum_len);
}


####
sub get_utr_span {
	my (@utr_segs) = @_;

	my @coords;
	foreach my $utr_seg (@utr_segs) {
		my ($end5, $end3) = @$utr_seg;
		push (@coords, $end5, $end3);
	}

	@coords = sort {$a<=>$b} @coords;

	my $lend = shift @coords;
	my $rend = pop @coords;

	return ($lend, $rend);
}
