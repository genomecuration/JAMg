#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use GFF3_utils;
use Carp;

$|++;

my $usage = "\n\nusage: $0 file.inx\n\n";

my $inx_file = $ARGV[0] or die $usage;

unless (-s $inx_file) {
    die $usage;
}

my $gene_obj_indexer = new Gene_obj_indexer( { "use" => $inx_file } );

my @gene_ids = $gene_obj_indexer->get_keys();
foreach my $gene_id (@gene_ids) {
    
    

    my $gene_obj = $gene_obj_indexer->get_gene($gene_id);
    
    my $contig_id = $gene_obj->{asmbl_id};
    
    my $strand = $gene_obj->get_orientation();
    my $model_id = $gene_obj->{Model_feat_name};
    my $source = $gene_obj->{source} || ".";
    
    my @exons = $gene_obj->get_exons();
    
	my @coord_features;
	foreach my $exon (@exons) {
		my ($lend, $rend) = sort {$a<=>$b} $exon->get_coords();
		push (@coord_features, [$lend, $rend, "exon"]);
	}

	my @UTR_coords = ($gene_obj->get_5prime_UTR_coords(), $gene_obj->get_3prime_UTR_coords());
	foreach my $utr_seg (@UTR_coords) {
		my ($lend, $rend) = sort {$a<=>$b} @$utr_seg;
		push (@coord_features, [$lend, $rend, "utr"]);
	}


	my $orient_char = ($strand eq '+') ? ">" : "<";

	my ($gene_lend, $gene_rend) = sort {$a<=>$b} $gene_obj->get_coords();

	@coord_features = sort {$a->[0]<=>$b->[0]} @coord_features;
	print "$orient_char $gene_lend $gene_rend $model_id\n";
	
	foreach my $coord_feature (@coord_features) {
		print join (" ", @$coord_feature) . "\n";
	}
	print "\n"; # spacer
}


exit(0);



