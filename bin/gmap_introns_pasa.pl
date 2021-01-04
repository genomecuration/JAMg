#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use GFF3_utils;
use Carp;
use Nuc_translator;


my $usage = "Provide a PASA genes GFF .inx file\n";

my $inx_file = shift or die $usage;

die $usage unless (-s $inx_file);

my $gene_obj_indexer = new Gene_obj_indexer( { "use" => $inx_file } ) || die $!;

my @gene_ids = $gene_obj_indexer->get_keys();

foreach my $gene_id (sort @gene_ids) {
    my $gene_obj = $gene_obj_indexer->get_gene($gene_id);

    my $contig_id = $gene_obj->{asmbl_id};

    my @intron_coords = $gene_obj->get_intron_coordinates();
    my $strand = $gene_obj->get_orientation();
    my $model_id = $gene_obj->{Model_feat_name};

    my $counter = 1;
    foreach my $intron (@intron_coords) {
        my ($intron_lend, $intron_rend) = sort {$a<=>$b} @$intron;
	#>Prudul26A000004T1.intron1 Pd01:22798934..22798574
	my $print =  ">$model_id.intron$counter $contig_id:";
	if ($strand eq '+'){
		$print .= $intron_lend.'..'.$intron_rend;
	}else{
		$print .= $intron_rend.'..'.$intron_lend;
	}
        print $print."\n";
	$counter++;
    }
}
