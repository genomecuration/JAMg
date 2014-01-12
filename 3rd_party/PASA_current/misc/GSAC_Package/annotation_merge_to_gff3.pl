#!/usr/local/bin/perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;
use CdbTools;
use GFF3_utils;
use Carp;

$|++;

my $usage = "\n\nusage: $0 merge_output genedb.inx\n\n";

my $merge_output = $ARGV[0] or die $usage;
my $inx_file = $ARGV[1] or die $usage;

my $gene_obj_indexer = new Gene_obj_indexer( { "use" => "$inx_file" } );

open (my $fh, "$merge_output") or die $!;
while (<$fh>) {
    chomp;
    my @x = split (/\t/);
    my $selected_models = $x[3];
    
    unless ($selected_models) {
        die "Error, no selected models! $_";
    }
    
    my @models = split (/,/, $selected_models);
    foreach my $model (@models) {
        my ($gene_acc, $model_acc, $org) = split (/\:+/, $model);
        
        my $gene_obj = $gene_obj_indexer->get_gene($model_acc);
        
        print $gene_obj->to_GFF3_format() . "\n";
    }
}

exit(0);


    
