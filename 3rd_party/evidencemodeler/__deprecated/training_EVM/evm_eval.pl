#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use GTF_utils;
use Gene_obj_indexer;
use Data::Dumper;
use SnSp::SnSp_analysis_manager;

my $usage = "\nusage: $0 reference.GTF predictions.GTF\n\n";
my $reference_GTF_file = $ARGV[0] or die $usage;
my $predictions_GTF_file = $ARGV[1] or die $usage;


my $FLANK_REGIONS_SnSp_BP = 500;

main: {
    
    my $ref_inx_file = "tmp.$$.ref.inx";
    my $pred_inx_file = "tmp.$$.pred.inx";

    my $ref_gene_obj_indexer = new Gene_obj_indexer( { create => $ref_inx_file } );
    my $ref_seqname_map_href = &GTF_utils::index_GTF_gene_objs($reference_GTF_file, $ref_gene_obj_indexer);
    

    my $pred_gene_obj_indexer = new Gene_obj_indexer( { create => $pred_inx_file } );
    my $pred_seqname_map_href = &GTF_utils::index_GTF_gene_objs($predictions_GTF_file, $pred_gene_obj_indexer);


    ## examine SnSp for each contig included in the reference set:
    my $snsp_analyzer = new SnSp::SnSp_analysis_manager($predictions_GTF_file);
    $snsp_analyzer->set_intergenic_included($FLANK_REGIONS_SnSp_BP);
    
    foreach my $contig (keys %$ref_seqname_map_href) {
        ## get the reference genes:
        
        my $ref_geneids_aref = $ref_seqname_map_href->{$contig};

        my $pred_geneids_aref = $pred_seqname_map_href->{$contig};
        
        foreach my $ref_gene_id (@$ref_geneids_aref) {
            my $ref_gene_obj = $ref_gene_obj_indexer->get_gene($ref_gene_id);
            
            my ($ref_lend, $ref_rend) = sort {$a<=>$b} $ref_gene_obj->get_model_span();
            
            print "Ref gene, $contig, $ref_gene_id\t$ref_lend-$ref_rend\n";


            my @pred_genes;
            ## get the pred genes within range:
            foreach my $pred_gene_id (@$pred_geneids_aref) {
                my $pred_obj = $pred_gene_obj_indexer->get_gene($pred_gene_id);
                
                my ($pred_lend, $pred_rend) = sort {$a<=>$b} $pred_obj->get_model_span();
                
                if ($pred_lend < $ref_rend && $pred_rend > $ref_lend) {
                    
                    print "\tPred gene: $pred_gene_id\t$pred_lend-$pred_rend\n";
                    
                    push (@pred_genes, $pred_obj);
                    
                }
            }
            
            my %ev_type_to_pred_obj = ( $predictions_GTF_file => \@pred_genes );
            
            $snsp_analyzer->add_analysis_entry($ref_gene_obj, \%ev_type_to_pred_obj, 1);
        }
    }
    
    unlink ($ref_inx_file, $pred_inx_file);

    
}
