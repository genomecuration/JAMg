#!/usr/bin/env perl

use FindBin;
use lib ("$FindBin::Bin/../PerlLib", "$FindBin::Bin/PerlLib");
use Annot_set_comparator;
use strict;
use Gene_obj;
use Gene_obj_indexer;
use GFF3_utils;
use Getopt::Long qw(:config no_ignore_case bundling);

our $SEE;
my ($gff3_A, $token_A, $gff3_B, $token_B, $help, $genome_db, $contig_id);

my $allow_opposite_strand_comparisons = 0;

my $percent_overlap;

&GetOptions ('h' => \$help,
             'v' => \$SEE, 
             'gff3_A=s' => \$gff3_A,
             'gff3_B=s' => \$gff3_B,
             'token_A=s' => \$token_A,
             'token_B=s' => \$token_B,
             'contig_id=s' => \$contig_id,
			 'allow_opposite_strand_comparisons' => \$allow_opposite_strand_comparisons,
			 'percent_overlap=f' => \$percent_overlap,
			 
             );


my $usage = <<_EOUSAGE_;

#########################################################
#
#   Gene Structure Annotation Comparison
#
#########################################################
#   --gff3_A   annotation set A in gff3 format
#   --token_A  name for annotation set A
#   
#   --gff3_B   annotation set B in gff3 format
#   --token_B  name for annotation set B
#
#   --contig_id  limit analysis to a specific contig id (1st field of gff3)
#
#   --allow_opposite_strand_comparisons  (by default, performing strand comparisons separately)
#   --percent_overlap      min requirement for gene mapping (default 20% either gene length)
#
#  -h this help menu
#  -v verbose
#########################################################

_EOUSAGE_

    ;

$|++;



print STDERR "allow_opposite_strand_comparisons: $allow_opposite_strand_comparisons\n";

if ($help) { 
    die $usage;
}
unless ($gff3_A && $gff3_B && $token_A && $token_B) {
    die $usage;
}


if ($percent_overlap) {
	$Annot_set_comparator::MIN_OVERLAP = $percent_overlap;
}

main: {
    
    ## print header:
    print "#GFF3_A: $gff3_A\tTOKEN_A: $token_A\t"
        . "GFF3_B: $gff3_B\tTOKEN_B: $token_B\n";
    

	my @strand_comparisons;
		
	if ($allow_opposite_strand_comparisons) {
		@strand_comparisons = ('?');
	}
	else {
		@strand_comparisons = ('+', '-');
	}
	
	
    ## parse annotations
    my $tmp_obj_file = "gene_objs.$$.inx";
    my $gene_obj_indexer = new Gene_obj_indexer({ create => $tmp_obj_file } );
    
    my $asmbls_to_A_gene_ids_href = &parse_gff3($gff3_A, $gene_obj_indexer);    
    my $asmbls_to_B_gene_ids_href = &parse_gff3($gff3_B, $gene_obj_indexer);
    
    my %asmbls;
    foreach my $asmbl (keys %$asmbls_to_A_gene_ids_href, keys %$asmbls_to_B_gene_ids_href) {
        $asmbls{$asmbl} = 1;
    }

    
    ## Do comparisons
    foreach my $asmbl_id (sort keys %asmbls) {
        print STDERR "Processing contig: $asmbl_id\n";
    
        foreach my $strand (@strand_comparisons) {
            
            my @genes_A = &get_genes_on_asmbl_id($asmbls_to_A_gene_ids_href->{$asmbl_id}, $strand, $gene_obj_indexer);
            
            my @genes_B = &get_genes_on_asmbl_id($asmbls_to_B_gene_ids_href->{$asmbl_id}, $strand, $gene_obj_indexer);
                        
			
			
            &Annot_set_comparator::do_comparison($asmbl_id, $token_A, \@genes_A, $token_B, \@genes_B);
       
		}
        
    }
    
    unlink ($tmp_obj_file);
    
    exit(0);
    
}

####
sub update_pseudogene_spans {
	my @gene_objs = @_;

	foreach my $gene_obj (@gene_objs) {
		if ($gene_obj->is_pseudogene() || $gene_obj->{com_name} =~ /pseudogene/) {
			my @exons = $gene_obj->get_exons();
			foreach my $exon (@exons) {
				my ($end5, $end3) = $exon->get_coords();
				if (my $cds_exon_obj = $exon->get_CDS_exon_obj()) {
					$cds_exon_obj->{end5} = $end5;
					$cds_exon_obj->{end3} = $end3;
				}
				else {
					$exon->set_CDS_exon_obj(new CDS_exon_obj($end5, $end3));
				}
			}
			$gene_obj->refine_gene_object();
			
			print STDERR "updated pseudogene coordinates: " . $gene_obj->toString();
		}
	}

	

	return;
}



####
sub parse_gff3 {
    my ($gff_filename, $gene_obj_indexer) = @_;
     
    return (&GFF3_utils::index_GFF3_gene_objs($gff_filename, $gene_obj_indexer, $contig_id));
    
}


####
sub get_genes_on_asmbl_id {
    my ($gene_id_list_aref, $strand, $gene_obj_indexer) = @_;

    my @gene_objs;
   
    foreach my $gene_id (@$gene_id_list_aref) {
        my $gene_obj = $gene_obj_indexer->get_gene($gene_id);
        my $orient = $gene_obj->get_orientation();
				
        if ($strand eq '?' || $orient eq $strand) {
            push (@gene_objs, $gene_obj);
        }
    }

	&update_pseudogene_spans(@gene_objs); ## pseudogenes should span their entire length. Don't trust any protein-coding CDS annotation here.
	foreach my $gene_obj (@gene_objs) {
		$gene_obj->trim_UTRs();  ## note, UTR regions are being ignored
	}
	
    return (@gene_objs);
}
