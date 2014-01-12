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

my $min_bp_overlap = 40;

&GetOptions ('h' => \$help,
             'v' => \$SEE, 
             'gff3_A=s' => \$gff3_A,
             'gff3_B=s' => \$gff3_B,
			 'contig_id=s' => \$contig_id,
             );


my $usage = <<_EOUSAGE_;


#########################################################
#   --gff3_A   annotation set A in gff3 format
#   
#   --gff3_B   annotation set B in gff3 format
#
#  optional:
#     --contig_id  limit analysis to a specific contig id (1st field of gff3)
#
#  -h this help menu
#  -v verbose
#########################################################

_EOUSAGE_

    ;

$|++;

if ($help) { 
    die $usage;
}
unless ($gff3_A && $gff3_B) {
    die $usage;
}

main: {
	
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
                
		my @genes_A = &get_genes_on_asmbl_id($asmbls_to_A_gene_ids_href->{$asmbl_id},$gene_obj_indexer);
		
		my @genes_B = &get_genes_on_asmbl_id($asmbls_to_B_gene_ids_href->{$asmbl_id}, $gene_obj_indexer);
		
		
		foreach my $gene (@genes_A) {
			my $TU_feat = $gene->{TU_feat_name};
			my $model_feat = $gene->{Model_feat_name};
			
		}
		
		
		&compare_A_to_B(\@genes_A, \@genes_B);
		
		
	}
	
        
    unlink ($tmp_obj_file);
    
    exit(0);
    
}

####
sub parse_gff3 {
    my ($gff_filename, $gene_obj_indexer) = @_;
     
    return (&GFF3_utils::index_GFF3_gene_objs($gff_filename, $gene_obj_indexer, $contig_id));
    
}


####
sub get_genes_on_asmbl_id {
    my ($gene_id_list_aref, $gene_obj_indexer) = @_;

    my @gene_objs;
   
	

    foreach my $gene_id (@$gene_id_list_aref) {
        my $gene_obj = $gene_obj_indexer->get_gene($gene_id);
        my $orient = $gene_obj->get_orientation();
		
		## examine each isoform independently:
		foreach my $isoform ($gene_obj, $gene_obj->get_additional_isoforms()) {
			
			push (@gene_objs, $isoform);
		}
	}
	
	return (@gene_objs);
}


####
sub compare_A_to_B {
	my ($genes_A_href, $genes_B_href) = @_;
	
	foreach my $gene_A (@$genes_A_href) {

		my $asmbl_id_A = $gene_A->{asmbl_id};

		my ($gene_A_lend, $gene_A_rend) = sort {$a<=>$b} $gene_A->get_coords();
		my $gene_A_orient = $gene_A->get_orientation();
		
		my $gene_ID_A = $gene_A->{TU_feat_name};
		my $gene_model_A = $gene_A->{Model_feat_name};
		
		my $comparison_text = "";


		foreach my $gene_B (@$genes_B_href) {
			
			my ($gene_B_lend, $gene_B_rend) = sort {$a<=>$b} $gene_B->get_coords();
			my $gene_B_orient = $gene_B->get_orientation();

			my $gene_ID_B = $gene_B->{TU_feat_name};
			my $gene_model_B = $gene_B->{Model_feat_name};
			

			# check for overlap
			
			if ($gene_A_lend <= $gene_B_rend - $min_bp_overlap && $gene_A_rend >= $gene_B_lend + $min_bp_overlap) {
				# acceptable overlap
				
				my $compare_IDs = "$gene_ID_A\t$gene_model_A\t$gene_ID_B\t$gene_model_B";
				

				# options for comparisons:
				# A and B are the same
				# A and B have overlapping exon
				# A is in B's intron
				# B is in A's intron

				my $compare_token = "";

				if (&identical_exons($gene_A, $gene_B)) {
					$compare_token = "identical";
				}
				
				elsif (&is_intronic_gene($gene_A, $gene_B)) {
					$compare_token = "A_in_INTRON_of_B";
				}
				elsif (&is_intronic_gene($gene_B, $gene_A)) {
					$compare_token = "B_in_INTRON_of_A";
				}
				elsif(&have_exon_overlap($gene_A, $gene_B)) {
					$compare_token = "exon_overlap";
				}
				else {
					die "Error, cannot classify comparison of geneA: " . $gene_A->toString() . " to geneB: " . $gene_B->toString();

				}

				my $strand_text = ($gene_A_orient eq $gene_B_orient) ? "SAME_strand" : "OPPOSITE_strand";
				
				$comparison_text .= "$compare_IDs\t($asmbl_id_A:$gene_A_lend-$gene_A_rend)\t$compare_token\t$strand_text\n";


			}
		} # end of foreach geneB


		if ($comparison_text) {
			print $comparison_text;
		}
		else {
			# must be in intergenic region
			print "$gene_ID_A\t$gene_model_A\tnone\tnone\t($asmbl_id_A:$gene_A_lend-$gene_A_rend)\tin_intergenic_region\n";
		}
	}
	

	return;
}


####
sub identical_exons {
	my ($geneA, $geneB) = @_;

	my %exons;
	
	foreach my $exon ($geneA->get_exons()) {
		my ($end5, $end3) = $exon->get_coords();
		$exons{"$end5;$end3"} = 1;
	}

	foreach my $exon ($geneB->get_exons()) {
		my ($end5, $end3) = $exon->get_coords();
		
		my $exon_token = "$end5;$end3";
		if (exists ($exons{$exon_token})) {
			delete $exons{$exon_token};
		}
		else {
			return (0); # not identical
		}
	}

	if (%exons) {
		return (0); # not identical
	}
	else {
		return (1); # yes, identical
	}

}



####
sub have_exon_overlap {
	my ($geneA, $geneB) = @_;

	my @exonsA = $geneA->get_exons();

	my @exonsB = $geneB->get_exons();

	foreach my $exonA (@exonsA) {

		my ($exonA_lend, $exonA_rend) = sort {$a<=>$b} $exonA->get_coords();

		foreach my $exonB (@exonsB) {
			
			my ($exonB_lend, $exonB_rend) = sort {$a<=>$b} $exonB->get_coords();

			if ($exonA_lend <= $exonB_rend || $exonA_rend >= $exonB_lend) {
				# got overlap
				return (1);
			}
		}
	}

	return (0); # no overlap
}

#### 
sub is_intronic_gene {
	my ($geneA, $geneB) = @_;

	## Checks to see if geneA is in the intron of geneB

	my ($geneA_lend, $geneA_rend) = sort {$a<=>$b} $geneA->get_coords();

	my @geneB_introns = $geneB->get_intron_coordinates();

	foreach my $geneB_intron (@geneB_introns) {
		my ($intron_lend, $intron_rend) = sort {$a<=>$b} @$geneB_intron;

		if ($geneA_lend >= $intron_lend && $geneA_rend <= $intron_rend) {

			return (1); # yes intronic gene
		}
	}

	return (0); # not found within any intron of B
}


	
			
