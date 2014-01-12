#!/usr/local/bin/perl

package Annot_set_comparator;

use strict;
use Gene_obj;
use Gene_obj_comparator;
use SingleLinkageClusterer;
require "overlapping_nucs.ph";
use CDNA::Gene_obj_alignment_assembler;
use CDNA::PASA_alignment_assembler;

$|++;

my $SEE = 0;

our $MIN_OVERLAP = 20;

our (%feat_name_to_gene_obj);
my ($SEQNAME, $TYPEA, $TYPEB); # just handy as globals

####
sub do_comparison {
    my ($seqname, $type_A, $gene_objs_list_A_aref, $type_B, $gene_objs_list_B_aref) = @_;
    
    ## =================================================================================
    ## Requirements: each gene_obj must have the TU_feat_name and Model_feat_name attributes populated with the gene_id and model_id, respectively.  
    ## Also, alternative splicing isoforms MUST be included in each gene_obj using the gene_obj->add_isoform() method.
    ## If genes are partial, the is_[53]prime_partial() method must provide accurate info.
    ##
    ## Finally, note that we change the TU and Model_feat_name attributes to be the following: 
    ##      gene_id::model_id:type
    ## =================================================================================
    
    $SEQNAME = $seqname;
    $TYPEA = $type_A;
    $TYPEB = $type_B;
    local %feat_name_to_gene_obj;  # by doing this, the memory will be cleared upon exit of this routine. 
    
    &prepareIDs ($type_A, $gene_objs_list_A_aref);
    &prepareIDs ($type_B, $gene_objs_list_B_aref);
    
    ## do comparison (all-vs-all):
    
    my %feat_names;
    foreach my $gene_obj (@$gene_objs_list_A_aref, @$gene_objs_list_B_aref) {
        my $feat_name = $gene_obj->{TU_feat_name};
        #print "FEAT: $feat_name\n";
        $feat_names{$feat_name} = 1;
        
        $feat_name_to_gene_obj{$feat_name} = $gene_obj;
    
    }
    
    
    my @pairs;
    foreach my $type_B_model (@$gene_objs_list_B_aref) {
        
        # print $type_B_model->toString();
        
        my ($type_B_lend, $type_B_rend) = sort {$a<=>$b} $type_B_model->get_gene_span();
        my $type_B_name = $type_B_model->{TU_feat_name};
        
        my $type_B_length = $type_B_rend - $type_B_lend + 1;
        
        my $B_orient = $type_B_model->get_orientation();


        foreach my $type_A_model (@$gene_objs_list_A_aref) {
            
            #print $type_A_model->toString();
            
            my ($type_A_lend, $type_A_rend) = sort {$a<=>$b} $type_A_model->get_gene_span();
            
            my $type_A_length = $type_A_rend - $type_A_lend + 1;
            
            my $type_A_name = $type_A_model->{TU_feat_name};
            
            my $A_orient = $type_A_model->get_orientation();
            
			#################################################################################
			## Note, restrictions based on same or opposite orientation are now
			## enforced by pre-selection of the genes to be compared.  Not Here!!! 
			##
            # unless ($B_orient eq $A_orient) { next; } ## require the same transcribed strand.
			#
            ####################################################################################
			
            if ($type_B_lend < $type_A_rend && $type_B_rend > $type_A_lend) {
                # got overlap
                
                unless (&have_exon_overlap($type_A_model, $type_B_model)) { 
                    # print "Sorry, no exon overlap\n";
                    next; 
                }
                

                my $nucs_common = &nucs_in_common ($type_B_lend, $type_B_rend, $type_A_lend, $type_A_rend);
                # print "nucs_common: $nucs_common\n";
                
                if ($nucs_common/$type_B_length * 100 > $MIN_OVERLAP ||
                    $nucs_common/$type_A_length * 100 >$MIN_OVERLAP ) {
                    
                    
                    push (@pairs, [$type_B_name, $type_A_name]);
                }
            }
        }
    }
    
    my @clusters = &SingleLinkageClusterer::build_clusters(@pairs);
    
        
    foreach my $cluster (@clusters) {
        
        my $ele_list = join (",", @$cluster);
        
        print "cluster: $ele_list\n" if $SEE;
        
        my @type_B_models;
        my @type_A_models;
        foreach my $ele (@$cluster) {
            if ($ele =~ /$type_A/i) {
                push (@type_A_models, $ele);
            } 
            elsif ($ele =~ /$type_B/i) {
                push (@type_B_models, $ele);
            }
            else {
                die "Error, don't recognize ele: $ele\n";
            }
        }
        
        
        my $num_type_B = scalar (@type_B_models);
        my $num_type_A = scalar (@type_A_models);
        
        
        eval {
            ## Check one-to-one mapping
            if ($num_type_B == 1 && $num_type_A == 1) {
                
                &compare_1_to_1_mapping ($feat_name_to_gene_obj{$type_B_models[0]}, 
                                         $feat_name_to_gene_obj{$type_A_models[0]}, 
                                         $ele_list);
            }
            
            ## Check TYPE_B-merging:
            elsif ($num_type_B > 1 || $num_type_A > 1) {
                my @type_A_objs;
                foreach my $feat (@type_A_models) {
                    push (@type_A_objs, $feat_name_to_gene_obj{$feat});
                }
                
                my @type_B_objs;
                foreach my $feat (@type_B_models) {
                    push (@type_B_objs, $feat_name_to_gene_obj{$feat});
                }
                
                &compare_complex_cases(\@type_A_objs, \@type_B_objs);
            }
        };

        if ($@) {
            print STDERR "Problem comparing genes: " . $@;
        }

            
        foreach my $feat_name (@type_B_models, @type_A_models) {
            delete $feat_names{$feat_name};
        }
        
        
        
        
    }
    
    
    
    
    ## Check for those feat_names that didn't cluster
    foreach my $feat_name (keys %feat_names) {
        #print "Processing feat: $feat_name\n";
        my $gene_obj_ref = $feat_name_to_gene_obj{$feat_name};
        
        if ($feat_name =~ /$type_A/) {
            foreach my $gene_obj  ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
                my $model_feat = $gene_obj->{Model_feat_name};
                print "${type_A}-NOMAP\t$SEQNAME\t$model_feat\n";
            }
        }
        elsif ($feat_name =~ /$type_B/i) {
            foreach my $gene_obj ($gene_obj_ref, $gene_obj_ref->get_additional_isoforms()) {
                my $model_feat = $gene_obj->{Model_feat_name};
                print "${type_B}-NOMAP\t$SEQNAME\t$model_feat\n";
            }
        }
        
        else {
            die "Error, don't recognize feat_name: $feat_name\n";
        }
    }
}



####
sub prepareIDs {
    my ($type, $gene_obj_list_aref) = @_;
    
    foreach my $gene_obj_ele (@$gene_obj_list_aref) {
        foreach my $gene_obj ($gene_obj_ele, $gene_obj_ele->get_additional_isoforms()) {
            my $new_feat_name = $gene_obj->{TU_feat_name} . "::" 
                . $gene_obj->{Model_feat_name} . ":" . $type;
        
            ## Replace model and TU feat_names with new format:
            $gene_obj->{Model_feat_name} = $new_feat_name;
            $gene_obj->{TU_feat_name} = $new_feat_name;
                    
        }
    }
}


####
sub are_compatible {
    my ($gene_obj_A, $gene_obj_B) = @_;
    
    my $assembler = new CDNA::Gene_obj_alignment_assembler(); 
    
    $assembler->assemble_genes($gene_obj_A, $gene_obj_B);
    my @assembly = $assembler->get_assemblies();
    
    #foreach my $assembly (@assembly) {
    #    print $assembly->toToken() . "\n";
    #}
    # 
    # my $assembler = new CDNA::PASA_alignment_assembler;
    # $assembler->{incoming_alignments} = \@assembly;
    # print $assembler->toAlignIllustration();
    # 
    #
    
    if (scalar (@assembly) == 1) {
        return (1);
    }
    else {
        return(0);
    }

}



####
sub compare_complex_cases {
    my ($type_A_models_aref, $type_B_models_aref) = @_;
    
    my $num_type_A_genes = scalar (@$type_A_models_aref);
    my $num_type_B_genes = scalar (@$type_B_models_aref);

    my @feats;
    foreach my $gene_obj_ele (@$type_A_models_aref, @$type_B_models_aref) {
        
        foreach my $gene_obj ($gene_obj_ele, $gene_obj_ele->get_additional_isoforms()) {
            push (@feats, $gene_obj->{Model_feat_name});
        }
    
    }
    

    my $ele_list = join (",", @feats);
    
    
    if ($num_type_A_genes == 1 && $num_type_B_genes > 1) {
        print "SPLIT\t$SEQNAME\t$ele_list\n";
    }

    elsif ($num_type_A_genes > 1 && $num_type_B_genes == 1) {
        print "MERGE\t$SEQNAME\t$ele_list\n";
    }

    elsif ($num_type_A_genes > 1 && $num_type_B_genes > 1) {
        print "COMPLEX\t$SEQNAME\t$ele_list\n";
    }

    else {
        die "Error, not sure what to do here: num type_A genes = $num_type_A_genes, num type_B genes = $num_type_B_genes\n";
    }

}



####
sub compare_isoforms {
    my ($type_A_models_aref, $type_B_models_aref) = @_;
    
    my %classified;


    my %type_A;
    my %type_B;


    ## Do SAME comparisons:

  COMPARE_A_SAME:
    foreach my $type_A_gene_obj (@$type_A_models_aref) {
        
        my $type_A_model = $type_A_gene_obj->{Model_feat_name};
        $type_A{$type_A_model} = 1;
        
        
        if ($classified{$type_A_model}) { next; }
        
      COMPARE_B_SAME:
        foreach my $type_B_gene_obj (@$type_B_models_aref) {
            
            my $type_B_model = $type_B_gene_obj->{Model_feat_name};
            
            $type_B{$type_B_model} = 1;
            
            if ($classified{$type_B_model}) { next; }
            
            &compare_genes($type_A_gene_obj, $type_B_gene_obj);
            if (&are_CDS_same()) {
                print "ISOFORM_SAME\t$SEQNAME\t$type_A_model,$type_B_model\n";
                $classified{$type_A_model} = 1;
                $classified{$type_B_model} = 1;
                last COMPARE_B_SAME;
            }
            
        }
    }
    
    


    ## Do compatible comparisons:
    
  COMPARE_A_COMPATIBLE:
    foreach my $type_A_gene_obj (@$type_A_models_aref) {
        
        my $type_A_model = $type_A_gene_obj->{Model_feat_name};
        $type_A{$type_A_model} = 1;
        
        
        if ($classified{$type_A_model}) { next; }
        
      COMPARE_B_COMPATIBLE:
        foreach my $type_B_gene_obj (@$type_B_models_aref) {
            
            my $type_B_model = $type_B_gene_obj->{Model_feat_name};
            
            $type_B{$type_B_model} = 1;
            
            if ($classified{$type_B_model}) { next; }
            
            if (&are_compatible($type_A_gene_obj, $type_B_gene_obj)) {
                
                $classified{$type_A_model} = 1;
                $classified{$type_B_model} = 1;
                
                
                #print "COMPAT\t$type_A_model,$type_B_model\n";
                classify_compatible($type_A_gene_obj, $type_B_gene_obj);
                
                last COMPARE_B_COMPATIBLE;
                
            }
            
        }
    }
    
    
    ## Examine those unclassified:
    
    my @unclassified_type_A;
    my @unclassified_type_B;
    
    foreach my $gene_obj (@$type_A_models_aref, @$type_B_models_aref) {
        my $feat = $gene_obj->{Model_feat_name};
        unless ($classified{$feat}) {
            if ($type_A{$feat}) {
                push (@unclassified_type_A, $feat);
            } elsif ($type_B{$feat}) {
                push (@unclassified_type_B, $feat);
            }
        }
    }
    
    if (@unclassified_type_A || @unclassified_type_B) {
        my $num_unclassified_type_A = scalar (@unclassified_type_A);
        my $num_unclassified_type_B = scalar (@unclassified_type_B);
        
        ## take the smallest number, make diff
        my $min = ($num_unclassified_type_A < $num_unclassified_type_B) ? $num_unclassified_type_A : $num_unclassified_type_B;
        my @diffs;
        for (1 .. $min) {
            my $type_A_model = shift @unclassified_type_A;
            my $type_B_model = shift @unclassified_type_B;
            print "ISOFORM_DIFF\t$SEQNAME\t" . join (",", $type_A_model, $type_B_model) . "\n";
        }

        
        my $rest = join (",", @unclassified_type_A, @unclassified_type_B);
        
        print "ISOFORM_NOMAP\t$SEQNAME\t$rest\n" if $rest;
    }
}


####
sub classify_compatible {
    my ($type_A_gene_obj, $type_B_gene_obj) = @_;
    
    my ($type_A_lend, $type_A_rend) = sort {$a<=>$b} $type_A_gene_obj->get_model_span();
    my ($type_B_lend, $type_B_rend) = sort {$a<=>$b} $type_B_gene_obj->get_model_span();
    

    my $type = "";
    
    if (
        ($type_B_lend > $type_A_lend && $type_B_rend < $type_A_rend) 
        ||
        ($type_A_lend > $type_B_lend && $type_A_rend < $type_B_rend) 
        ) {
        
        # encapsulation
        $type = "ISOFORM_COMPAT-encaps";
    
    }
    
    elsif ( ($type_B_lend < $type_A_lend && $type_B_rend < $type_A_rend) 
            ||
            ($type_B_lend > $type_A_lend && $type_B_rend > $type_A_rend) 
            )   {
    
        $type = "ISOFORM_COMPAT-staggered";
    }
    
    elsif ( ($type_B_lend == $type_A_lend && $type_B_rend != $type_A_rend) 
            ||
            ($type_B_lend != $type_A_lend && $type_B_rend == $type_A_rend) 
            ) {
        $type = "COMPAT-end_agree";
    }
    
    else {
        $type = "ISOFORM_COMPAT-other";
    }
    
    my $ele_list = $type_A_gene_obj->{Model_feat_name} . "," . $type_B_gene_obj->{Model_feat_name};
    
    print "$type\t$SEQNAME\t$ele_list\n";
}


    
####
sub compare_1_to_1_mapping {
    my ($type_A_model_obj, $type_B_obj, $ele_list) = @_;
    
    
    my @type_A_models;
    foreach my $type_A_model ($type_A_model_obj, $type_A_model_obj->get_additional_isoforms()) {
        push (@type_A_models, $type_A_model);
    }

    my @type_B_models;
    foreach my $type_B_model ($type_B_obj, $type_B_obj->get_additional_isoforms()) {
        push (@type_B_models, $type_B_model);
    }
    
    ## make sure at least two of the isoforms encode a cds in the same frame:
        
    my $got_same_frame_flag = 0;
    
  
  FRAME_ANALYSIS:  
    foreach my $typeA_model (@type_A_models) {
        
        foreach my $typeB_model (@type_B_models) {
            
            my $compare_struct = &Gene_obj_comparator::analyze_cds_frames($typeA_model, $typeB_model);
            
            if ($compare_struct->{coding_coding_same}) {
                $got_same_frame_flag = 1;
            }
            
        }
    }

    
    if ($got_same_frame_flag) {
        &compare_isoforms(\@type_A_models, \@type_B_models);
    }
    else {
        ## No isoform comparison indicated a CDS in the same frame overlapped.
        my @feat_list;
        foreach my $gene_obj (@type_A_models, @type_B_models) {
            my $model_feat = $gene_obj->{Model_feat_name};
            push (@feat_list, $model_feat);
        }
        print "MAP_EXTREME_DIFF\t$SEQNAME\t" . join (",", @feat_list) . "\n";
        
    }
}



####
sub have_exon_overlap {
    my ($model_A, $model_B) = @_;


    foreach my $model_A_exon ($model_A->get_exons()) {
        
        my ($A_lend, $A_rend) = sort {$a<=>$b} $model_A_exon->get_coords();


        foreach my $model_B_exon ($model_B->get_exons()) {
            
            my ($B_lend, $B_rend) = sort {$a<=>$b} $model_B_exon->get_coords();
            
            if ($B_lend < $A_rend && $B_rend > $A_lend) {
                ## good enough. They overlap
                return (1);
            }
        }
    }


    return (0); # no overlap of exons
}


1; #EOM


    
