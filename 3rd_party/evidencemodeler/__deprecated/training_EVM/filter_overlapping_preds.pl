#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use Gene_obj;
use Gene_obj_indexer;

our $SEE = 1;

my $usage = "usage: $0 file.inx\n\n";

my $gene_inx = $ARGV[0] or die $usage;

my $gene_obj_indexer = new Gene_obj_indexer( { "use" => $gene_inx } );

my %asmbl_id_to_gene_list;

foreach my $gene_id ($gene_obj_indexer->get_keys()) {
    my $gene_obj = $gene_obj_indexer->get_gene($gene_id);

    my $asmbl_id = $gene_obj->{asmbl_id};
    my $orient = $gene_obj->get_orientation();
    $asmbl_id .= $orient;

    my ($lend, $rend) = sort {$a<=>$b} $gene_obj->get_model_span();

    # print "parsing: $asmbl_id\t$lend-$rend\n";

    my $gene_struct = [$gene_id, $lend, $rend];

    my $gene_list_aref = $asmbl_id_to_gene_list{$asmbl_id};
    unless (ref $gene_list_aref) {
        $gene_list_aref = $asmbl_id_to_gene_list{$asmbl_id} = [];
    }

    push (@$gene_list_aref, $gene_struct);

}

# print "\n\nFiltering for overlaps.\n";


foreach my $asmbl_id (keys %asmbl_id_to_gene_list) {
    my @gene_structs = sort {$a->[1]<=>$b->[1]} @{$asmbl_id_to_gene_list{$asmbl_id}};
    
    my @nonoverlapping_structs;

    for (my $i = 0; $i <= $#gene_structs; $i++) {
        
        my $struct_A = $gene_structs[$i];
        my $struct_B = $gene_structs[$i+1];
        
        if ($i == $#gene_structs || !&overlaps($struct_A, $struct_B)) {
            push (@nonoverlapping_structs, $struct_A);
        }
        else {
            # overlapping structs.
            # choose either A or B
            if (&encapsulates($struct_A, $struct_B)) {
                push (@nonoverlapping_structs, $struct_A);
            }
            elsif (&encapsulates($struct_B, $struct_A) ) {
                push (@nonoverlapping_structs, $struct_B);
            }
            else {
                # they better overlap!
                unless (&overlaps($struct_A, $struct_B) ) {
                    die "Error, supposed to overlap but don't! ";
                }

                if (&get_length($struct_A) >= &get_length($struct_B)) {
                    push (@nonoverlapping_structs, $struct_A);
                }
                else {
                    push (@nonoverlapping_structs, $struct_B);
                }
            }
            $i++; # skip past B
        }
    }

    foreach my $nonoverlapping_struct (@nonoverlapping_structs) {
        # print join ("\t", @$nonoverlapping_struct) . "\n";
        my ($acc, $lend, $rend) = @$nonoverlapping_struct;
        my $gene_obj = $gene_obj_indexer->get_gene($acc);
        print $gene_obj->to_GFF3_format() . "\n";
        
    }
    
}



exit(0);

####
sub overlaps {
    my ($struct_A, $struct_B) = @_;

    my ($acc_A, $lend_A, $rend_A) = @$struct_A;
    my ($acc_B, $lend_B, $rend_B) = @$struct_B;

    if ($lend_A < $rend_B && $rend_A > $lend_B) {
        return (1);
    }
    else {
        return (0);
    }
}

###
sub get_length {
    my ($struct) = @_;
    my ($acc, $lend, $rend) = @$struct;
    return ($rend - $lend + 1);
}

###
sub encapsulates {
    my ($struct_A, $struct_B) = @_;

    
    my ($acc_A, $lend_A, $rend_A) = @$struct_A;
    my ($acc_B, $lend_B, $rend_B) = @$struct_B;

    if ($lend_A <= $rend_A && $lend_B >= $rend_B) {
        return (1);
    }
    else {
        return (0);
    }
}


