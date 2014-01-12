#!/usr/local/bin/perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES});
use Data::Dumper;
use Overlap_piler;


my $usage = "usage: $0 assemblies_described.txt\n\ndownloaded text from pasa website\n\n";
my $asmbls_text = $ARGV[0] or die $usage;


my $MIN_INTRON_LENGTH = 5000;

my %contig_to_asmbls;

## Parse the alignments for pasa assemblies:
open (my $fh, $asmbls_text) or die $!;
while (<$fh>) {
    chomp;
    my $line = $_;
    my ($pasa_asmbl, $cdna_listing, $genome_acc, $alignment_string) = split (/\t/);
    
    
    my $alignment_struct = &build_alignment_struct ($pasa_asmbl, $genome_acc, $alignment_string);
    
    unless ($alignment_struct) {
        next;
    }
    
    my $contig_to_asmbls_listref = $contig_to_asmbls{$genome_acc};
    unless (ref $contig_to_asmbls_listref) {
        $contig_to_asmbls_listref = $contig_to_asmbls{$genome_acc} = [];
    }
    push (@$contig_to_asmbls_listref, $alignment_struct); 
}
close $fh;


## Determine which ones are in introns of other assemblies:
foreach my $contig (keys %contig_to_asmbls) {
    
    my $assembly_listref = $contig_to_asmbls{$contig};
    
    ## cluster the overlapping assemblies
    my $overlapper = new Overlap_piler();
    my %acc_to_alignment_struct;
    foreach my $alignment_struct (@$assembly_listref) {
        
        my ($acc, $lend, $rend) = ($alignment_struct->{acc},
                                   $alignment_struct->{lend},
                                   $alignment_struct->{rend});
        
        $overlapper->add_coordSet($acc, $lend, $rend);
        
        $acc_to_alignment_struct{$acc} = $alignment_struct;
    }

    my @clusters = $overlapper->build_clusters();
    foreach my $cluster (@clusters) {
        my @clustered_accs = @$cluster;
        if (scalar @clustered_accs > 1) {
            
            ## examine the alignments for intron encapsulation:
            my @alignments_to_examine;
            foreach my $acc (@clustered_accs) {
                my $alignment_struct = $acc_to_alignment_struct{$acc};
                push (@alignments_to_examine, $alignment_struct);
            }
            
            &analyze_genes_in_introns(@alignments_to_examine);
        }
    }
}

exit(0);

####
sub analyze_genes_in_introns {
    my @alignment_structs = @_;
    
    my $num_overlaps = scalar (@alignment_structs);
    #print "NUM_OVERLAP = $num_overlaps\n";
    
    my $cluster_reference = "CLUSTER_REFERENCE";
    
    foreach my $alignment (@alignment_structs) {
        if ($alignment->{has_long_intron_flag}) {
            my $acc = $alignment->{acc};
            my $genome_acc = $alignment->{genome_acc};
            my $spliced_orient = $alignment->{spliced_orient};

            my @long_introns = @{$alignment->{long_introns}};
            foreach my $long_intron (@long_introns) {
                my ($intron_lend, $intron_rend) = @$long_intron;
                
                foreach my $other_alignment (@alignment_structs) {
                    my $other_acc = $other_alignment->{acc};
                    my $other_spliced_orient = $other_alignment->{spliced_orient};

                    if ($acc eq $other_acc) {
                        next; # no self comparisons
                    }
                    
                    my $orient_compare = ($spliced_orient eq $other_spliced_orient) ? "SAME_ORIENT" : "DIFF_ORIENT";

                    my ($lend, $rend) = ($other_alignment->{lend}, $other_alignment->{rend});
                    if ($lend > $intron_lend && $rend < $intron_rend) {
                        # found in intron:
                        print "$genome_acc\t$other_acc ($lend-$rend,$spliced_orient) within intron of $acc ($intron_lend-$intron_rend,$other_spliced_orient) $orient_compare $cluster_reference\n";
                        $cluster_reference = "";
                    }
                }
            }
        }
    }
}



####
sub build_alignment_struct {
    my ($pasa_asmbl, $genome_acc, $alignment_string) = @_;
    
    ## only analyzing multi-exon assemblies


    ## get spliced orientation
    $alignment_string =~ /orient\(a.\/s(.)/;
    my $spliced_orient = $1 or die "Error, no spliced orient for $alignment_string ";
    
    if ($spliced_orient eq '?') {
        return (undef);
    }
    

    ## get alignment segments
    my @coordsets;
    while ($alignment_string =~ /(\d+)\(\d+\)-(\d+)\(/g) {
        my ($lend, $rend) = ($1, $2);
        push (@coordsets, [$lend,$rend]);
    }
    
    if (! @coordsets) {
        die "Error, no coordsets extracted from $alignment_string ";
    }

    if (scalar @coordsets == 1) {
        return (undef); # want more than one exon
    }

    my $min_coord = $coordsets[0]->[0];
    my $max_coord = $coordsets[$#coordsets]->[1];
    

    ## get intron coordinates
    my @long_introns;

    my $first_coordset_ref = $coordsets[0];
    my $rend = $first_coordset_ref->[1];
    
    my $has_long_intron_flag = 0;
    my $has_intron_flag = 0;

    for (my $i=1; $i <= $#coordsets; $i++) {
        
        $has_intron_flag = 1;
        
        my $next_coordset_ref = $coordsets[$i];
        my $next_lend = $next_coordset_ref->[0];
        
        my $intron_lend = $rend+1;
        my $intron_rend = $next_lend-1;
        
        my $intron_length = $intron_rend - $intron_lend + 1;
        if ($intron_length >= $MIN_INTRON_LENGTH) {
            $has_long_intron_flag = 1;
            push (@long_introns, [$intron_lend, $intron_rend] );
        }
        $rend = $next_coordset_ref->[1]; #reset for next comparison
    }
    
    if (! $has_intron_flag) {
        die "Error, no introns for $alignment_string ";
    }
    
    
    my $alignment_struct = { acc => $pasa_asmbl,
                             genome_acc => $genome_acc,
                             lend => $min_coord,
                             rend => $max_coord,
                             spliced_orient => $spliced_orient,
                             has_long_intron_flag => $has_long_intron_flag,
                             long_introns => \@long_introns,
                         };
    
    return ($alignment_struct);
    

}





