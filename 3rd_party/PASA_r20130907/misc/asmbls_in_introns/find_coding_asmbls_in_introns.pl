#!/usr/local/bin/perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES});
use Data::Dumper;
use Overlap_piler;


my $usage = "usage: $0 pasa_training.gff_file\n\n";
my $gffFile = $ARGV[0] or die $usage;


my $MIN_INTRON_LENGTH = 5000;

my %contig_to_asmbls; #alignment structs
my %contig_to_accs; #just accs
my %asmbl_to_coords;
## Parse the alignments for pasa assemblies:
open (my $fh, $gffFile) or die $!;
while (<$fh>) {
    chomp;
    my @x = split (/\t/);
    if ($x[2] !~ /CDS/) { 
        next;
    }
    my ($genome_acc, $end5, $end3, $orient, $pasa_asmbl) = ($x[0], $x[3], $x[4], $x[6], $x[8]);
    
    if ($orient eq '-') {
        ($end5, $end3) = ($end3, $end5);
    }

    my $coords_aref = $asmbl_to_coords{$pasa_asmbl};
    unless (ref $coords_aref) {
        $coords_aref = $asmbl_to_coords{$pasa_asmbl} = [];
    }
    push (@$coords_aref, [$end5,$end3]);

    $contig_to_accs{$genome_acc}->{$pasa_asmbl} = 1;
}
close $fh;


foreach my $contig (keys %contig_to_accs) {
    
    my @pasa_accs = keys %{$contig_to_accs{$contig}};
    foreach my $pasa_acc (@pasa_accs) {
        my $exons_aref = $asmbl_to_coords{$pasa_acc};
        
        my $alignment_struct = &build_alignment_struct ($pasa_acc, $contig, $exons_aref);
        
        unless ($alignment_struct) {
            next;
        }

        #print Dumper ($alignment_struct);
        
        my $contig_to_asmbls_listref = $contig_to_asmbls{$contig};
        unless (ref $contig_to_asmbls_listref) {
            $contig_to_asmbls_listref = $contig_to_asmbls{$contig} = [];
        }
        push (@$contig_to_asmbls_listref, $alignment_struct); 
    }
}


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
            my $num_accs = scalar (@clustered_accs);
            #print "CONTIG: $contig\tNUM_ACCS: $num_accs\t@clustered_accs\n";
            
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


            #print "$acc has long intron.\n";
            
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
    my ($pasa_asmbl, $genome_acc, $exons_aref) = @_;
    my $exons = Dumper ($exons_aref);
    $exons =~ s/\n//g;
    
    #print "Analyzing: $pasa_asmbl\t$exons\n"; 
    
    ## only analyzing multi-exon assemblies
    my @coordsets;
    my $spliced_orient = "";
    foreach my $coords_aref (@$exons_aref) {
        my ($end5, $end3) = @$coords_aref;
        if ($end5 < $end3) {
            $spliced_orient = '+';
        } 

        elsif ($end5 > $end3) {
            $spliced_orient = '-';
            ($end5, $end3) = ($end3, $end5);
        }
        push (@coordsets, [$end5, $end3]);
    }
    
    @coordsets = sort {$a->[0]<=>$b->[0]} @coordsets;
    
    #print "$pasa_asmbl ($spliced_orient)\n";
    
    if ($spliced_orient eq "") {
        return (undef); 
    }
        
    if (! @coordsets) {
        die "Error, no coordsets extracted from $pasa_asmbl ";
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
        
        #print "intron: $pasa_asmbl\t$intron_length\n";

        if ($intron_length >= $MIN_INTRON_LENGTH) {
            $has_long_intron_flag = 1;
            push (@long_introns, [$intron_lend, $intron_rend] );
        }
        $rend = $next_coordset_ref->[1]; #reset for next comparison
    }
    
    if (! $has_intron_flag) {
        die "Error, no introns for $pasa_asmbl ";
    }
    
    #print "INTRONS: $pasa_asmbl\t$has_long_intron_flag\n";
    
    
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





