#!/usr/local/bin/perl

use FindBin;
use lib ($FindBin::Bin);
use Pasa_init;
use lib ("../PerlLib");
use strict;
use warnings;
use SingleLinkageClusterer;
require "overlapping_nucs.ph";

my $MIN_PERCENT_OVERLAP = 40;


$/ = "######################################################\n ";

while (<STDIN>) {
    
    if (/^cluster: (\d+)/) {
        my $cluster_no = $1;
        &parse_assemblies($_);
    }
}

exit(0);


####
sub parse_assemblies {
    my ($text) = @_;
    
    my @x = split (/\n/, $text);
    my @assemblies;

    foreach my $line (@x) {
        if ($line =~ /^Assembly/) {
            push (@assemblies, $line);
        }
    }


    # print "ASSEMBLIES:\n" . join ("\n", @assemblies) . "\n\n";
    &cluster_assemblies_via_strand_and_overlap(@assemblies);
    
    return;
}


####
sub cluster_assemblies_via_strand_and_overlap {
    my @assemblies = @_;

    my @assembly_structs;
    
    
    my $counter = 0;
    my %entry_to_struct;
    

    foreach my $assembly (@assemblies) {
        
        $assembly =~ m|orient\(a[\+\-]/s(\S)|;
        my $spliced_orient = $1 or die "Error, cannot extract spliced orientation from $assembly";
        
        my @genome_coords;
        
        while ($assembly =~ m|(\d+)\(\d+\)-(\d+)\(\d+\)|g) {
            my ($lend, $rend) = ($1, $2);
            push (@genome_coords, $lend, $rend);
        }

        @genome_coords = sort {$a<=>$b} @genome_coords;

        my $lend = shift @genome_coords;
        my $rend = pop @genome_coords;

        # print "assembly: $lend, $rend, $spliced_orient\n\n";
        
        $counter++;

        my $struct = { spliced_orient => $spliced_orient,
                       lend => $lend,
                       rend => $rend, 
                       entry => $counter,
                       line => $assembly, # store for later printing.
                   };
        
        push (@assembly_structs, $struct);
        $entry_to_struct{$counter} = $struct;
        
    }
    
    
    ## cluster them.
    if (scalar @assembly_structs == 1) {
        &report_clusters ([@assembly_structs]);
    }
    else {
        my @pairs;
        for (my $i = 0; $i < $#assembly_structs; $i++) {
            
            my $struct_i = $assembly_structs[$i];
            my $i_orient = $struct_i->{spliced_orient};
            my $i_lend = $struct_i->{lend};
            my $i_rend = $struct_i->{rend};
            my $i_length = $i_rend - $i_lend + 1;
            my $i_entry = $struct_i->{entry};

            for (my $j = $i + 1; $j <= $#assembly_structs; $j++) {
                
                my $struct_j = $assembly_structs[$j];
                my $j_orient = $struct_j->{spliced_orient};
                my $j_lend = $struct_j->{lend};
                my $j_rend = $struct_j->{rend};
                my $j_length = $j_rend - $j_lend + 1;
                my $j_entry = $struct_j->{entry};
                

                unless ($i_orient eq $j_orient) { next; }
                
                my $overlapping_nucs = &nucs_in_common($i_lend, $i_rend, $j_lend, $j_rend);

                my $percent_i = $overlapping_nucs / $i_length * 100;
                my $percent_j = $overlapping_nucs / $j_length * 100;

                print "PERCENTAGES: $percent_i, $percent_j\n";
                
                if ($percent_i >= $MIN_PERCENT_OVERLAP || $percent_j >= $MIN_PERCENT_OVERLAP) {

                    push (@pairs, [$i_entry, $j_entry]);
                }
            }
        }
        
        my @assembly_clusters;
        my %seen;
        if (@pairs) {
            my @clusters = &SingleLinkageClusterer::build_clusters(@pairs);
            
            foreach my $cluster (@clusters) {
                my @eles = @$cluster;
                my @struct_cluster;
                foreach my $ele (@eles) {
                    $seen{$ele} = 1;
                    my $struct = $entry_to_struct{$ele} or die "Error, no entry retrieved for $ele\n";
                    push (@struct_cluster, $struct);
                }
                if (@struct_cluster) {
                    push (@assembly_clusters, [@struct_cluster]);
                }
            }
        }
        foreach my $struct (@assembly_structs) {
            my $entry = $struct->{entry};
            unless ($seen{$entry}) {
                push (@assembly_clusters, [$struct]);
            }
        }
    
        &report_clusters(@assembly_clusters);
        
    }
    
    return;
}


####
sub report_clusters {
    my @clusters = @_;

    
    foreach my $cluster (@clusters) {
        my @structs = @$cluster;
        my $cluster_size = scalar @structs;
        print "CLUSTERS SIZE: $cluster_size ";
        my @lines;
        
        foreach my $struct (@structs) {
            my $line = $struct->{line};
            push (@lines, $line);
        }

        print "{ " . join (";", @lines) . "} ";
    }

    print "\n";
}



