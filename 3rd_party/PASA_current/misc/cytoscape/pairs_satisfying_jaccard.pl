#!/usr/local/bin/perl

use strict;
use warnings;

my $usage = "usage: $0 pairs clusters\n\n";

my $pairs_file = $ARGV[0] or die $usage;
my $clusters = $ARGV[1] or die $usage;

my %OK;

{ # read pairs
    my @clusters = `cat $clusters`;
    chomp @clusters;
    foreach my $cluster (@clusters) {
        my @entries = split (/\s+/, $cluster);
        
        foreach my $eA (@entries) {
            
            foreach my $eB (@entries) {
                unless ($eA eq $eB) {
                    $OK{$eA}->{$eB} = 1;
                    $OK{$eB}->{$eA} = 1;
                }
            }
        }
    }
}

{ # process pairs, only those in jaccard clusters are passed thru
    open (my $fh, $pairs_file) or die "Error, cannot open $pairs_file";
    while (<$fh>) {
        my $line = $_;
        chomp;
        my ($eA, $eB) = split (/\s+/);
        if ($OK{$eA}->{$eB}) {
            print $line;
        }
    }
}


exit(0);


             
        
