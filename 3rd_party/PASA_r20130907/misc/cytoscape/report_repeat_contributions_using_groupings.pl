#!/usr/local/bin/perl

use strict;
use warnings;


my %entry_to_jaccard_cluster;
my @jaccard_clusters;
my %entry_to_repeat_content;

my $usage = "usage: $0 repeat_sums jaccard_clusters\n\n";

my $repeat_sums_file = $ARGV[0] or die $usage;
my $jaccard_clusters_file = $ARGV[1] or die $usage;

{  # read jaccard clusters
    open (my $fh, $jaccard_clusters_file) or die "Error, cannot open $jaccard_clusters_file";
    while (<$fh>) {
        chomp;
        
        my @entries = split (/\s+/);
        my $jaccard_count = scalar @jaccard_clusters;

        push (@jaccard_clusters, \@entries);
        foreach my $entry (@jaccard_clusters) {
            $entry_to_jaccard_cluster{$entry} = $jaccard_count;
        }

    }

}

open (my $fh, $repeat_sums_file) or die "Error, cannot open $repeat_sums_file";
while (<$fh>) {
    my $line = $_;
    chomp;
    my ($repeat_acc, $sum_len) = split (/\t/);
    $entry_to_repeat_content{$repeat_acc} = $sum_len;

    unless ($entry_to_jaccard_cluster{$repeat_acc}) {
        print $line;
    }

}

# report the jaccard clusters
my $count = 0;
foreach my $jaccard_cluster (@jaccard_clusters) {
    
    $count++;

    my @entries = @$jaccard_cluster;

    my $sum_length = 0;

    my $annot_text = "";

    foreach my $entry (@entries) {
        my $length = $entry_to_repeat_content{$entry}||0;
        $sum_length += $length;
        $annot_text .= "$entry\[$length], ";
    }
        
    print "Jaccard_$count\t$sum_length\t$annot_text\n";
}


exit(0);


