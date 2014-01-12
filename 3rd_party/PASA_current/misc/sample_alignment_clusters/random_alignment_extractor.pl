#!/usr/local/bin/perl

use strict;
use warnings;

srand();

my $usage = "usage: $0 alignment_clusters_file target_number\n\n";

my $alignment_cluster_file = $ARGV[0] or die $usage;
my $target_number_of_entries = $ARGV[1] or die $usage;

my $total_num_entries = &count_alignments($alignment_cluster_file);
my $target_ratio = $target_number_of_entries / $total_num_entries;

print STDERR "total entries: $total_num_entries, target_ratio = $target_ratio\n";

my $cluster_id = "";
my @alignments;


main: {
 
    open (my $fh, "$alignment_cluster_file") or die $!;
    while (<$fh>) {
        unless (/\w/) { next; }
        
        if (/^\//) {
            &process_alignments() if @alignments;
            $cluster_id = $_;
            @alignments = (); # clear it.
        }
        else {
            push (@alignments, $_);
        }
    }
    
    &process_alignments() if @alignments;

    exit(0);
}

####
sub process_alignments {
    
    my $alignment_text = "";
    
    foreach my $alignment (@alignments) {
        my $prob = rand(1);
        if ($prob < $target_ratio) { 
            $alignment_text .= $alignment;
        }
    }

    if ($alignment_text) {
        print "$cluster_id" . $alignment_text . "\n\n";
    }

    return;
}


####
sub count_alignments {
    my $file = shift;
    
    my $count = 0;

    open (my $fh, $file) or die $!;
    while (<$fh>) {
        if (/^\//) { next; }
        unless (/\w/) { next; }
        $count++;
    }
    close $fh;

    return ($count);
}
        
        
