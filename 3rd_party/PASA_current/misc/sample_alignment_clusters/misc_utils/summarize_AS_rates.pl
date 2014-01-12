#!/usr/local/bin/perl

use strict;
use warnings;

main: {
    foreach my $file (@ARGV) {
        
        my $num_clusters = 0;
        my $num_AS_clusters = 0;

        $file =~ /(\d+)k\./;
        my $num_sequences = $1;
        $num_sequences *= 1000;


        open (my $fh, $file) or die $!;
        while (<$fh>) {
            if (  /CLUSTERS SIZE: (\d+)/) {
                my $cluster_size = $1 or die "Error, cannot extract cluster size from $_" ;
                $num_clusters++;
                if ($cluster_size > 1) {
                    $num_AS_clusters ++;
                }
            }
        }
        close $fh;

        
        print $file . "\tNum_seqs: $num_sequences\tRate AS:\t" . sprintf ("%.2f", $num_AS_clusters / $num_clusters * 100) . "\n";
        
    }

    
    exit(0);
    

}

