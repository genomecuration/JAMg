#!/usr/local/bin/perl

use strict;
use warnings;

my $usage = "\n\nusage: $0 cluster_id\n\n";

my $cluster_id = $ARGV[0] or die $usage;

while (<STDIN>) {
    if (m|// cluster: (\d+)|) {
        if ($1 == $cluster_id) {

            my $line = <STDIN>;
            while ($line !~ m|^//|) {
                print $line;
                $line = <STDIN>;
            }
            last;
        }
    }
}


