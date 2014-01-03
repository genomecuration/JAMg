#!/usr/local/bin/perl

use strict;
use warnings;

while (<STDIN>) {
    chomp;
    my ($eA, $eB) = split (/\s+/);
    print "$eA\tpp\t$eB\n";
}

exit(0);

