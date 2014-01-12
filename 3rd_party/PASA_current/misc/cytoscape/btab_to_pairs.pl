#!/usr/local/bin/perl

use strict;
use warnings;

while (<STDIN>) {
    my @x = split (/\t/);
    print "$x[0]\t$x[5]\n" if ($x[0] && $x[5]);
}

exit(0);

