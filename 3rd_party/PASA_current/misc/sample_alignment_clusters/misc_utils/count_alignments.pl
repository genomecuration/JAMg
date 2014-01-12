#!/usr/local/bin/perl

use strict;
use warnings;

foreach my $file (@ARGV) {
    my $cmd = "egrep -v \"//\" $file | whitespaceremoval.pl | wc -l";
    my $count = `$cmd`;
    chomp $count;
    print "$file\t$count\n";
}

exit(0);

