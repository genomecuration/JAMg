#!/usr/bin/env perl

use strict;
use warnings;

my $fh;
my $prev_ev_type;

while (<STDIN>) {
    unless (/\w/) { next; }
    my @x = split (/\t/);
    my $ev_type = $x[1];
    if ($prev_ev_type ne $ev_type) {
        $prev_ev_type = $ev_type;
        if ($fh) { close $fh; }
        print "writing $ev_type\n";
        open ($fh, ">>$ev_type.gff3") or die $!;
    }
    print $fh $_;
}

exit(0);

