#!/usr/bin/env perl

use strict;
use warnings;
my $file=shift || die;

open (IN,$file);
while (my $ln=<IN>){
	print "\n" if $ln=~/\tgene\t/;
	print $ln;
}

close IN;
