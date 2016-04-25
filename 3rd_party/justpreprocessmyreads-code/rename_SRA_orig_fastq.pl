#!/usr/bin/env perl

use strict;
use warnings;
my @dirs = glob("./?RR*");

foreach my $dir (@dirs){
	next unless -d $dir;
	my $file1 = $dir . "/pass/1/fastq";
	my $file2 = $dir . "/pass/2/fastq";
	rename($file1, "./$dir" . "_1_fastq" ) if (-s $file1);
	rename($file2, "./$dir" . "_2_fastq" ) if (-s $file2);

}
