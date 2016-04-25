#!/usr/bin/env perl

use strict;
use warnings;

my @dirs = glob("?RR*");

foreach my $dir (@dirs){
	next unless -d $dir;
	my $file1 = $dir . "/pass/1/fastq.trimmomatic.bz2";
	my $file2 = $dir . "/pass/2/fastq.trimmomatic.bz2";
	my $file1u = $dir . "/pass/1/fastq.unpaired.bz2";
	my $file2u = $dir . "/pass/2/fastq.unpaired.bz2";
	rename($file1, "./$dir" . "_1_fastq.trimmomatic.bz2" ) if (-s $file1);
	rename($file2, "./$dir" . "_2_fastq.trimmomatic.bz2" ) if (-s $file2);
	rename($file1u, "./$dir" . "_left_unpaired_fastq.trimmomatic.bz2" ) if (-s $file1u);
	rename($file2u, "./$dir" . "_right_unpaired_fastq.trimmomatic.bz2" ) if (-s $file2u);

	$file1 = $dir . "/pass/1/fastq.trimmomatic";
	$file2 = $dir . "/pass/2/fastq.trimmomatic";
	$file1u = $dir . "/pass/1/fastq.unpaired";
	$file2u = $dir . "/pass/2/fastq.unpaired";
	rename($file1, "./$dir" . "_1_fastq.trimmomatic" ) if (-s $file1);
	rename($file2, "./$dir" . "_2_fastq.trimmomatic" ) if (-s $file2);
	rename($file1u, "./$dir" . "_left_unpaired_fastq.trimmomatic" ) if (-s $file1u);
	rename($file2u, "./$dir" . "_right_unpaired_fastq.trimmomatic" ) if (-s $file2u);
}
