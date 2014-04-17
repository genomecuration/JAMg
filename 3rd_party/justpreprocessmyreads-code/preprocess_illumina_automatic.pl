#!/usr/bin/env perl

# example perl script to run it automatically.

use strict;
use warnings;

# change this to find all the files from left pairs
my @files = glob("*_1_sequence.fastq");

foreach my $f (sort @files){
	my $pair = $f;
	# change this to grab the pair's filename by substituting something 
	$pair=~s/_1_sequence/_2_sequence/;
	# change the cmd to something you like
	my $cmd = "preprocess_illumina.pl -cdna -paired $f $pair";
	system($cmd);
	sleep(1);
}

