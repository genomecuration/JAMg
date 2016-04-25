#!/usr/bin/env perl

# example perl script to run it automatically.

use strict;
use warnings;
use FindBin qw/$RealBin/;
$ENV{PATH} .= ":$RealBin";

# change this to find all the files from left pairs
my @files = glob("*_1_sequence.fastq");
#my @files = glob("*R1_001.fastq");

foreach my $f (sort @files){
	my $pair = $f;
	# change this to grab the pair's filename by substituting something 
	$pair=~s/_1_sequence/_2_sequence/;
#	$pair=~s/R1_001.fastq/R2_001.fastq/;
	next if !-s $pair || $pair eq $f;
	# change the cmd to something you like
	my $cmd = "preprocess_illumina.pl -paired '$f' '$pair'";
	system($cmd);
	sleep(1);
}

