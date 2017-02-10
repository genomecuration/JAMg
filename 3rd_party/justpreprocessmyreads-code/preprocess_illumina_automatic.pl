#!/usr/bin/env perl

# example perl script to run it automatically.

use strict;
use warnings;
use FindBin qw/$RealBin/;
$ENV{PATH} .= ":$RealBin";

my ($is_sra,$is_seqcenter1,$is_seqcenter2) = (1,0,0);

my $option = shift;
my $extras = join(' ',@ARGV);


$is_sra = 1 if $option && $option=~/sra/i;
$is_seqcenter1 = 1 if $option && $option=~/seqcenter1/i;
$is_seqcenter2 = 1 if $option && ($option=~/seqcenter2/i || $option=~/rama/i  ) ;

# change this to find all the files from left pairs
my @files;
@files = glob("*_1_fastq*") if $is_sra; #SRA
@files = glob("*_1_sequence.fastq*") if $is_seqcenter1;
@files = glob("*R1_*.fastq*") if $is_seqcenter2;


foreach my $f (sort @files){
	my $cmd = "$RealBin/preprocess_illumina.pl ";
	my $pair = $f;
	# change this to grab the pair's filename by substituting something 
	$pair=~s/_1_fastq/_2_fastq/ if $is_sra; #SRA
	$pair=~s/_1_sequence/_2_sequence/ if $is_seqcenter1;
	$pair=~s/R1_/R2_/ if $is_seqcenter2;
	if (!-s $pair || $pair eq $f){
		$cmd .= " '$f'";
	}else{
		$cmd .= " -paired $extras '$f' '$pair'";
	}
	# can change the cmd to something you like
	print "CMD: $cmd\n";
	sleep(3);
	system($cmd);
	sleep(1);
}

