#!/usr/bin/env perl

use strict;
use warnings;
my @dirs = glob("./?RR*");

foreach my $dir (@dirs){
	next unless -d $dir;
	my @files = glob($dir . "/pass/?/fastq");
	foreach my $file (@files){
		if ($file=~/\/pass\/(\d)\//){
			my $id = $1 || die;
			my $new_name = "./$dir" . "_$id"."_fastq";
			warn "SKIP: File $new_name already exists\n" if -s $new_name;
			rename($file, $new_name ) if !-s $new_name;
		}
	}
}
