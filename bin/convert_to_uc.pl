#!/usr/bin/perl -w

use strict;

my $infile=$ARGV[0] || die ();
my $outfile=$infile.".uc";

open (IN,$infile);
open (OUT,">$outfile");

while (my $line=<IN>){
	if ($line=~/^>/ || $line=~/^#/){print OUT $line;}
	else {
		$line=uc($line);
		print OUT $line;
	}
}
close (IN);
close (OUT);
