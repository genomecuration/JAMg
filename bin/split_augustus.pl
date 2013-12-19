#!/usr/bin/env perl

use strict;
use warnings;
my $in = shift||die;
mkdir('hints') || die ("hints directory already exists\n");
open (IN,$in)||die;

my $prev_id = '';
while (my $ln=<IN>){
	$ln=~/^(scaffold\S+)\b/;
	next unless $1;
	my $outfile = "hints/$1.hints";
	if ($1 ne $prev_id){
		close OUT;
		open (OUT,'>>'.$outfile);
	}
	print OUT $ln;
        $prev_id = $1;
}
close IN;
close OUT;
