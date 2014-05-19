#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";

my $in = shift||die;
my $prefix = shift;
$prefix='scaffold\S+' if !$prefix;
mkdir('hints') || die ("hints directory already exists\n");
open (IN,$in)||die;

my $prev_id = '';
while (my $ln=<IN>){
	$ln=~/^($prefix)\b/;
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
