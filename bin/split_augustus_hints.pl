#!/usr/bin/env perl

use strict;
use warnings;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";

my $in = shift||die ("Please provide a hint file\n");

mkdir('hints') || die ("hints directory already exists\n");
open (IN,$in)||die;
my $prev_id = '';
while (my $ln=<IN>){
	next if $ln=~/^#/ || $ln=~/^\s*$/;
	my @data = split("\t",$ln);
	next unless $data[8];
	my $id = $data[0];
	my $outfile = "hints/$id.hints";
	if ($id ne $prev_id){
		close OUT;
		open (OUT,'>>'.$outfile);
	}
	print OUT $ln;
       	$prev_id = $id;
}
close IN;
close OUT;


my @outputs = glob("hints/*");
print "Produced hints for ".$#outputs. " reference sequences\n";
