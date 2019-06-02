#!/usr/bin/env perl

use strict;
use warnings;
my $list = shift;
my $bg = shift||die;
my %hash;
open (IN,$list)||die;
while (my $ln=<IN>){
	chomp($ln);
	my @data = split("\t",$ln);
	$hash{$data[0]} = $data[1] if $data[1];
}
close IN;

open (BG,$bg)||die;
open (OUT,">$bg.renamed")||die;

while (my $ln=<BG>){
	my @data = split("\t",$ln);
	if ($hash{$data[0]} && $data[3]){
		$data[0] = $hash{$data[0]};
		print OUT join("\t",@data);
	}
}
close BG;
close OUT;
