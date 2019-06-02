#!/usr/bin/env perl

use strict;
use warnings;
my $list = shift;
my $gff = shift||die;
my %hash;
open (IN,$list)||die;
while (my $ln=<IN>){
	chomp($ln);
	my @data = split("\t",$ln);
	$hash{$data[0]} = $data[1] if $data[1];
}
close IN;

open (GFF,$gff)||die;
open (OUT,">$gff.renamed")||die;

while (my $ln=<GFF>){
	my @data = split("\t",$ln);
	if ($hash{$data[0] && $data[8]){
		$data[8] =~s/$data[0]\b/$hash{$data[0]}/g;
		$data[0] = $hash{$data[0]};
		print OUT join("\t",@data);
	}elsif ($data[8]){
		next;
	}else{
		print OUT $ln;
	}
}
close GFF;
