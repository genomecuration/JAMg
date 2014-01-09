#!/usr/bin/env perl

use strict;
use warnings;
my $gff = shift;
my $predgff = shift||die;

my %hash;
open (IN,$predgff)||die;
while (my $ln=<IN>){
	next if $ln=~/^#/;
	$ln=~/^(\S+)/;
	my $id=$1;
	$hash{$id} = 1;
}
close IN;
open (IN,$gff);
open (OUT,">$predgff.valid");
while (my $ln=<IN>){
	if ($ln=~/^#/){
		print OUT $ln;next;
	}
	$ln=~/^(\S+)/;
        my $id=$1;
	print OUT $ln if $hash{$id};
}
close IN;
close OUT;

