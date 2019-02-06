#!/usr/bin/env perl

use strict;
use warnings;

my $file = shift || die;

$/ = ">";

open (IN,$file);
open (OUT,">$file.x");

while (my $record=<IN>){
	chomp($record);
	next unless $record;
	my @lines = split("\n",$record);
	my $id = shift(@lines);
	my $seq = uc(join('',@lines));
	if ($seq=~/^N+/){
		die "Sequence $id has Ns in the start of the sequence\n";
	}
	$seq=~s/N+$//;
#	$seq=~s/N{5,}[ATCG]{0,10}$//g;
#	$seq=~s/N{15,}[ATCG]{0,35}$//g;
#	$seq=~s/N{150,}[ATCG]{0,200}$//g; # ~40%
#	#assume not and much quicker
	while ( ( substr($seq,-350)=~ tr/N/N/) > 150 ){
                $seq = substr($seq,0,length($seq)-351);
        }
	while ( ( substr($seq,-50)=~ tr/N/N/) > 15 ){
		$seq = substr($seq,0,length($seq)-51);
	}
	while ( ( substr($seq,-15)=~ tr/N/N/) > 5 ){
                $seq = substr($seq,0,length($seq)-16);
	}

	next if (length($seq) < 2000);
        my $Ns = ($seq =~ tr/N//);
	next if ($Ns / length($seq) > 0.4) && length($seq) < 20000;
	$seq =~s/(\S{80})/$1\n/g;
	print OUT ">$id\n$seq\n";

}
close IN;
close OUT;
