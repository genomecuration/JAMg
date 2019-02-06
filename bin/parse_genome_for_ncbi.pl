#!/usr/bin/env perl

use strict;
use warnings;

my $file = shift || die;

$/ = ">";

open (IN,$file);
open (OUT,">$file.x");
open (INFO,">$file.x.info");
print INFO "ID\tLength prior 5' trim\tLength after 5' trim\tPadding needed for GFF\n";
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
	$seq=~s/N+$//; # if 5 Ns left at end of seq

	# now the 5' but sometimes we may a gff
	my $orig_size = length($seq);
	$seq=~s/^N+//;
	while ( ( substr($seq,0,350)=~ tr/N/N/) > 150 ){
                $seq = substr($seq,350);
        }
	while ( ( substr($seq,0,50)=~ tr/N/N/) > 15 ){
		$seq = substr($seq,50);
	}
	while ( ( substr($seq,0,15)=~ tr/N/N/) > 5 ){
                $seq = substr($seq,15);
	}
	$seq=~s/^N+//; # 5 Ns left at beginning

	my $new_size = length($seq);
	next if ($new_size < 2000);
        my $Ns = ($seq =~ tr/N//);
	next if ($Ns / $new_size > 0.4) && $new_size < 20000;
	$seq =~s/(\S{80})/$1\n/g;
	print OUT ">$id\n$seq\n";
	if ($new_size != $orig_size){
		print INFO "$id\t$orig_size\t$new_size\t".($orig_size-$new_size)."\n";
	}

}
close IN;
close OUT;
close INFO;
