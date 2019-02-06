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
	my $seq = join("\n",@lines);
	if ($seq=~/^N+/){
		die "Sequence $id has Ns in the start of the sequence\n";
	}
	$seq=~s/N+$//;
	next if (length($seq) < 2000);
        my $Ns = ($seq =~ tr/N//);
	next if ($Ns / length($seq) > 0.4) && length($seq) < 20000;
	print OUT ">$id\n$seq\n";

}
close IN;
close OUT;
