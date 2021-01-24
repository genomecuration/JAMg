#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
my $file = shift;
my $ids = shift ||die $!;


my %hash;
# "51" "ID=gene1066;Name=5F0FF05B6F1B67DBBC90E77B9BCBD4B7"
open (IN,"$ids")||die;
while (my $ln=<IN>){
	chomp($ln);
	next unless $ln;
	my @data = split(" ",$ln);
	if ($data[-1] && $data[-1]=~/ID=([^;]+)/){
		$hash{$1}++;
	}
}
close IN;
$/ = "\n\n";

open (IN, $file) || die;
while (my $record =<IN>){
	chomp($record);
	next unless $record;
	my @lines = split("\n",$record);
	my $gene = shift(@lines);
	my @data = split("\t",$gene);
	if ($data[8] && $data[8]=~/ID=([^;]+)/){
		my $id = $1;
		$id=~s/\s+$//;
		if (!$hash{$id}){
			print "\n\n".$record;
		}
	}
}
close IN;
