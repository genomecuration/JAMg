#!/usr/bin/env perl

use strict;
use warnings;

my $list = shift;
my $gff = shift;

die "Provide a LIST and a gff to split into stdout and stderr\n" unless $list && -s $list && $gff && -s $gff;
my %hash;
open (IN,$list);
while (my $ln=<IN>){
	chomp($ln);
	$ln=~s/^>//;
	$ln=~s/^\s+//;
	$ln=~s/\s+.+//;
	next unless $ln;
	$hash{$ln}++;
}
close (IN);

my $orig_sep = $/;
$/ = "\n\n";

open (GFF,$gff);
while (my $record=<GFF>){
	chomp($record);
	next unless $record;
	my @lines = split($orig_sep,$record);
	my @check = split("\t",$lines[0]);
	if ($hash{$check[0]}){
		print $record.$/;
	}else{
		warn $record.$/;
	}
}
close (GFF);

