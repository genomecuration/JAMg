#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $file = shift||die "Give tRNAScan output file";
my $cutoff = shift;
my $out = $file.'.gff3';
$cutoff = 22 if !$cutoff;

open (IN,$file);
open (OUT,">$out");
print OUT "gff-version 3\n";
my $discard = <IN>.<IN>.<IN>;
my %hash;
while (my $ln=<IN>){
	chomp($ln);
	my @data = split("\t",$ln);
	next unless $data[8];
	my ($start,$end,$strand) = ($data[3] > $data[2]) ? ($data[2],$data[3],'+') : ($data[3],$data[2],'-');
	my $trna = $data[4];
	next if $trna eq 'Undet';
	my $anticodon = $data[5];
 	$hash{$trna}++;
	my $id = $trna.':tRNAscan:'. $hash{$trna};
	my $score = $data[8];
	next if $score < $cutoff;
	my $ref = $data[0];
	# ideally we would like to also capture introns with match/match part, but another day..
	print OUT "$ref\ttRNAscan\tsimilarity\t$start\t$end\t$score\t$strand\t.\tID=$id;Name=tRNA-$trna;anticodon=$anticodon\n";
}
close IN;
close OUT;

