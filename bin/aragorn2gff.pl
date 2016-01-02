#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $file = shift||die "Give aragorn output file";
my $out = $file.'.gff3';
open (IN,$file);
open (OUT,">$out");
my $orig_sep = $/;
$/=">";
my $discard = <IN>;
my %pseudo_counter;
my %trna_counter;
print OUT "gff-version 3\n";
while (my $record=<IN>){
	chomp($record);
	my @lines = split("\n",$record);
	my $id = $lines[0];
	my $print;
	for (my $i=2;$i<scalar(@lines);$i++){
		my @data = split(/\s+/,$lines[$i]);
		my ($counter,$pseudo);
		my $hit = $data[1];
		my $type = 'tRNA';
		my $hit_id;
		if ($hit=~/\*$/){
			$type = 'pseudogene_related_to_tRNA';
			$pseudo = 1;
			chop($hit);
			$hit .= ":pseudogene";
			$pseudo_counter{$hit}++;
			$counter=$pseudo_counter{$hit};
			$hit_id = $hit;
		}else{
			$trna_counter{$hit}++;
			$counter=$trna_counter{$hit};
			$hit_id = $hit;
		}
		my $strand = $data[2]=~/^c/ ? '-' : '+';
		my $score = int($data[3]);
		$data[2]=~/(\d+),(\d+)/;
		my $start =$1;
		my $end = $2;
   		$data[4]=~/^\(([a-z]{3})\)/;
		my $anticodon = $1;
		$data[4] =~/i\((\d+),(\d+)\)/;
		my $intron_start = $1 if $1;
		my $intron_length = $2 if $2;
		$print .= $id."\tARAGORN\tsimilarity\t$start\t$end\t$score\t$strand\t.\tID=$hit_id.$counter;Alias=$hit";

		$print .= ";anticodon=$anticodon" if $anticodon;
		$print .= ";is_pseudo=yes" if $pseudo;
		$print .= "\n";
	}
	print OUT $print if $print;

}
close IN;
close OUT;

$/ = $orig_sep;
print "\nCompleted. Found these tRNAs:\n";
print Dumper \%trna_counter;
print "\nAnd these pseudogenes related to tRNAs:\n";
print Dumper \%pseudo_counter;
print "\n";
