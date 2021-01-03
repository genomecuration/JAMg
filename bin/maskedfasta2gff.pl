#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $fasta = shift || die $!;
my %hash;

my $orig_sep = $/;

$/ = '>';
open (IN,$fasta) || die $!;
open (OUTH,">$fasta.hardmasked");
open (OUTS,">$fasta.softmasked");
open (OUTN,">$fasta.normal");
open (BED,">$fasta.repeats.bed");
open (GFF3,">$fasta.repeats.gff3");
while (my $record = <IN>){
	chomp($record); next unless $record;
	my @lines = split("\n",$record);
	next unless $lines[1];
	my $id = shift(@lines);
	my $softmasked_seq = join('',@lines);
	$softmasked_seq=~s/\s+$//;
	$softmasked_seq=~s/^\s+//;
	$softmasked_seq=~tr/n/N/;
	my $hardmasked_seq = $softmasked_seq =~ tr/atgc/NNNN/r;
	my $normal_seq = $softmasked_seq =~ tr/atgc/ATGC/r;	
	print OUTH ">$id\n$hardmasked_seq\n";
	print OUTN ">$id\n$normal_seq\n";
	print OUTS ">$id\n$softmasked_seq\n";

	my $counter;
	while ($softmasked_seq =~ /([a-z]+)/g){
		my $repeat = $1; 
		my $region_end = pos($softmasked_seq);
		$counter++;
		my $region_start = $region_end - length($repeat);
		print BED "$id\t$region_start\t$region_end\n";
		print GFF3 "$id\tmaskedfasta2gff\trepeat\t$region_start\t$region_end\t.\t.\t.\tID=repeat.$counter\n";
	}
}
close IN;
close BED;
close GFF3;
close OUTH;
close OUTS;
close OUTN;
