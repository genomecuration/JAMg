#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;


my $tablefile = shift;
my $exonfile = shift;

my $score_cutoff = 10.0;


die "Please provide the TABLE output of HHMSEARCH\n" unless $tablefile && -s $tablefile;
die "Please provide the input to HHMSEARCH\n" unless $exonfile && -s $exonfile;

my %hash;

open (EXON,$exonfile);
print "Building lookup from $exonfile\n";

while (my $ln = <EXON>){
	next unless $ln=~/^>/;
	my $rev = $ln=~/REVERSE SENSE/ ? 1 : 0;
	if ($ln=~/^>(\S+)\s\[(\d+)\s\-\s(\d+)\]/){
		my $id = $1;
		my $start = $rev ? $2 : $3;
		my $end = $rev ? $3 : $2;
		$hash{$id} = $start."\t".$end;
	}
}
close EXON;

print "Processing TABLE HMMSEARCH file $tablefile\n";
open (IN, $tablefile)||die $!;


open (OUT,">$tablefile.gff3") ||die;


while (my $ln=<IN>){
	next if $ln=~/^#/;
	chomp($ln);
	my @data = split(/\s+/,$ln);
	next unless $data[5];
	next unless $data[5] >= $score_cutoff;
	my $name = $data[0];
	my $hit = $data[1];
	$hit = $name if $hit eq '-';
	my $score = $data[5];
	my $id = $data[2] || next;
	my $start_end = $hash{$id} || next;
	$id =~s/_\d+$//;

	print OUT "$id\thmmsearch\tprotein_match\t".$start_end
	."\t$score\t.\t.\tID=$hit;Name=$name\n";
}
close IN;
close OUT;
