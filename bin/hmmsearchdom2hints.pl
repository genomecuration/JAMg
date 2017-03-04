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
	my $strand = $ln=~/REVERSE SENSE/ ? '-' : '+';
	if ($ln=~/^>(\S+)\s\[(\d+)\s\-\s(\d+)\]/){
		my $id = $1;
		my $start = ($strand eq '+') ? $2 : $3;
		my $end = ($strand eq '+') ? $3 : $2;
		$hash{$id} = "$strand ".$start."\t".$end;
	}
}
close EXON;

print "Processing TABLE HMMSEARCH file $tablefile\n";
open (IN, $tablefile)||die $!;


open (GFF3,">$tablefile.gff3") ||die;
open (HINTS,">$tablefile.hints") ||die;


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
	my $strand;
	if ($start_end=~s/^(\S)\s//){
		$strand = $1;
	}else{next;}
	$id =~s/_\d+$//;

	print GFF3 "$id\thmmsearch\tprotein_match\t".$start_end
	."\t$score\t$strand\t.\tID=$hit;Name=$name\n";

	print HINTS "$id\thmmsearch\tCDSpart\t".$start_end
	."\t$score\t$strand\t.\tsrc=HU;grp=$hit;prio=2\n";
}
close IN;
close GFF3;
close HINTS;
