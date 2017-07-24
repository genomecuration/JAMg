#!/usr/bin/env perl

use strict;
use warnings;

my $file = shift;
die "Please give a FASTA file to adjust single Ns\n" unless $file && -s $file;

my $orig_sep = $/;
$/ = '>';

my @nucs = qw/A T C G/;
my $random_nuc = $nucs [rand @nucs];

open (IN,$file) ||die;
open (OUT,">$file.x") ||die;
while (my $record = <IN>){
	chomp($record);
	next unless $record;

	my @lines = split("\n",$record);
	my $id = shift(@lines);
	my $seq = join('',@lines);
	$seq=~s/\s+//g;
	while ($seq=~s/([ATCG])N{1}([ATCG])/$1$random_nuc$2/g){
		$random_nuc = $nucs [rand @nucs];
	}
	$seq = &wrap_text($seq);
	print OUT ">$id\n" . $seq;

}

close IN;
close OUT;


#########################################

sub wrap_text() {
 my $string      = shift;
 my $wrap_length = shift;
 $wrap_length = 120 if !$wrap_length;
 $string =~ s/(.{0,$wrap_length})/$1\n/g;
 $string =~ s/\n{2,}/\n/;
 return $string;
}

