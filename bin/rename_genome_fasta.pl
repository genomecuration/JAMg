#!/usr/bin/env perl

use strict;
use warnings;
my $list = shift;
my $fasta = shift||die ("Provide a tab delimited rename list and FASTA file\n");
my %hash;
open (IN,$list)||die;
while (my $ln=<IN>){
	chomp($ln);
	my @data = split("\t",$ln);
	$hash{$data[0]} = $data[1] if $data[1];
}
close IN;


my $orig_sep = $/;
$/ = '>';
open (FASTA,$fasta)||die;
open (OUT,">$fasta.renamed")||die;

while (my $record=<FASTA>){
	chomp($record);
	next unless $record;
	my @lines = split($orig_sep,$record);
	my $id = shift(@lines);
	my $seq = join('',@lines);
	my $new_id = $hash{$id};
	if (!$new_id){
		$new_id = $id;
		warn "$id no renaming entry\n";
	}
	print OUT ">$new_id\n".&wrap_text($seq)."\n";
}
close FASTA;
close OUT;


#################
sub wrap_text() {
# there seems to be an error sometimes.
 my $string      = shift;
 my $wrap_length = shift;
 $wrap_length = 80 if !$wrap_length;
 $string =~ s/(.{0,$wrap_length})/$1\n/g;
 $string =~ s/\n{2,}/\n/g;
 $string =~ s/\s+$//;
 return $string;
}

