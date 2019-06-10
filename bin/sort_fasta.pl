#!/usr/bin/env perl

use strict;
use warnings;

my $orig_sep = $/;
$/ = '>';

my %hash;

while (my $record=<STDIN>){
	chomp($record);next unless $record;
	my @lines = split($orig_sep,$record);
	my $id = shift(@lines);
	my $seq = join('',@lines);
	next unless $seq;
	$hash{$id}{'length'} = length($seq);
	$hash{$id}{'seq'} = $seq;
}

$/ = $orig_sep;

foreach my $id (sort {$hash{$b}{'length'} <=> $hash{$a}{'length'}} keys %hash){
	print ">".$id.$/.&wrap_text($hash{$id}{'seq'}).$/;

}



#########################################

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

