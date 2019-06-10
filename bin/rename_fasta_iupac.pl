#!/usr/bin/env perl

use strict;
use warnings;

my $code = shift;
$code ='' if !$code;

my $orig_sep = $/;
$/ = ">";
my $counter =1;
while (my $record=<STDIN>){
	chomp($record);next unless $record;
	my @lines = split($orig_sep,$record);
	my $id = shift (@lines);
	$id = ">".$code.$counter;
	my $seq = join("",@lines);
	$seq = &iupac_replace($seq);
	print $id.$orig_sep.&wrap_text($seq).$orig_sep;
	$counter++;
}

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

sub iupac_replace($) {
    my $string = shift;
    $string =~ s/\s+//g;
    $string =~ s/S/G/g;
    $string =~ s/R/G/g;
    $string =~ s/Y/T/g;
    $string =~ s/M/C/g;
    $string =~ s/K/T/g;
    $string =~ s/W/T/g;
    $string =~ s/B/T/g;
    $string =~ s/D/T/g;
    $string =~ s/H/T/g;
    $string =~ s/V/G/g;
    $string =~ s/X/N/g;

    $string =~ s/s/G/g;
    $string =~ s/r/G/g;
    $string =~ s/y/T/g;
    $string =~ s/m/C/g;
    $string =~ s/k/T/g;
    $string =~ s/w/T/g;
    $string =~ s/b/T/g;
    $string =~ s/d/T/g;
    $string =~ s/h/T/g;
    $string =~ s/v/G/g;
    $string =~ s/x/N/g;

    return $string;
}
