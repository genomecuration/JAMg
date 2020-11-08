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
    $string =~ s/U/T/ig;
    $string =~ s/S/G/ig;
    $string =~ s/R/G/ig;
    $string =~ s/Y/T/ig;
    $string =~ s/M/C/ig;
    $string =~ s/K/T/ig;
    $string =~ s/W/T/ig;
    $string =~ s/B/T/ig;
    $string =~ s/D/T/ig;
    $string =~ s/H/T/ig;
    $string =~ s/V/G/ig;
    $string =~ s/X/N/ig;

    return $string;
}
