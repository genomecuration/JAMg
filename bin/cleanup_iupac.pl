#!/usr/bin/perl -w

use strict;

my $file = shift;
$file = &cleanupiupac($file);


##########################################
sub cleanupiupac ($) {
    my $infile  = shift;
	die "No input!\n" if !$infile || !-s $infile;
    my $outfile = "$infile.x";
    die "Outfile already exists\n" if ( -s $outfile );
    my $orig_sep = $/;
    $/ = '>';
    open (IN,$infile);
    open (OUT,">".$outfile);
    while (my $record=<IN>){
       chomp($record);
       next unless $record;
       my @lines = split("\n",$record);
       my $id = shift(@lines);
       my $seq=uc(join("",@lines));
       $seq=~s/\s+//g;
       next unless $seq;
	$seq = &iupac_replace($seq);
	$seq = &wrap_text($seq);
	print OUT ">$id\n$seq\n";
    }
    $/ = $orig_sep;    
    return ($outfile);
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
    return $string;
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

