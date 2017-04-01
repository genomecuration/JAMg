#!/usr/bin/env perl

use strict;
use warnings;

my $hhindex = shift;
my @hhblitsfiles = @ARGV;

die "Provide INDEX and RESULT file from FFINDEX and HHblits\n" unless $hhindex && -s $hhindex && $hhblitsfiles[0]  && -s $hhblitsfiles[0] ;

my %hash;

open (INDEX,$hhindex)||die;

while (my $record = <INDEX>){
	chomp($record); next unless $record;
	my @data = split("\t",$record);
	$hash{$data[0]}=$data[1] if $data[0] && $data[1];
}
close INDEX;

foreach my $hhblitsfile (@hhblitsfiles){
	open (RESULTS, $hhblitsfile)||die;
	while (my $ln = <RESULTS>){
		if ($ln=~/^Query\s+(\S+)\s*/){
			delete($hash{$1});
		}
	}
	close RESULTS;
}

open (OUT,">$hhindex.notdone");
foreach my $id (keys %hash){
	print OUT "$id\t".$hash{$id}."\n";
}
my $count = `wc -l < $hhindex.notdone`;
close OUT;
print "Done!	";
print "See $hhindex.notdone : ".$count if -s "$hhindex.notdone";
print "All queries processed\n" if !-s "$hhindex.notdone";
