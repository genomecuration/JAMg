#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

my $file = shift; 

die "Please provide a FASTA file to split when N{100,} gaps are found\n" if !$file;
die "File $file not found\n" if !-s $file;

## >fkdljssdflkdsfkl kdsjlfjklsdfdsf
## [ATGCatgc]+

my $orig_sep = $/;
$/ = '>';

open (IN, $file) || die $!;
open (OUT,">$file.split");


while (my $line = <IN>){
	chomp($line);
	next if !$line;
	
	my @data = split("\n", $line);

	my $id = shift(@data);
	my $desc;
	if ($id=~s/^(\S+)\s(.+)/$1/){
		$desc = $2;
	}

	my $sequence = join('',@data);	

	my @subsequences = split(/N{100,}/i,$sequence);
	my $counter = int(0);

	for (my $i=0; $i < scalar(@subsequences); $i++){
		print OUT ">$id.$i";
		# find 
		my $pos = index($sequence,$subsequences[$i],$counter) + 1;

		print OUT " coord=$pos:" . ($pos + length($subsequences[$i]) );
		print OUT " $desc" if $desc;
		print OUT "\n" . $subsequences[$i] . "\n";
		$counter += length($subsequences[$i]);
	}
	

}



close IN;
close OUT;
