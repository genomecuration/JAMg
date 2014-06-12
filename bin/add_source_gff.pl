#!/usr/bin/env perl

use strict;
use warnings;

my $file = shift||die ("Give a GFF file and the source tag\n");
my $source = shift ||die;
open (IN,$file);
open (OUT,">$file.out");


while (my $ln=<IN>){
	my @data = split("\t",$ln);
	if ($data[1]){
		$data[1] = $source;
		$ln = join("\t",@data);
	}
	print OUT $ln;
}

close IN;
close OUT;
