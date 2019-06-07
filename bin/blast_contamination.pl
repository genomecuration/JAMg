#!/usr/bin/env perl

use strict;
use warnings;

my $tabfile = shift;
die unless $tabfile && -s $tabfile;

open (IN,$tabfile);
open (OUT1,">$tabfile.bacterial");
open (OUT2,">$tabfile.dunno");
my $size_total1 = int(0);
my $size_total2 = int(0);
my $header = <IN>;
while (my $ln=<IN>){
	chomp($ln); next unless $ln;
	my @data = split("\t",$ln);
	my $scaffold = shift (@data);
	next unless $data[4];
	next if ($data[4]=~/Insecta|Mammalia|Aves|Amphibia/ );
	if ($data[4] =~/bacteria|Bacilli|Clostridia|Chlamydiia/i){
		$size_total1 += $data[1];
		print OUT1 $ln."\n";
	}else{
		next if $data[5] && ($data[5]=~/Insecta|Mammalia|Aves|Amphibia/);
		$size_total2 += $data[1];
		print OUT2 $ln."\n";
	}
}

close IN;
close OUT1;
close OUT2;

print "$tabfile: Found $size_total1 bp of putative bacterial data\n";
print "$tabfile: Found $size_total2 bp of suspect non-Insect/non-Mammalian data I couldn't decide\n";
