#!/usr/bin/env perl

use strict;
use warnings;

my $file = shift;
die "Give me a GFF" unless $file && -s $file;

open (GFF,$file) || die $!;
while (my $ln=<GFF>){
	next unless $ln;
	my @data = split("\t",$ln);
# from
#Pd01    CNAG    intron  22798574        22798934        .       -       .       ID=Prudul26A000004T1.intron1;Parent=Prudul26A000004T1
# to
#>Prudul26A000004T1.intron1 Pd01:22798934..22798574
	next unless $data[8];
	if ($data[2] eq 'intron' && $data[8]=~/ID=([^;]+)/){
		if ($data[6] eq '+'){
			print ">$1 ".$data[0].':'.$data[3].'..'.$data[4]."\n";
		}elsif ($data[6] eq '-'){
			print ">$1 ".$data[0].':'.$data[4].'..'.$data[3]."\n";
		}
	}

}

close GFF;
