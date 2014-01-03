#!/usr/local/bin/perl
#Copyright (c) 2003  by  Mihaela Pertea.

use strict;
use FileHandle;

return 1;

sub splitiso{
    my ($seqs,$beg,$end,$newseqsname)=@_;

    open(F,$seqs);
    open(O,">$newseqsname");

    while(<F>) {
	chomp;
	if($_) {
	    my ($name,$seq)=split;
	
	    my $acgt=($seq =~ tr/ACGTacgt//);
	    my $cg=($seq =~ tr/CGcg//);

	    my $cgperc=$cg*100/$acgt;

	    if($beg <= $cgperc && $cgperc<=$end) { print O $_,"\n"; }

	}
    }

    close(F);
    close(O);
}


