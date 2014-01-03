#!/usr/local/bin/perl

#Copyright (c) 2003  by  Mihaela Pertea.

# this library takes the files seqs and exons.dat and extracts from them
# two output files:
# one with the orfs
# the other one with introns --- commented

use strict;
use FileHandle;

return 1;

sub genorf{ 
    my ($f,$g) = @_;

    my $f1="orfs";
#    my $f2="introns";

    my (%isexon, %isanum, %isintron, $end, $start, $start1, $end1, $anum);

    open(F,$f) or die "ERROR 131: Couldn't open file $f for reading!\n";

    while(<F>){
	chomp;
	($anum,$start,$end) = split;

	if($isanum{$anum}) {
	    $isexon{$anum}.="$start $end ";
	    $start1=$start-1;
	    $isintron{$anum}.="$end1 $start1 ";	
	}

	else {
	    if($anum ne "") {
		$isanum{$anum}=1;
		$isexon{$anum}="$start $end ";
	    }
	}

	$end1=$end+1;
	
    }
    close(F);
	

    open(G,$g) or die "ERROR 132: Couldn't open file $g for reading!\n";

    open(F1,">$f1") or die "ERROR 133: Couldn't open file $f1 for writing!\n";
#    open(F2,">$f2");

    my $seq;

    while(<G>){
	chomp;
	($anum,$seq)=split;
    
	my @e=split(/\s+/,$isexon{$anum});
	my @i=split(/\s+/,$isintron{$anum});

	my $orf="";

	for(my $k=0;$k<=$#e;$k+=2){
	    $orf.=substr($seq,$e[$k]-1,$e[$k+1]-$e[$k]+1);
	}

	# my comment begin
	# if you want the whole coding gene use:
	# print F1 "$anum $orf\n";

	# otherwise: skip start condon and last codon:
	
	my $l=length($orf);
	my $neworf=substr($orf,3,$l-6);

	print F1 "$anum $neworf\n";

	# my comment end

#	for(my $k=0;$k<$#i;$k+=2){
#	    my $intron=substr($seq,$i[$k]-1,$i[$k+1]-$i[$k]+1);
#	    print F2 "$anum $intron\n";
#	}
    }

    close(G);
    
    F1->autoflush(1);
    close(F1);

    STDOUT->autoflush;
}

sub gencodorf {
    my ($f,$g,$f1)=@_;

    open(F,$f);

    my %isexon;
    while(<F>){
	chomp;
	if($_) {
	    my ($anum,$start,$end) = split;
	    if($start<$end) { $isexon{$anum}.="$start $end ";}
	}
    }
    close(F);
	
    open(G,$g);

    open(F1,">$f1");
    
    while(<G>){
	chop;
	my ($anum,$seq)=split;
	my @e=split(/\s+/,$isexon{$anum});
	my $orf="";

	for(my $k=0;$k<=$#e;$k+=2){
	    $orf.=substr($seq,$e[$k]-1,$e[$k+1]-$e[$k]+1);
	}

	chop($orf);chop($orf);chop($orf); # to skip stop codon
	print F1 "$anum $orf\n";

    }

    close(G);
}


sub gennoncod {
    my ($exons,$seqs,$outpf)=@_;

    open(F,$exons);
    my %noncod;
    my $lastend;
    my $name;
    
    while(<F>) {
	chomp;
	if($_) {
            my ($beg,$end);
	    ($name,$beg,$end)=split;
	    if(!$noncod{$name}) { 
		$noncod{$name}=$beg;
	    }
	    else {
		$noncod{$name}.=" ".$lastend." ".$beg;
	    }
	    $lastend=$end;
	}
	else {
	    $noncod{$name}.=" ".$lastend;
	}
    }

    close(F);

    open(F,$seqs);
    open(O,">$outpf");
   
    while(<F>) {
    
        my $seq;
	($name,$seq)=split;

	$seq =~ tr/\n//d;
	$seq =~ tr/\r//d;

	my $len=length($seq);

	my @c=split(/\s+/,$noncod{$name});

	my $n=0;

	print O $name,"_",++$n," ",substr($seq,0,$c[0]-1),"\n";
	for(my $i=1;$i<$#c-1;$i+=2) {
	    print O $name,"_",++$n," ",substr($seq,$c[$i],$c[$i+1]-1-$c[$i]),"\n";
	}
	print O $name,"_",++$n," ",substr($seq,$c[$#c]),"\n";

    }
    close(F);
}
