#!/usr/bin/perl

# this program takes a multifasta file and computes the predicted 
# genes for each fasta sequence in the file

use strict;

my ($program,$file,$traindir,$options)=@ARGV; 

if($#ARGV<2) {
    print STDERR "Usage: glimmhmm.pl <glimmerhmm_program> <fasta_file> <train_dir> <options>\n";
    exit;
}

open(F,$file)
    || die "can't open $file.\n";

$/=">";
my $first=<F>;

$first=1;

while(<F>){

    chomp;

    my ($name)=/^(\S+)\s+/;
    my $pos=index($_,"\n");
    my $seq=substr($_,$pos+1);
    
    die "ERROR: Wrong FASTA format: .$_." if (!$name || !$seq);
    
    $seq =~ tr/\n//d;
    $seq =~ tr/\r//d;
    

    my $temp="glimmhmm.temp.fasta";

    open(T,">$temp");
    print T ">$name\n$seq";
    close(T);

    print STDERR "Process $name\n";

    my $out=`$program $temp $traindir $options`;

    if(!($options =~ s/-g/-g/)) { 
	my @args = ("echo", "\n\n\n######################### Analyzing contig $name #########################\n\n\n");
	system(@args) == 0
	or die "system @args failed: $?";
	print $out,"\n";
    }
    else {
	if($first) { print $out; $first=0;}
	else {
	    my @line=split(/\n/,$out);
	    for(my $i=1;$i<=$#line;$i++) {
		print $line[$i],"\n";
	    }
	}
    }

    unlink($temp);
    print STDERR "Done $name\n";

}
close(F);
   
