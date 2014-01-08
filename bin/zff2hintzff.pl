#!/usr/bin/env perl

use strict;
use warnings;
my $file = shift||die "<file> <out directory> <command ADJ|SET> <amount>\n Defaults to ./evidence ADJ +10\n";
die "Didn't find $file\n<file> <out directory> <command ADJ|SET> <amount>\n Defaults to ./evidence ADJ +10\n" unless -s $file;
my $directory = shift;
$directory = 'evidence' unless $directory;
my $command = shift;
$command = 'ADJ' unless $command;
my $amount = shift;
$amount = '10' unless $amount;
if ($amount =~/^\d+/){
	$amount = '+'.$amount;
}
mkdir $directory unless -d $directory;
open (IN,$file);

my $master_id;
while (my $ln=<IN>){
	if ($ln=~/^>(\S+)/){
		close OUT if $master_id;
		$master_id = $1;
		open (OUT,">$directory/$master_id.snap.evidence") || die "Can't create output $directory/$master_id.snap.evidence";
		print OUT $ln;
	}else{
		chomp($ln);
		my @data=split("\t",$ln);
		my $start = $data[1] > $data[2] ? $data[2] : $data[1];
		my $end   = $data[1] > $data[2] ? $data[1] : $data[2];
		my $strand= $data[1] > $data[2] ? '-' : '+';
		print OUT "Coding\t$start\t$end\t$strand\t$amount\t.\t.\t.\t$command\n";

	}
}
close IN;
close OUT;
print "Done, see $directory/\n";
