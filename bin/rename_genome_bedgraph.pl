#!/usr/bin/env perl

use strict;
use warnings;
my $list = shift;
my $bg = shift||die;
my %hash;
open (IN,$list)||die;
while (my $ln=<IN>){
	chomp($ln);
	my @data = split("\t",$ln);
	$hash{$data[0]} = $data[1] if $data[1];
}
close IN;

my $out_bg_file = "$bg.renamed";
open (BG,$bg)||die;
open (OUT,">$out_bg_file")||die;

my $bedGraphToBigWig_exec = `which bedGraphToBigWig`;
chomp ($bedGraphToBigWig_exec) if $bedGraphToBigWig_exec;


while (my $ln=<BG>){
	my @data = split("\t",$ln);
	if ($hash{$data[0]} && $data[3]){
		$data[0] = $hash{$data[0]};
		print OUT join("\t",@data);
	}
}
close BG;
close OUT;
my $bw_file = $out_bg_file;
$bw_file=~s/\.bg\.renamed$/.bw.renamed/;
system("$bedGraphToBigWig_exec $out_bg_file ".$ENV{'GENOME_PATH'}.".fai $bw_file") if $bedGraphToBigWig_exec && $ENV{'GENOME_PATH'};
