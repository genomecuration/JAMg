#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
my $fasta_file = shift;
my $ncbi_contam_file = shift;


die "FASTA and ID file\n" unless $fasta_file && -s $fasta_file && $ncbi_contam_file && -s $ncbi_contam_file;
my (%hash);
my $hash_ref =\%hash;
#scaffold_1159   299557  104372..104416,181894..181935,251543..251577    adaptor:NGB00843.1
open (NCBI,$ncbi_contam_file);
while (my $ln = <NCBI>){
	chomp($ln);
	next unless $ln=~/\t/;
	my @data = split("\t",$ln);
	next unless $data[3];
	my @contams = split(',',$data[2]);
	foreach my $contam (@contams){
		push(@{$hash{$data[0]}},$contam);
	}
}
close NCBI;

#die Dumper \%hash;

$/='>';

open (FASTA,$fasta_file);
open (OUT,">$fasta_file.x");
while (my $record=<FASTA>){
	chomp($record);
	next unless $record;
	my @lines = split("\n",$record);
	my $id = shift(@lines);
	$id=~s/^\s+//g;
	$id=~s/\s+$//g;
	my $seq = join('',@lines);
	if ($hash_ref->{$id}){
		my @seq_array = split('',$seq);
		foreach my $gap (@{$hash{$id}}){
			if ($gap =~/^(\d+)\.\.(\d+)$/){
				my $start = $1;
				my $end = $2;
				next unless $start && $end;
print "Changing $id for $start to $end\n";
				# count from 0
				for (my $i=$start-1;$i<$end;$i++){
					$seq_array[$i] = 'N';
				}
			}
		}
		$seq = join('',@seq_array);
	}
	print OUT ">$id\n".$seq."\n";
	
}
close FASTA;
close OUT;
