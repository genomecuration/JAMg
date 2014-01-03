#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;
use BHStats;

use Data::Dumper;

my $usage = "usage: $0 sequences.fasta KmerSize NumTopKmers\n\n";

my $seqFile = $ARGV[0] or die $usage;
my $kmerSize = $ARGV[1] or die $usage;
my $numTopKmers = $ARGV[2] or die $usage;


my $num_kmers_per_round = 100;

## find all Kmers on top strand
my %Kmers;

{
	my $fasta_reader = new Fasta_reader($seqFile);
	while (my $seq_obj = $fasta_reader->next()) {
		
		my $sequence = uc $seq_obj->get_sequence();
		
		for (my $i = 0; $i <= length($sequence) - $kmerSize; $i++) {
			my $kmer = substr($sequence, $i, $kmerSize);
			$Kmers{$kmer}++;
		}
	}
	$fasta_reader->finish();
}


## find position abundance.

my %Kmer_to_deviation;


my @K = sort keys %Kmers;
for (my $i = 0; $i <=$#K; $i+= $num_kmers_per_round) {
	my $i2 = $i + $num_kmers_per_round-1;
	if ($i2 > $#K) {
		$i2 = $#K;
	}
	my @kmers_to_analyze = @K[$i..$i2];
	
	my %kmer_to_position_abundance = &search_Kmer_position_abundance(@kmers_to_analyze);
	
	## report
	foreach my $kmer (keys %kmer_to_position_abundance) {
		my @vals;
		#print "$kmer";
		my $seq_length = $#{$kmer_to_position_abundance{$kmer}} + 1;
		for (my $i = 0; $i < $seq_length; $i++) {
			my $count = $kmer_to_position_abundance{$kmer}->[$i] || 0;
			#print "\t$count";
			push (@vals, $count);
		}
		#print "\n";
		
		#print join (" ", @vals) . "\n";
		
		my $median = &BHStats::median(@vals);
		my $max = &BHStats::max(@vals);
		
		my $delta = $max - $median;
		$Kmer_to_deviation{$kmer} = $delta;
		
		print STDERR "Kmer: $kmer, Delta: $delta\n" #, Max: $max, Median: $median\n";
	}
	
}

## Report the top most deviant Kmers
my @deviant_Kmers = reverse sort {$Kmer_to_deviation{$a}<=>$Kmer_to_deviation{$b}} keys %Kmer_to_deviation;

{
	open (my $ofh, ">Kmer.deltas.$$") or die $!;
	foreach my $kmer (@deviant_Kmers) {
		print $ofh join ("\t", $kmer, $Kmer_to_deviation{$kmer}) . "\n";
	}
	close $ofh;
}


if ($numTopKmers > scalar(@deviant_Kmers)) {
	$numTopKmers = scalar(@deviant_Kmers);
}

my @topKmers = @deviant_Kmers[0..$numTopKmers-1];
my %positionAbundance = &search_Kmer_position_abundance(@topKmers);
&report_Kmer_abundance(%positionAbundance);




exit(0);

####
sub search_Kmer_position_abundance {
	my (@kmers_to_analyze) = @_;

	my %Kmers = map { + $_ => 1 } @kmers_to_analyze;
	
		
	my %kmer_to_position_abundance;	
	
	my $seq_length;
	my $fasta_reader = new Fasta_reader($seqFile);
	while (my $seq_obj = $fasta_reader->next()) {
		
		my $sequence = uc $seq_obj->get_sequence();
		
		if ($seq_length && length($sequence) != $seq_length) {
			die "Error, inconsistent sequence lengths: $seq_length vs. " . length($sequence) . " ";
		}
		elsif (! $seq_length) {
			$seq_length = length($sequence);
			## init kmer positions
			for (my $i = 0; $i <= $seq_length - $kmerSize; $i++) {
				foreach my $kmer (keys %Kmers) {
					$kmer_to_position_abundance{$kmer}->[$i] = 0;
				}
			}
		}
		
		## iterate Kmers in sequence
		for (my $i = 0; $i <= $seq_length - $kmerSize; $i++) {
			my $kmer = substr($sequence, $i, $kmerSize);
			#print "$kmer\n";
			if ($Kmers{$kmer}) {
				$kmer_to_position_abundance{$kmer}->[$i]++;
			}
		}
	}
	
	return(%kmer_to_position_abundance);
}


####
sub report_Kmer_abundance {
	my (%kmer_to_position_abundance) = @_;
	
	foreach my $kmer (keys %kmer_to_position_abundance) {
		print "$kmer";
		my $seq_length = $#{$kmer_to_position_abundance{$kmer}} + 1;
		for (my $i = 0; $i < $seq_length; $i++) {
			my $count = $kmer_to_position_abundance{$kmer}->[$i];
			print "\t$count";
		}
		print "\n";
	}
	return;
}
