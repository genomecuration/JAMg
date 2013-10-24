#!/usr/bin/perl
use strict; use warnings;
use IK;
use FAlite;
use DataBrowser;

die "usage: $0 <fasta> <k-mer> <Intron|Inter>" unless @ARGV == 3;
my ($file, $K, $type) = @ARGV;
die unless $ARGV[2] =~ /Intron|Inter/;


my $model = IK::blank_table($K, 1);
open(IN, $file) or die;
my $fasta = new FAlite(\*IN);
while (my $entry = $fasta->nextEntry) {
	my $dna = uc $entry->seq;
	for (my $i = 0; $i < length($dna) -$K +1; $i++) {
		$model->{substr($dna, $i, $K)}++;
	}
}
close IN;

my $table = IK::blank_table($K -1);
my @alph = qw(A C G T);
printf "%s LUT %d %d 4 0 0.000\n", $type, $K, $K-1;

foreach my $kmer (sort keys %$table) {
	my $total = 0;
	foreach my $nt (@alph) {
		$total += $model->{"$kmer$nt"};
	}
	
	foreach my $nt (@alph) {
		my $obs = $model->{"$kmer$nt"} / $total;
		printf "\t%.3f", log($obs/0.25)/log(2);
	}
	print "\n";
}
