#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Overlap_piler;


my $usage = "usage: $0 btab_file min_perID[=98]\n\n";

my $btab_file = $ARGV[0] or die $usage;
my $min_per_id = $ARGV[1] || 98;

my %data;

open (my $fh, $btab_file) or die "Error, cannot open file $btab_file";
while (<$fh>) {
	chomp;
	my @x = split (/\t/);
	my ($contig_acc, $lend, $rend, $per_id) = ($x[0], $x[6], $x[7], $x[10]);

	if ($per_id >= $min_per_id) {
		($lend, $rend) = sort {$a<=>$b} ($lend, $rend);

		push (@{$data{$contig_acc}}, [$lend, $rend]);
	}
}
close $fh;

foreach my $contig (keys %data) {
	my @coord_pairs = @{$data{$contig}};

	my @clusters = &Overlap_piler::simple_coordsets_collapser(@coord_pairs);

	foreach my $cluster (@clusters) {
		my ($lend, $rend) = @$cluster;

		print "$contig\t$lend-$rend\n";
		
	}

}

exit(0);


