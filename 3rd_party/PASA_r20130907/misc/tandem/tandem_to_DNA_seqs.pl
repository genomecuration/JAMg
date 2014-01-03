#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use CdbTools;

my %contigs_to_tandem_seqs;

my $usage = "usage: $0 tandem_results genomeFasta outCoreName\n\n";

my $tandem_results = $ARGV[0] or die $usage;
my $genomeFasta = $ARGV[1] or die $usage;
my $out_core_name = $ARGV[2] or die $usage;

main: {
	
	open (my $fh, $tandem_results) or die "Error, cannot open file $tandem_results";
	while (<$fh>) {
		chomp;
		my $line = $_;
		if (/^\#Tandem/) {
			my ($tandem, $contig_info, $tandem_info) = split (/\t/);
			
			my @contig_stuff = split (/\s+/, $contig_info);
			my $contig = $contig_stuff[0];
			
			$tandem_info =~ s/,//g;
			my @x = split (/ /, $tandem_info);
			my ($coords_range, $period_len) = ($x[1], $x[3]);
			
			push (@{$contigs_to_tandem_seqs{$contig}}, { 
				info => $line,
				range => $coords_range,
				period_len => $period_len,
				line => $line,
			} );

		}
	}
	close $fh;

	my $tandem_counter = 0;
	

	open (my $regions_fh, ">$out_core_name.tandem_regions") or die $!;
	open (my $periods_fh, ">$out_core_name.period_rep") or die $!;

	foreach my $contig (keys %contigs_to_tandem_seqs) {
	   
		my $genome_seq = &cdbyank_linear($contig, $genomeFasta);

		foreach my $tandem_href (@{$contigs_to_tandem_seqs{$contig}}) {
			
			$tandem_counter++;
			
			my ($info, $range, $period_len, $line) = ($tandem_href->{info}, $tandem_href->{range}, $tandem_href->{period_len}, $tandem_href->{line});

			my ($lend, $rend) = split (/-/, $range);

			my $tandem_region = substr($genome_seq, $lend - 1, $rend - $lend + 1);
			my $period_seq = substr($tandem_region, 0, $period_len);

			print $regions_fh ">T.$tandem_counter.region $line\n$tandem_region\n";

			print $periods_fh ">T.$tandem_counter.period $line\n$period_seq\n";
		}
	}

	close $regions_fh;
	close $periods_fh;
	
	print "Done. See output files $out_core_name.\*\n";

	exit(0);
}
