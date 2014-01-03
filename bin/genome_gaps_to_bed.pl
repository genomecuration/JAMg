#!/usr/bin/env perl

=pod

=head1 NAME

 genome_gaps_to_bed
 
=head1 USAGE

genome_gaps_to_bed.pl <fastafile> [sequence ID]

 Create a BED file that shows where that gaps (Ns) are in a FASTA file.
 Requires Kent's src code

 git clone git://genome-source.cse.ucsc.edu/kent.git
 
=head1 AUTHORS

 Alexie Papanicolaou

        CSIRO Ecosystem Sciences
        alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

Copyright 2012-2014 the Commonwealth Scientific and Industrial Research Organization. 
This software is released under the Mozilla Public License v.2.

It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.mozilla.org/MPL/2.0.
 
=cut

use strict;
use warnings;
use Pod::Usage;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";


my $cwd = `pwd`;
chomp($cwd);
my ($wig_to_bw_exec) = &check_program('wigToBigWig');

my $genome = shift || pod2usage("Need a FASTA file\n");
my $target_id = shift;

my $wigfile = $target_id ? "$genome.gap.$target_id.wig" : "$genome.gap.wig";
my $chrsizes =
  $target_id ? "$genome.chrom.$target_id.sizes" : "$genome.chrom.sizes";
unless ( -s $wigfile || -s "$wigfile.bw" ) {
	open( WIG,   ">$wigfile" );
	open( STATS, ">$chrsizes" );
	print WIG
"track name='gaps' description='gaps' type='wiggle_0' color=50,50,200 yLineMark=0.0 yLineOnOff=on visibility=full maxHeightPixels=40:40:2 priority=1\n";
	my $orig_sep = $/;
	$/ = ">";
	open( IN, $genome ) || die;
	<IN>;

	while ( my $record = <IN> ) {
		my @lines = split( "\n", $record );
		my $id = shift @lines;
		$id =~ /^(\S+)/ || next;
		$id = $1;
		next if ( $target_id && $target_id ne $id );
		my %hash;
		chomp(@lines);
		my $seq = join( '', @lines );
		$seq =~ s/>$//;
		$seq =~ s/\s+//g;
		my $length = length($seq);
		print STATS "$id\t$length\n";
		my $offset = int(0);
		my $gap_search = index( $seq, 'N', $offset );

		while ( $gap_search != -1 ) {
			$hash{$gap_search} = 1;
			$offset = $gap_search + 1;
			$gap_search = index( $seq, 'N', $offset );
		}
		print WIG "fixedStep  chrom=" . $id . "  start=1  step=1\n";
		for ( my $i = 0 ; $i < $length ; $i++ ) {
			if ( $hash{$i} ) {
				print WIG "1\n";
			}
			else {
				print WIG "0\n";
			}

		}
	}
	close WIG;
	close STATS;
	$/ = $orig_sep;

}

&process_cmd("$wig_to_bw_exec $wigfile $chrsizes $wigfile.bw")
  unless -s "$wigfile.bw";
unlink($wigfile) if -s "$wigfile.bw";

sub process_cmd {
	my ( $cmd, $dir ) = @_;
	chdir($dir) if $dir;
	print "$cmd\n";
	my $ret = system($cmd);
	if ( $ret && $ret != 256 ) {
		chdir($cwd) if $dir;
		die "Error, cmd died with ret $ret\n";
	}
	chdir($cwd) if $dir;
	return;
}

###
sub check_program() {
	my @paths;
	foreach my $prog (@_) {
		my $path = `which $prog`;
		pod2usage
		  "Error, path to a required program ($prog) cannot be found\n\n"
		  unless $path =~ /^\//;
		chomp($path);
		$path = readlink($path) if -l $path;
		print "Found $path\n";
		push( @paths, $path );
	}
	return @paths;
}

