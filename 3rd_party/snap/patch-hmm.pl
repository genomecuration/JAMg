#!/usr/bin/perl
use strict; use warnings;

die "usage: $0 <hmm> <model> <model> <etc>" unless @ARGV > 1;
my ($hmm, @file) = @ARGV;


my %model;
foreach my $file (@file) {
	my $contents = `cat $file`;
	my ($type) = $contents =~ /^(\S+)/;
	$model{$type} = $contents;
}

open(IN, $ARGV[0]) or die;
while (<IN>) {
	print;
	last if /SEQUENCE_MODELS/;
}
while (<IN>) {
	if (/^(\S+)/) {
		my $type = $1;
		if (defined $model{$type}) {
			print $model{$type};
			while (<IN>) {last unless /\S/}
		} else {
			print;
			while (<IN>) {
				print;
				last unless /\S/;
			}
		}
	}
	print;
}