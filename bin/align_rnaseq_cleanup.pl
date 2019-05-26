#!/usr/bin/env perl

use strict;
use warnings;

my @logfiles = glob("gsnap.*log");
exit(0) if !$logfiles[0];

print "These can be deleted with these commands:\n";
foreach my $log (sort @logfiles){
	my $bg_file = $log;
	$bg_file=~s/\.log$/.concordant_uniq.coverage.bg/;
	if (!-s $bg_file){
		$bg_file=~s/.concordant_uniq.coverage.bg$/.unpaired_uniq.coverage.bg/;
	}
	if (-s $log && -s $bg_file){
		open (LOG, $log);
		my @lines =  <LOG>;
		close LOG;
		foreach my $ln (@lines){
			if ($ln && $ln =~ /^GSNAP Completed/){
				if ($log =~/^gsnap\.(.+)_vs_/){
					my $input_readset = $1;
					if (glob($input_readset."*fastq")){
						print "\trm -f ".$input_readset."*.fastq\n";
					}
				}
			}
		}
	}
}
