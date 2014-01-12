#!/usr/bin/env perl

use strict;
use warnings;

main: {

	

	## Purge the current sample mysql database if it exists from a previous run of this pipeline. Start fresh.
	
	my $cmd = "../scripts/drop_mysql_db_if_exists.dbi -c alignAssembly.config";
	&process_cmd($cmd);
	
	
	## run RNA-Seq inchworm assembly alignment to genome with PASA alignment assembly:
	$cmd = "../scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -C -R -g top100k.genome -t Trinity.fasta --transcribed_is_aligned_orient --stringent_alignment_overlap 30 --ALIGNERS blat ";
	&process_cmd($cmd);
	
	
	## Compare to existing annotations and provide annotation updates:
	$cmd = "../scripts/Launch_PASA_pipeline.pl -c alignAssembly.config -g top100k.genome -t Trinity.fasta  -A -L --annots_gff3 top100k.genes.gff3";
	&process_cmd($cmd);



	print "\n\n\nDone.\n\n";
	exit(0);
	
}


####
sub process_cmd {
	my ($cmd) = @_;

	print "CMD: $cmd\n";
	
	my $ret = system($cmd);

	if ($ret) {
		die "Error, cmd: $cmd died with ret ($ret)";
	}

	return;
}


