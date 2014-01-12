#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;


## we delete all files we don't need in this directory. Be careful in case users try running it somewhere else, outside this dir.
chdir $FindBin::Bin or die "error, cannot cd to $FindBin::Bin";



my @files_to_keep = qw (
						top100k.Left.fq
						top100k.Left.fq.gz
						top100k.Right.fq
						top100k.Right.fq.gz
						top100k.bam
						top100k.genome
                        Trinity.fasta
						top100k.genes.gff3
						cleanme.pl
						run_sample_pipeline.pl
						alignAssembly.config
                        cufflinks.transcripts.gtf
						
                        __run_test_use_cufflinks.pl

						);

my %keep = map { + $_ => 1 } @files_to_keep;

if (-d "gmap_db_dir") {
    `rm -rf gmap_db_dir`;
}
if (-d "blat_out_dir") {
    `rm -rf blat_out_dir/`;
}

foreach my $file (<*>) {
	
	if (-f $file && ! $keep{$file}) {
		print STDERR "-removing file: $file\n";
		unlink($file);
	}
}

`rm -rf pasa_run.*.log`;

exit(0);
