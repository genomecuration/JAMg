#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/AccompLib");

use Fasta_reader;

my $usage = "\n\nusage: $0 Fasta_file\n\n";

my $fasta_file = $ARGV[0] or die $usage;

main: {
  my $fasta_reader = new Fasta_reader($fasta_file);
  while (my $seq_obj = $fasta_reader->next()) {
    my $accession = $seq_obj->get_accession();
    my $sequence = $seq_obj->get_sequence();

    while ($sequence =~ /([^N]+)/gi) {
      my $start = $-[0] + 1;
      my $end = $+[0];
	  
	  my $start_retrieval = $start;
	  my $end_retrieval = $end;

	  {
		  # include flanking N chars as sanity check
		  if ($start > 1) { $start_retrieval--; }
		  if ($end < length($sequence)) {$end_retrieval++;}
	  }
	  
	  my $len = $end - $start + 1;
	  
	  my $subseq = substr($sequence, $start_retrieval - 1, $end_retrieval - $start_retrieval + 1);
	  
      print "$accession\t$start\t$end\t$len\t$subseq\n";
    }
  }

  exit(0);

}


