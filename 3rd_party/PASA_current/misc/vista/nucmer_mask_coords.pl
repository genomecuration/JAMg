#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/AccompLib");

use Fasta_reader;

my $usage = "\nusage: $0 FastaFile (Query|Reference) CoordsFile\n\n";

my $fasta_file = $ARGV[0] or die $usage;
my $seq_type = $ARGV[1] or die $usage;
my $coords_file = $ARGV[2] or die $usage;

unless ($seq_type =~ /Query|Reference/i) { die $usage; }

main: {

  ## parse the coords file:
  my %acc_to_coords_list;
  {
    open (my $fh, $coords_file) or die $!;
    while (<$fh>) {
		s/\s+$//; # trim trailing ws
		my @x = split (/\t/);
      unless (scalar @x == 9) {
	die "Error, unexpected number of columns.  Should only be 9 fields in nucmer output. Be sure to run show-coords -H -T to generate the out.coords file";
      }
      my ($acc, $end5, $end3) = ($seq_type =~ /Reference/) 
	? ($x[7], $x[0], $x[1])
	  : ($x[8], $x[2], $x[3]);
	  	  
      push (@{$acc_to_coords_list{$acc}}, [$end5, $end3]);
    }
  }

  my $fasta_reader = new Fasta_reader($fasta_file);
  while (my $seq_obj = $fasta_reader->next()) {
    my $acc = $seq_obj->get_accession();
    my $sequence = $seq_obj->get_sequence();

    if (my $coords_aref = $acc_to_coords_list{$acc}) {
      my @chars = split (//, $sequence);
      foreach my $coordpair (@$coords_aref) {
	my ($lend, $rend) = sort {$a<=>$b} @$coordpair;
	for (my $i = $lend; $i <= $rend; $i++) {
	  $chars[$i-1] = 'N';
	}
      }
      $sequence = join ("", @chars);
    }

    $sequence =~ s/(\S{60})/$1\n/g;
    
    print ">$acc\n$sequence\n";
    
  }
  
  exit(0);
}
