#!/usr/local/bin/perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;

my $usage = "usage: $0 fastaFile\n\n";

my $fasta_file = $ARGV[0] or die $usage;

my $fasta_reader = new Fasta_reader($fasta_file);

while (my $seq_obj = $fasta_reader->next()) {
    my $accession = $seq_obj->get_accession();

    my $sequence = $seq_obj->get_sequence();

    if ($accession =~ /^gi\|/) {
        my @acc_parts = split (/\|/, $accession);
        my $acc = $acc_parts[3];
        $acc =~ s/\.\d+$//;
        $accession = $acc;
    }

    # put sequence back in fasta format
    $sequence =~ s/(\S{60})/$1\n/g;
    chomp $sequence;
    
    unless ($accession) {
        die "Error, no accession extracted from entry.\n";
    }
    
    print ">$accession\n$sequence\n";
}

exit(0);

