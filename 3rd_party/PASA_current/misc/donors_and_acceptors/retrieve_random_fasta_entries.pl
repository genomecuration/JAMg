#!/usr/local/bin/perl

use strict;
use warnings;
use List::Util qw (shuffle);
use lib ($ENV{EUK_MODULES});
use Fasta_reader;

my $usage = "usage: $0 fasta_file num_rand_entries\n";

my $fasta_file = $ARGV[0] or die $usage;
my $num_rand_entries = $ARGV[1] or die $usage;

srand;

unless ($num_rand_entries >= 1) {
    die $usage;
}

main: {
    my $num_fasta_entries = `grep '>' $fasta_file | wc -l `;
    $num_fasta_entries =~ s/\s//g;
    
    my %rand_indices;
    my $count = 0;
    while ($count <= $num_rand_entries) {
        my $index = int (rand($num_fasta_entries) );
        $index++; #start at one, not zero
        unless ($rand_indices{$index}) {
            $count++;
            $rand_indices{$index} = 1;
        }
    }

    my $fasta_reader = new Fasta_reader($fasta_file);
    
    ## stream the fasta file and print the randomly selected entries:
    $count = 0;
    while (my $seq_obj = $fasta_reader->next()) {
        $count++;
        if ($rand_indices{$count}) {
            print $seq_obj->get_FASTA_format();
            delete $rand_indices{$count}; # entry accounted for.
        }
        
        unless (%rand_indices) {
            last; # don't keep reading the file if nothing left to print.
        }
    }
    
    if (%rand_indices) {
        die "Error, not all random entries were accounted for in the output.\n";
    }
}

exit(0);
        
