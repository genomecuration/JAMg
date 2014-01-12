#!/usr/local/bin/perl

use strict;
use warnings;
use lib ($ENV{EGC_SCRIPTS});
use Egc_library;


my $SEE = 0;

my %data;

while (<STDIN>) {
    chomp;
    my ($pasa_acc, $cdna_list, $genome_acc, $alignment) = split (/\t/);
    
    $pasa_acc =~ s/alignment_assembly://;
    $genome_acc =~ s/genomic_acc://;

    $alignment =~ /orient\(a[+-]\/s([+-\?])/;
    my $spliced_orient = $1;

    if ($spliced_orient eq '?') {
        next;
    }

    print "$pasa_acc\t$spliced_orient\n" if $SEE;
    
    while (m|(\d+)\(\d+\)>(\w{2})\.\.\.\.(\w{2})<(\d+)|g) {
        my ($c1, $s1, $s2, $c2) = ($1, $2, $3, $4);

        print "$c1 $s1 $s2 $c2\n" if $SEE;
    
        if ($spliced_orient eq '+') {
            
            my $donor_lend = $c1 + 1;
            my $donor_rend = $c1 + 2;

            my $donor_entry = "DON: $genome_acc-$donor_lend-$donor_rend-[$spliced_orient]-$s1";
            
            $data{$donor_entry} = 1;
            

            my $acceptor_lend = $c2 - 2;
            my $acceptor_rend = $c2 - 1;
            
            my $acceptor_entry = "ACC: $genome_acc-$acceptor_lend-$acceptor_rend-[$spliced_orient]-$s2";
            
            $data{$acceptor_entry} = 1;
            
        }
        
        
        else {
            # minus strand
            
            my $acceptor_lend = $c1 + 1;
            my $acceptor_rend = $c1 + 2;
            $s1 = &reverse_complement($s1);

            my $acceptor_entry = "ACC: $genome_acc-$acceptor_lend-$acceptor_rend-[$spliced_orient]-$s1";
            
            $data{$acceptor_entry} = 1;
            

            my $donor_lend = $c2 - 2;
            my $donor_rend = $c2 - 1;
            $s2 = &reverse_complement($s2);
            
            my $donor_entry = "DON: $genome_acc-$donor_lend-$donor_rend-[$spliced_orient]-$s2";
            
            $data{$donor_entry} = 1;
        }

    }
    
}

foreach my $splice_site (keys %data) {
    print "$splice_site\n";
}



exit(0);
    
