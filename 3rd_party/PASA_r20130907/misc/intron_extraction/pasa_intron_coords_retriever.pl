#!/usr/local/bin/perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});

use CdbTools;

my $usage = "usage: $0 alignments_described_file\n";


my $alignments_described = $ARGV[0] or die $usage;

my %data;

my $intron_counter = 0;

open (my $fh, $alignments_described) or die $!;
while (<$fh>) {
    chomp;
    my ($pasa_acc, $cdna_list, $genomic_acc, $alignment_text) = split (/\t/);
    
    $pasa_acc =~ s/alignment_assembly://;
    $genomic_acc =~ s/genomic_acc://;
    
    $alignment_text =~ m|orient\(a(\S)/s(\S)\)|;
    my $spliced_orient = $2;
    if ($spliced_orient eq '?') {
        next;
    }
    
	&parse_introns($genomic_acc, $alignment_text);
    
	print join("\n", keys %data) . "\n";
	
	
}


exit(0);

####
sub parse_introns {
    my ($genomic_acc, $alignment_text) = @_;
    
    ## orient(a-/s-) align: 73652(685)-74094(243)>CT....AC<74286(242)-74527(1)

    while ($alignment_text =~ m|\-(\d+)\(\d+\)\>\w{2}\.\.\.\.\w{2}\<(\d+)\(\d+\)|g) {
        my $exon_A_rend = $1;
        my $exon_B_lend = $2;
        
        my $intron_lend = $exon_A_rend + 1;
        my $intron_rend = $exon_B_lend - 1;
        
        

        my $intron_length = $intron_rend - $intron_lend + 1; 
        
        my $key = "$genomic_acc\t$intron_lend\t$intron_rend";
        $data{$key} = 1;
    }

	return;

}


    
    
