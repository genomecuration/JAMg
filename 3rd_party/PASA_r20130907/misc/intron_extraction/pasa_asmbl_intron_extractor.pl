#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use CdbTools;
use Nuc_translator;

my $usage = "usage: $0 alignments_described_file   genome_fasta_file\n";


my $alignments_described = $ARGV[0] or die $usage;
my $genome_fasta_file = $ARGV[1] or die $usage;


my $FLANK_LENGTH =10;
my %data;

my $intron_counter = 0;

open (my $fh, $alignments_described) or die $!;
while (<$fh>) {
    if (/^\#/) { next; }
    chomp;
    my ($genomic_acc, $subcluster_id, $pasa_acc, $cdna_list, $alignment_text) = split (/\t/);
    
    $alignment_text =~ m|orient\(a(\S)/s(\S)\)|;
    my $spliced_orient = $2;
    if ($spliced_orient eq '?') {
        next;
    }
    
    my @introns = &parse_introns($alignment_text);
    if (@introns) {
        $data{$genomic_acc}->{"$pasa_acc$;$spliced_orient"} = \@introns;
    }
}
close $fh;


foreach my $genomic_acc (keys %data) {
    ## get list of introns and attribute each to the list of pasa accs:
    
    my %intron_to_pasa_list;

    my $pasa_to_introns_href = $data{$genomic_acc};

    foreach my $pasa_key (keys %$pasa_to_introns_href) {
        my $introns_aref = $pasa_to_introns_href->{$pasa_key};

        my ($pasa_acc, $spliced_orient) = split (/$;/, $pasa_key);
        
        foreach my $intron_coordset (@$introns_aref) {
            
            my ($intron_lend, $intron_rend) = @$intron_coordset;

            my $intron_key = "${intron_lend}_${intron_rend}_${spliced_orient}";
           
            my $pasa_list_aref = $intron_to_pasa_list{$intron_key};
            unless (ref $pasa_list_aref) {
                $pasa_list_aref = $intron_to_pasa_list{$intron_key} = [];
            }

            push (@$pasa_list_aref, $pasa_acc);
            
        }
        
        
    }
    
    ## extract the introns:
    
    my $genome_seq = cdbyank_linear($genomic_acc, $genome_fasta_file);

    foreach my $pasa_key (keys %intron_to_pasa_list) {
        my @pasa_accs = @{$intron_to_pasa_list{$pasa_key}};

        $intron_counter++;

        my ($intron_lend, $intron_rend, $spliced_orient) = split (/_/, $pasa_key);

        ## wants 10 bp of nt on each end:
        
        $intron_lend -= $FLANK_LENGTH;
        $intron_rend += $FLANK_LENGTH;

        my $extract_length = $intron_rend - $intron_lend + 1;

        my $subseq = substr ($genome_seq, $intron_lend - 1, $extract_length);

        if ($spliced_orient eq '-') {
            $subseq = &reverse_complement($subseq);
        }
        $subseq = lc $subseq;
        my @chars = split (//, $subseq);
        for (my $i = 0; $i < $FLANK_LENGTH; $i++) {
            # make exon chars upcase
            $chars[$i] = uc $chars[$i];
            $chars[-1 * ($i + 1) ] = uc $chars[-1 * ($i + 1) ];
        }
        
        $subseq = join ("", @chars);
        
        my $pasa_listing = join (",", @pasa_accs);
        print "$genomic_acc\t$intron_lend-$intron_rend ($spliced_orient)\t$subseq\t$pasa_listing\n";
        
    }
    

    
}







exit(0);

####
sub parse_introns {
    my ($alignment_text) = @_;
    
    ## orient(a-/s-) align: 73652(685)-74094(243)>CT....AC<74286(242)-74527(1)

    my @introns;

    while ($alignment_text =~ m|\-(\d+)\(\d+\)\>\w{2}\.\.\.\.\w{2}\<(\d+)\(\d+\)|g) {
        my $exon_A_rend = $1;
        my $exon_B_lend = $2;
        
        my $intron_lend = $exon_A_rend + 1;
        my $intron_rend = $exon_B_lend - 1;

        push (@introns, [$intron_lend, $intron_rend]);

    }

    return (@introns);

}


    
    
