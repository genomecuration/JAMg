#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use CdbTools;
use Nuc_translator;

my $usage = "usage: $0 alignment_validations_file genome_fasta_file min_perID[=98]\n";

my $alignment_validations_file = $ARGV[0] or die $usage; 
my $genome_fasta_file = $ARGV[1] or die $usage;
my $min_perID = $ARGV[2] || 98;

my %data;

my $intron_counter = 0;

open (my $fh, $alignment_validations_file) or die $!;
<$fh>; # rid header
while (<$fh>) {
    chomp;
    unless (/\w/) { next; }
	my @x = split (/\t/);
	
	my $genomic_acc = $x[2];
	my $cdna_acc = $x[1];
	my $alignment_text = $x[10];
	    
	&parse_introns($genomic_acc, $alignment_text);
}

my %contig_to_intron_list;
foreach my $intron_key (keys %data) {
	my ($genomic_acc, $orient, $lend, $rend) = split (/\t/, $intron_key);
	
	push (@{$contig_to_intron_list{$genomic_acc}}, [$orient, $lend, $rend]);
}

foreach my $contig_acc (keys %contig_to_intron_list) {
	my @introns = @{$contig_to_intron_list{$contig_acc}};
	
	my $genome_seq = &cdbyank_linear($contig_acc, $genome_fasta_file);
	
	foreach my $intron (@introns) {
		my ($orient, $lend, $rend) = @$intron;
		
		my $intron_seq = substr($genome_seq, $lend - 1, $rend - $lend + 1);
		if ($orient eq '-') {
			$intron_seq = &reverse_complement($intron_seq);
		}
		
		my $intron_length = length($intron_seq);
		
		print "$contig_acc\t$orient\t$lend-$rend\t$intron_length\t$intron_seq\n";
	}
}


exit(0);

####
sub parse_introns {
    my ($genomic_acc, $alignment_text) = @_;
    
    ## orient(a-/s-) align: 73652(685)-74094(243)>CT....AC<74286(242)-74527(1)

    my @introns;

    while ($alignment_text =~ m|\-(\d+)\(\d+\)\>(\w{2})\.\.\.\.(\w{2})\<(\d+)\(\d+\)|g) {
		my $exon_A_rend = $1;
		my $splice_chars_A = $2;
		my $splice_chars_B = $3;
        my $exon_B_lend = $4;
    
		my $orient;
		if ($splice_chars_A eq 'GT' && $splice_chars_B eq 'AG') {
			$orient = '+';
		}
		elsif ($splice_chars_A eq 'CT' && $splice_chars_B eq 'AC') {
			$orient = '-';
		}
		else {
			# not a consensus splice site containing intron.
			next;
		}


        my $intron_lend = $exon_A_rend + 1;
        my $intron_rend = $exon_B_lend - 1;
        
        

        my $intron_length = $intron_rend - $intron_lend + 1; 
        
        my $key = "$genomic_acc\t$orient\t$intron_lend\t$intron_rend";
		$data{$key} = 1;
    }
	
	return;
}


    
    
