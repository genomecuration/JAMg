#!/usr/local/bin/perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use Fasta_reader;


my $usage = "usage: $0 genome_fasta annot_coords lend-rend\n\n";

my $genome_file = $ARGV[0] or die $usage;
my $annot_coords = $ARGV[1] or die $usage;
my $seq_range = $ARGV[2] or die $usage;

my ($seq_lend, $seq_rend) = sort {$a<=>$b} split (/-/, $seq_range);
unless ($seq_lend && $seq_rend) {
    die "Error, coords range $seq_range is not as expected.\n\n$usage\n";
}

my $len = $seq_rend - $seq_lend + 1;


main: {


    ## create the subsequence

    my $fasta_reader = new Fasta_reader($genome_file);
    my $seq_obj = $fasta_reader->next();
    
    my $sequence = $seq_obj->get_sequence();
    my $subseq = substr($sequence, $seq_lend - 1, $len);

    my $accession = $seq_obj->get_accession();

    &write_sequence_file("$accession.$seq_range.seq", $subseq);

    
    ## get the annotation coordinates and adjust them.
    &parse_annotation_coordinates("$accession.$seq_range.coords", $annot_coords, $seq_lend, $seq_rend);
    

    exit(0);

}

####
sub write_sequence_file {
    my ($accession, $sequence) = @_;

    $sequence =~ s/(\S{60})/$1\n/g;

    open (my $fh, ">$accession") or die "Error, cannot write file $accession";
    print $fh ">$accession\n$sequence\n";
    close $fh;
    
    return;
}
    

####
sub parse_annotation_coordinates {
    my ($filename, $coords_file, $range_lend, $range_rend) = @_;
    
    my $text = `cat $coords_file`;
    my @genes = split (/\n\n/, $text);

    
    my $adj_text = "";


    foreach my $gene (@genes) {
        my @x = split (/\n/, $gene);
        my $header = shift @x;
        my ($orient, $lend, $rend, $annot) = split (/\s+/, $header);
        if ($lend >= $range_lend && $rend < $range_rend) {
            ## dump gene w/ adjusted coordinates:
            $lend = $lend - $range_lend + 1;
            $rend = $rend - $range_lend + 1;
            $adj_text .= "$orient $lend $rend $annot\n";
            foreach my $exon (@x) {
                my ($lend, $rend, $exon_text) = split (/\s+/, $exon);
                $lend = $lend - $range_lend + 1;
                $rend = $rend - $range_lend + 1;
                $adj_text .= "$lend $rend $exon_text\n";
            }
            $adj_text .= "\n"; # add separator between gene entries.
        }
    }
    
    open (my $ofh, ">$filename") or die "Error, cannot write file $filename";
    print $ofh $adj_text;
    close $ofh;


    return;
}
