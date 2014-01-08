#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;


my $usage = "usage: $0 < file_listing \n\n";

main: {
    my @files = <STDIN>;
    chomp @files;

    my @structs;

    foreach my $file (@files) {
        my $struct = &parse_myEval($file);
        
        # print Dumper ($struct);

        push (@structs, $struct);
    }

    ## report
    @structs = sort {$a->{file} cmp $b->{file}} @structs;

    foreach my $struct (@structs) {
        my ($gene_count, $tpg, 
            $nuc_Sn, $nuc_Sp, 
            $exon_Sn, $exon_Sp, 
            $transcript_Sn, $transcript_Sp, 
            $gene_Sn, $gene_Sp) = ($struct->{'gene count'}, $struct->{TpG},
                                   $struct->{'Nucleotide Sensitivity'}, $struct->{'Nucleotide Specificity'},
                                   $struct->{'Exon Sensitivity'}, $struct->{'Exon Specificity'},
                                   $struct->{'Transcript Sensitivity'}, $struct->{'Transcript Specificity'},
                                   $struct->{'Gene Sensitivity'}, $struct->{'Gene Specificity'});
        
        my $filename = $struct->{file};
        
        print "$filename\t$gene_count\t$tpg\t$nuc_Sn\t$nuc_Sp\t$exon_Sn\t$exon_Sp\t$transcript_Sn\t$transcript_Sp\t$gene_Sn\t$gene_Sp\n";
    }
    


    exit(0);
}

    
####
sub parse_myEval {
    my ($file) = @_;

    my %data;

    open (my $fh, $file) or die "Error, cannot open $file";
    while (<$fh>) {
        chomp;
        if (/^(Gene|Transcript|Exon|Nucleotide)/) {
            my @x = split (/\s+/);
            my $percent = pop @x;
            $percent =~ s/\%//;
            my $token = join (" ", @x);
            $data{$token} = $percent;
        }
        elsif (/Predicted genes: (\d+)/) {
            $data{'gene count'} = $1;
        }
        elsif (/Predicted transcripts: \d+ \(([^\)]+)\)/) {
            $data{'TpG'} = $1;
        }
    }

    $data{file} = $file;
    
    return (\%data);
}


