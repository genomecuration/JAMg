#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");

use Fasta_reader;
use GTF_to_geneobjs;
use GFF3_to_geneobjs;
use Accuracy_Eval;
use Getopt::Long qw(:config no_ignore_case bundling);
use Data::Dumper;

umask(0000);

my $usage = <<_EOUSAGE_;

####################################################################################################
#  
#  --genome                  genome sequence in multi-fasta format
#  --format                  GTF or GFF3
#  --reference_genes         reference genes file
#  --predicted_genes         predicted genes file
#  -v                        verbose
#####################################################################################################

_EOUSAGE_

    ;

my ($reference_genes_file, $predicted_genes_file, $genome_file, $format, $help);

our $SEE;

&GetOptions ( 
              "reference_genes=s" => \$reference_genes_file,
              "predicted_genes=s" => \$predicted_genes_file,
              "genome=s" => \$genome_file,
              "help|h" => \$help,
              "format=s" => \$format,
              "v" => \$SEE,
              );

if ($help) { die $usage; }
unless ($reference_genes_file && $predicted_genes_file && $genome_file && $format =~ /^(GTF|GFF3)$/) { 
    die $usage;
}

main: {
        
    ## get the genes:
    my %ref_gene_id_to_gene;
    my $ref_seqname_map_href;
    my %other_gene_id_to_gene;
    my $other_seqname_map_href;
    
    if ($format eq 'GTF') {
        print STDERR "-processing GTF files\n" if $SEE;
        $ref_seqname_map_href = &GTF_to_geneobjs::parse_file($reference_genes_file, \%ref_gene_id_to_gene);
        $other_seqname_map_href = &GTF_to_geneobjs::parse_file($predicted_genes_file, \%other_gene_id_to_gene);
    }
    else {
        print STDERR "-processing GFF3 files\n" if $SEE;
        $ref_seqname_map_href = &GFF3_to_geneobjs::parse_file($reference_genes_file, \%ref_gene_id_to_gene);
        $other_seqname_map_href = &GFF3_to_geneobjs::parse_file($predicted_genes_file, \%other_gene_id_to_gene);
    }
    
    # get the sequence lengths:
    my $fasta_reader = new Fasta_reader($genome_file);
    my %contig_lengths;
    while (my $seq_obj = $fasta_reader->next()) {
        my $accession = $seq_obj->get_accession();
        my $sequence = $seq_obj->get_sequence();
        my $seq_length = length($sequence);
        $contig_lengths{$accession} = $seq_length;
        print "$accession length: $seq_length\n" if $SEE;
    }
    

    ## perform SnSp analysis for each contig.
    my $accuracy_obj = new Accuracy_Eval();
    
    foreach my $contig (keys %contig_lengths) {
        
        print "// processing $contig\n" if $SEE;
        
        my $seq_length = $contig_lengths{$contig} or die "Error, cannot find length for contig $contig from sequence in $genome_file";
        
        my $ref_gene_ids_aref = $ref_seqname_map_href->{$contig};
        my @ref_genes;
        if (ref $ref_gene_ids_aref) {
            foreach my $ref_gene_id (@$ref_gene_ids_aref) {
                my $gene_obj = $ref_gene_id_to_gene{$ref_gene_id} or die "Error, no gene retrieved from reference via $ref_gene_id";
                push (@ref_genes, $gene_obj);
            }
        }

        my @other_genes;
        my $other_gene_ids_aref = $other_seqname_map_href->{$contig};
        if (ref $other_gene_ids_aref) {
            foreach my $other_gene_id (@$other_gene_ids_aref) {
                my $other_gene = $other_gene_id_to_gene{$other_gene_id} or die "Error, no gene retreived from others via $other_gene_id";
                push (@other_genes, $other_gene);
            }
        }

        $accuracy_obj->add_entry(1, $seq_length, \@ref_genes, \@other_genes);

    }


    ## generate report.
    my ($Nuc_Sn, $Nuc_Sp) = $accuracy_obj->compute_nucleotide_SnSp();
    my ($Exon_Sn, $Exon_Sp) = $accuracy_obj->compute_exon_SnSp();
    my ($Transcript_Sn, $Transcript_Sp) = $accuracy_obj->compute_transcript_SnSp();
    my ($Gene_Sn, $Gene_Sp) = $accuracy_obj->compute_gene_SnSp();

    
    print "\nReference: $reference_genes_file\tPredictions: $predicted_genes_file\n"
        . "--------------------------------------------------------------------\n";  
    print $accuracy_obj->get_summary_stats();

    # print Dumper ($accuracy_obj);
    
    exit(0);
}

