#!/usr/bin/env perl

use strict;
use warnings;
use List::Util qw (shuffle);
use Getopt::Long qw(:config no_ignore_case bundling);
use Cwd;
use FindBin;
use File::Basename;

my $param_string = <<__PARAMS__;

################# Evidence Modeler ###############################################################
#
#
# Opts for evidence modeler:
#
#  --genome            | -G   genome sequence in fasta format
#  --gene_predictions       | -g   gene predictions file
#  --protein_alignments | -p  protein alignments file
#  --est_alignments   | -e    est alignments file
#  --repeats          | -r    repeats masked from genome file
#  --terminalExons    | -t   supplementary file of additional terminal exons to consider
#  
#
####################################################################################################


__PARAMS__

    ;


my ($genomicSeqFile, $genePredictionsFile, $estAlignmentsFile, 
    $proteinAlignmentsFile, $repeatsFile, $terminalExonsFile);


&GetOptions ( 
              ## evm opts
              "genome|G=s"=>\$genomicSeqFile,
              "gene_predictions|g=s"=>\$genePredictionsFile,
              "est_alignments|e=s"=>\$estAlignmentsFile,
              "protein_alignments|p=s" => \$proteinAlignmentsFile,
              "repeats|r=s"=>\$repeatsFile,
              "terminalExonsFile|t=s"=>\$terminalExonsFile,
              );


unless ($genomicSeqFile) {
    die $param_string;
}


my $classified_file = "manually_classified_entries.txt";

my %seen;

my $viewer_cmd = $FindBin::Bin . "/../TkGFF3_viewer/TkGFF3_viewer.pl ";


if (-s $classified_file) {
    open (my $fh, $classified_file) or die $!;
    while (<$fh>) {
        my ($entry, $status) = split (/\s+/);
        $seen{$entry}=1;
    }
    close $fh;
}

open (my $outfh, ">>$classified_file") or die $!;

my @eles;
my $input_list = "all_train.entries";
open (my $fh, $input_list) or die "Error, cannot open file all_train.entries.  Make sure you're current directory is the root training directory and that the all_train.entries file exists" ;
while (<$fh>) {
    my ($entry, $rest) = split (/\s+/);
    $entry =~ s/\s+//g;
    if ($seen{$entry}) {
        next;
    }
    push @eles, $entry;
}

@eles = shuffle (@eles);

foreach my $entry (@eles) {
    
    unless (-s "$entry/template.gff3") {
        print STDERR "Sorry, cannot find $entry/template.gff3\n";
        next;
    }

    #my $cmd = "$viewer_cmd $entry/template.gff3 $entry/template.gff3 $entry/rice.genome $entry/gene_predictions.gff3 $entry/protein_alignments.gff3 $entry/transcript_alignments.gff3 $entry/pasa.terminal_exons.gff3 ";
    
    my $cmd = "$viewer_cmd --genome $entry/" . basename($genomicSeqFile);
    $cmd .= " --gene_predictions_listing $entry/template.gff3";
    if ($genePredictionsFile) {
        $cmd .= ",$entry/" . basename($genePredictionsFile);
    }
    if ($estAlignmentsFile || $repeatsFile || $terminalExonsFile) {
        $cmd .= " --alignments_listing ";
        $cmd .= " $entry/" . basename($estAlignmentsFile) . "," if $estAlignmentsFile;
        $cmd .= " $entry/" . basename($repeatsFile) . "," if $repeatsFile;
        $cmd .= " $entry/" . basename($terminalExonsFile) . "," if $terminalExonsFile;
        chop $cmd; #remove trailing comma
    }
    
    system $cmd;
    if ($?) {
        print STDERR "Error!  CMD: $cmd died with ret $?";
        die; next;
    }
    while (1) {
        print "Is $entry good?  [YN?]: ";
        my $answer = uc <STDIN>;
        chomp $answer;
        if ($answer =~/^[YN\?]$/) {
            print $outfh "$entry\t$answer\n";
            last; 
        }
    }
}

exit(0);

