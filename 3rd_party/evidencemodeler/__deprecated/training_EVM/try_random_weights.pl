#!/usr/bin/env perl

use warnings;
use strict;
use FindBin;
use Getopt::Long qw(:config no_ignore_case bundling);
use Cwd;
use Data::Dumper;
use Carp;
use lib ("$FindBin::Bin/../PerlLib");
use File::Basename;

my $bindir = $FindBin::Bin;

my $param_string = <<__PARAMS__;

################# Evidence Modeler ##############################
#
#
# Opts for evidence modeler:
#
#  --genome                 genome sequence in fasta format
#  --gene_predictions       gene predictions file
#  --protein_alignments     protein alignments file
#  --transcript_alignments     transcript alignments file
#  --repeats                   repeats masked from genome file
#  --terminalExons             supplementary file of additional terminal exons to consider
#  --stitch_ends             alignment types (prognames) to apply end stitching into existing exons (ie. 'genewise,alignAssembly')
#  --extend_to_terminal      extend terminal alignment segment into a terminal exon (ie. 'genewise')
#  --stop_codons             list of stop codons (default: TAA,TGA,TAG)
#                               *for Tetrahymena, set --stop_codons TGA
# Opts for training software:
#
#  required:
#  --weights_file  | -w         file containing weights to initialize to. 
#  --training_entries            file containing the list of entries to be used for training purposes.
#
#  optional:  
#  --num_rand_iterations      number of fully random weight settings to explore before doing gradient descent on the best scoring combo. (default 20)
#
#
# flags:
#  -S                    verbose flag
#  --debug               debug mode, writes lots of extra files.
#
#  --use_fixed_genefinder_weights  given initialized weights in weights file, or computed best relative weights, keep genefinder 
#                                       weights relatively fixed when computing other best weights
#  --just_score_existing_weights   file containing ev_type\\tweight, and scoring is performed.
#  --maximize_ev_type_weight       ev_type to find the maximum performing weight value given others in a weights file.
#                                  **formatted like so  ev_class,ev_type   ie.   TRANSCRIPT,alignAssembly-pasaRunName           
#
#  --just_measure_other_prediction_accuracies    flag to report Sn,Sp values for the gene predictions used as evidence
#
#  --use_grid_computing            ** tigr-only  (additional opt: --cmds_per_node, default's num_entries/100 || 1 )
#
#  --use_trained_gene_prediction_weights     file containing already trained genefinder weights in the weights file format. (resulting from phase 1 of training)  
#           ** assumes --use_fixed_genefinder_weights   
#  --use_best_individual_weights             file containing the best individual evidence weights from phase 2 of training process.
#           ** requires above --use_trained_gene_prediction_weights
# 
#  --NO_RECURSE                    no recursion into unpredicted regions or long introns.
#
#################################################################################################################################



__PARAMS__

    ;

$|++;

my ($genomicSeqFile, 
    $genePredictionsFile, 
    $estAlignmentsFile, 
    $proteinAlignmentsFile,
    $repeatsFile, 
    $terminalExonsFile, 
    $SEE, 
    $DEBUG, 
    $stop_codons_arg, 
    $help_flag, 
    $stitch_ends,
    $extend_to_terminal,
    $weights_file, 
    $train_listing_file, 
    $use_grid_computing_flag,
    
    $NO_RECURSE,
        
    );


my $num_rand_iterations = 20; #default

my $NO_RECURSE_FLAG = 0;
my $cmds_per_node;

&GetOptions ( 
              ## evm opts
              "genome=s"=>\$genomicSeqFile,
              "gene_predictions=s"=>\$genePredictionsFile,
              "transcript_alignments=s"=>\$estAlignmentsFile,
              "protein_alignments=s" => \$proteinAlignmentsFile,
              "repeats=s"=>\$repeatsFile,
              "terminalExonsFile=s"=>\$terminalExonsFile,
              "S"=>\$SEE,
              "debug" => \$DEBUG,
              "stop_codons=s" => \$stop_codons_arg,
              "help|h" => \$help_flag,
              "stitch_ends=s" => \$stitch_ends,
              "extend_to_terminal=s" => \$extend_to_terminal,
            
              
              ## training ops
              "num_rand_iterations=i" => \$num_rand_iterations,
              "weights_file|w=s" => \$weights_file,
              "training_entries=s" => \$train_listing_file, 
              
              ## misc
              "use_grid_computing" => \$use_grid_computing_flag,
              
              "NO_RECURSE" => \$NO_RECURSE_FLAG,
              'cmds_per_node=i' => \$cmds_per_node,
              

              );

if ($help_flag) {
    die "$param_string\n\n";
}


unless ($genomicSeqFile && $genePredictionsFile && $train_listing_file && $weights_file) {
    die $param_string;
}


main: {

    srand();
    umask(0000);
    
    ## build the command string:
    my $EVM_cmd_template = "$bindir/train_weights.pl"
        . " --genome " . basename($genomicSeqFile)
        . " --gene_predictions " . basename($genePredictionsFile) 
        . " --training_entries $train_listing_file ";
        
    if ($repeatsFile) {
        $EVM_cmd_template .= " --repeats " . basename($repeatsFile);
    }
    
    if ($stop_codons_arg) {
        $EVM_cmd_template .= " --stop_codons $stop_codons_arg";
    }
    
    if ($NO_RECURSE_FLAG) {
        $EVM_cmd_template .= " --NO_RECURSE";
    }
    
    my $JUST_GENEFINDERS_cmd_template = $EVM_cmd_template;
    
    if ($estAlignmentsFile) {
        $EVM_cmd_template .= " --transcript_alignments " . basename($estAlignmentsFile);
    }
    if ($proteinAlignmentsFile) {
        $EVM_cmd_template .= " --protein_alignments " . basename($proteinAlignmentsFile);
    }
    
    if ($stitch_ends) {
        $EVM_cmd_template .= " --stitch_ends " . basename($stitch_ends);
    }

    if ($extend_to_terminal) {
        $EVM_cmd_template .= " --extend_to_terminal " . basename($extend_to_terminal);
    }
    
    if ($terminalExonsFile) {
        $EVM_cmd_template .= " --terminalExons " . basename($terminalExonsFile);
    }

    if ($use_grid_computing_flag) {
        $EVM_cmd_template .= " --use_grid_computing ";
                
        if ($cmds_per_node) {
            $EVM_cmd_template .= " --cmds_per_node $cmds_per_node ";
        }
    }
    
    $EVM_cmd_template .= " --just_score_existing_weights ";


    &evaluate_random_weights($EVM_cmd_template, $weights_file);


    exit(0);
}


####
sub evaluate_random_weights {
    my ($EVM_cmd_template, $weights_file) = @_;

    # parse weights file
    my @evidences;
    open (my $fh, "$weights_file") or die "Error, cannot open weights file: $weights_file";
    while (<$fh>) {
        chomp;
        unless (/\w/) { next; }
        if (/^\#/) { next; }
        
        my ($ev_class, $ev_type, $weight) = split (/\s+/);
        push (@evidences, { ev_type => $ev_type,
                            ev_class => $ev_class,
                            weight => undef, } );
    }
    close $fh;

    for (1 .. $num_rand_iterations) {
        
        ## set random weight values
        my $sum_random = 0;
        foreach my $evidence (@evidences) {
            my $weight = rand(1);
            $evidence->{weight} = $weight;
            $sum_random += $weight;
        }
        # normalize to one
        foreach my $evidence (@evidences) {
            $evidence->{weight} /= $sum_random;
        }

        # create a weights file:
        my $tmp_weight_file = "weights.tmp.$$";
        open (my $fh, ">$tmp_weight_file") or die $!;
        foreach my $evidence (@evidences) {
            my ($ev_class, $ev_type, $weight) = ($evidence->{ev_class}, $evidence->{ev_type}, $evidence->{weight});
            print $fh "$ev_class\t$ev_type\t$weight\n";
        }
        close $fh;

        my $cmd = "$EVM_cmd_template --weights_file $tmp_weight_file";
        
        my $result = `$cmd`;
        
        if ($?) {
            die "$result\n\nError, cmd: $cmd died with ret $? ";
        }

        print $result;

    }
    

    return;
}


    
    
    
    
