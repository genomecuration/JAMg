#!/usr/bin/env perl

use warnings;
no warnings 'once';
use strict;
use FindBin;
use Getopt::Long qw(:config no_ignore_case bundling);
use Cwd;
use Data::Dumper;
use Carp;
use lib ("$FindBin::Bin/../PerlLib");
use GFF3_utils;
use Gene_obj_indexer;
use SnSp::SnSp_analysis_manager;

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
    $JUST_SCORE_GIVEN_WEIGHTS, 
    $USE_FIXED_GENEFINDER_WEIGHTS, 
    $weights_file, 
    $JUST_MEASURE_OTHER_PREDICTION_ACCURACIES, 
    $train_listing_file, 
    $maximize_ev_type_weight, 
    $use_grid_computing_flag,

    $use_trained_gene_prediction_weights,
    $use_best_individual_weights,

    $NO_RECURSE,
    $cmds_per_node,
    
    );


my $num_rand_iterations = 20; #default

my $gene_prediction_gradient_lambda = 0.2; # +- of existing value for walking weights along the gradient descent 
my $gene_prediction_max_lambda = 10; ## cap the weight adjustments for gradient descent.

my $other_ev_max_lambda = 100;
my $other_ev_gradient_lambda = 0.5;

my $LAMBDA_AMPLIFICATION_FACTOR = 1.5; # increase lambda value by this factor. Used to set at two, but doubling was too aggresive

my $required_gradient_improvement_points = 1;

my $FLANK_REGIONS_SnSp_BP = 500; # default for regions flanking the reference gene structure to be used in SnSp calculations.

my $NO_RECURSE_FLAG = 0;

my $INTERGENIC_SCORE_ADJUST_FACTOR = 1;


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
              "just_score_existing_weights" => \$JUST_SCORE_GIVEN_WEIGHTS,
              "weights_file|w=s" => \$weights_file,
              "use_fixed_genefinder_weights" => \$USE_FIXED_GENEFINDER_WEIGHTS,
              "just_measure_other_prediction_accuracies" => \$JUST_MEASURE_OTHER_PREDICTION_ACCURACIES,
              "training_entries=s" => \$train_listing_file, 
              "maximize_ev_type_weight=s" => \$maximize_ev_type_weight,
              
              ## misc
              "use_grid_computing" => \$use_grid_computing_flag,
              
              "use_trained_gene_prediction_weights=s" => \$use_trained_gene_prediction_weights,
              "use_best_individual_weights=s" => \$use_best_individual_weights,

              "NO_RECURSE" => \$NO_RECURSE_FLAG,
              'cmds_per_node=i' => \$cmds_per_node,
              

              );

if ($help_flag) {
    die "$param_string\n\n";
}


unless ($genomicSeqFile && $genePredictionsFile && $train_listing_file && $weights_file) {
    die $param_string;
}

if ($use_trained_gene_prediction_weights) {
    $USE_FIXED_GENEFINDER_WEIGHTS = 1;
}


my $MIN_EUCLIDIAN_DISTANCE = 0.1; # minimum euclidian distance between closest random weight set previously employed for valid new random weight set.
my $MAX_RAND_ITERATION_TRY_COUNT = 100; #max no of tries to find combination with distance above away from all other already tried combos.

my $JUST_SCORE_GIVEN_WEIGHTS_FLAG = 0;
if ($JUST_SCORE_GIVEN_WEIGHTS) {
    $JUST_SCORE_GIVEN_WEIGHTS_FLAG = 1;
}

my %EVIDENCE_CLASS;  # ev_type => class, where class is in PROTEIN,TRANSCRIPT,ABINITIO_PREDICTION, OTHER_PREDICTION

my $READY_FOR_GRID_SUBMISSION = 0; 

my $workdir = cwd();

main: {

    srand();
    umask(0000);
    
    my $evidence_modeler_prog = "$bindir/../evidence_modeler.pl";
    
    my @train_entries = &get_train_listing();
    
    my %weights = ();
    &populate_weights_from_file($weights_file, \%weights);
    
    ## build the command string:
    my $EVM_cmd_template = "$evidence_modeler_prog"
        . " --genome __PATH__/$genomicSeqFile"
        . " --gene_predictions __PATH__/$genePredictionsFile";
    
    if ($repeatsFile) {
        $EVM_cmd_template .= " --repeats __PATH__/$repeatsFile";
    }
    
    if ($stop_codons_arg) {
        $EVM_cmd_template .= " --stop_codons $stop_codons_arg";
    }


    if ($NO_RECURSE_FLAG) {
        $EVM_cmd_template .= " --NO_RECURSE";
    }
    
    my $JUST_GENEFINDERS_cmd_template = $EVM_cmd_template;
    
    if ($estAlignmentsFile) {
        $EVM_cmd_template .= " --transcript_alignments __PATH__/$estAlignmentsFile";
    }
    if ($proteinAlignmentsFile) {
        $EVM_cmd_template .= " --protein_alignments __PATH__/$proteinAlignmentsFile";
    }
    
    if ($stitch_ends) {
        $EVM_cmd_template .= " --stitch_ends $stitch_ends";
    }

    if ($extend_to_terminal) {
        $EVM_cmd_template .= " --extend_to_terminal $extend_to_terminal";
    }
    
    if ($terminalExonsFile) {
        $EVM_cmd_template .= " --terminalExons __PATH__/$terminalExonsFile";
    }
    
    my @prediction_prog_list;
    my @non_predictions_list;
    my @other_predictions;
    foreach my $ev_type (keys %EVIDENCE_CLASS) {
        if ( (my $class = $EVIDENCE_CLASS{$ev_type}) eq "ABINITIO_PREDICTION") {
            push (@prediction_prog_list, $ev_type);
        }
        else {
            push (@non_predictions_list, $ev_type);
            if ($class eq 'OTHER_PREDICTION') {
                push (@other_predictions, $ev_type);
            }
        }
    }
    
    if ($JUST_MEASURE_OTHER_PREDICTION_ACCURACIES) {
        &measure_other_prediction_accuracies(\@train_entries, [@prediction_prog_list, @other_predictions]);
        exit(0);
    }
    
    
    ## First, train weights for genefinders alone:
    if ($JUST_SCORE_GIVEN_WEIGHTS_FLAG) {

                
        my ($BJHscore, $percent_gene_accuracy, $percent_exon_accuracy) = &score_weight_combination(\@train_entries, $EVM_cmd_template, \%weights);
        exit(0);
        
    }
    
    if ($maximize_ev_type_weight) {
        my ($max_ev_class, $max_ev_type) = split (/,/, $maximize_ev_type_weight);
        unless ($max_ev_class && $max_ev_type) {
            die "Error, cannot extract evidence class and type from maximize_ev_type_weight parameter: $maximize_ev_type_weight";
        }
        $EVIDENCE_CLASS{$max_ev_type} = $max_ev_class;
        
        %weights = &find_maximum_weight_for_ev_type(\@train_entries, $EVM_cmd_template, \%weights, $max_ev_type);
    }
    else {
        
        ## RANDOMIZATION and/or GRADIENT DESCENT
        my %fixed_weights;    
        if ($USE_FIXED_GENEFINDER_WEIGHTS) {
            
            if (my $genefinder_weights_file = $use_trained_gene_prediction_weights) {
                
                ## Read apriori geneprediction weights.
                open (my $fh, $genefinder_weights_file) or die "Error, cannot open $genefinder_weights_file";
                
                my %genefinder_weights;
                my $weight_sum = 0;
                while (<$fh>) {
                    unless (/\w/) { next; }
                    if (/^\#/) { next; } # commented line
                    chomp;
                    my ($ev_class, $ev_type, $weight_val) = split (/\s+/);
                    if ($ev_class eq "ABINITIO_PREDICTION") {
                        $genefinder_weights{$ev_type} = $weight_val;
                        $weight_sum += $weight_val;
                    }
                }
                close $fh;
                
                ## normalize so sum to one.
                foreach my $genefinder (keys %genefinder_weights) {
                    my $normalized_weight_value = $genefinder_weights{$genefinder} / $weight_sum;
                    $fixed_weights{$genefinder} = $normalized_weight_value;
                }
                
                print "fixed weights normalized to : " . Dumper (\%fixed_weights);
                %weights = %fixed_weights;
            }
            
            else {
                                
                ## train them separately:
    
                ## find an appropriate intergenic region scaling factor.
                ## set each prediction weight equally, then find the scaling factor that optimizes accuracy
                
                # &compute_intergenic_score_adjust_factor (\@train_entries, $EVM_cmd_template, \@prediction_prog_list);
                
                print "\n\n-----------------------------------------------------\n"
                  . "MODE: RANDOM_THEN_GRADIENT\n\n";
                print "*-training fixed weights separately.\n";
                %fixed_weights = &train_weights(\@train_entries, 
                                                $JUST_GENEFINDERS_cmd_template, 
                                                \@prediction_prog_list, 
                                                {}, #start w/ no predefined weights.
                                                { gradient_lambda => $gene_prediction_gradient_lambda,
                                                  max_lambda => $gene_prediction_max_lambda,
                                              }
                                                );
                print "-relative genefinder fixed weights fixed to: " . Dumper (\%fixed_weights);
                
                ## store them in a separate file:
                &write_weights_file("trained_genefinder_weights.$$.txt", \%fixed_weights);
                
                %weights = %fixed_weights;
                

            }
        }
                
        my @ev_types; #all the ev_types that do not have fixed weights.
        if (%fixed_weights) {
            foreach my $ev_type (keys %EVIDENCE_CLASS) {
                if (! exists $fixed_weights{$ev_type}) {
                    push (@ev_types, $ev_type);
                }
            }
        }
        else {
            @ev_types = keys %EVIDENCE_CLASS;
        }
        
        if (@ev_types) {
                        
            my %best_individual_weights;
        
            if ($use_best_individual_weights) {
                ## parse values from existing data:
                open (my $fh, $use_best_individual_weights) or die "Error, cannot open file $use_best_individual_weights";
                while (<$fh>) {
                    chomp;
                    unless (/\w/) { next; }
                    my ($ev_class, $ev_type, $weight) = split (/\s+/);
                    unless ($weight =~ /\d/) { die "Error, no numeric value in weight for entry: $_";}
                    $best_individual_weights{$ev_type} = $weight;
                }
                close $fh;
            }
            else {
                ## compute them:
                
                %best_individual_weights = &optimize_individual_evidence_weights(\@train_entries, 
                                                                                 $EVM_cmd_template, 
                                                                                 \@ev_types, 
                                                                                 \%fixed_weights,
                                                                                 { gradient_lambda => $other_ev_gradient_lambda,
                                                                                   max_lambda => $other_ev_max_lambda
                                                                                   }
                                                                                 );
            }
            
            ## Now do regular grad descent with all of them:
            print "-----------------------------------------------------------------\n"
                . "MODE: GRADIENT_ALL_EVIDENCE [Round 1, heavy lambda]\n\n";
            
            %weights = &perform_gradient_descent(\@train_entries, $EVM_cmd_template, \%best_individual_weights,
                                                 \%fixed_weights, { gradient_lambda => $other_ev_gradient_lambda,
                                                                    max_lambda => $other_ev_max_lambda });
            
            
            
            
            
            my %var_weights;
            foreach my $ev_type (keys %weights) {
                unless (exists $fixed_weights{$ev_type}) {
                    $var_weights{$ev_type} = $weights{$ev_type};
                }
            }
            
            print "-----------------------------------------------------------------\n"
                . "MODE: GRADIENT_ALL_EVIDENCE [Round 2, light lambda]\n\n";
            ## perform gradient descent using a lighter step (as w/ gene predictions).
            %weights = &perform_gradient_descent(\@train_entries,
                                                 $EVM_cmd_template,
                                                 \%var_weights,
                                                 \%fixed_weights,
                                                 { gradient_lambda => $gene_prediction_gradient_lambda,
                                                   max_lambda => $gene_prediction_max_lambda,
                                               } 
                                                 );
            
        }
    }
        
    my $final_weights_txt = "";
    
    ## print final weights:
    print "Trained weights are as follows:\n";
    foreach my $ev_type (reverse sort {$weights{$a}<=>$weights{$b}} keys %weights) {
        my $ev_class = $EVIDENCE_CLASS{$ev_type};
        $final_weights_txt .= "$ev_class\t$ev_type\t" . $weights{$ev_type} . "\n";
    }
    
    open (my $final_weights_fh, ">final_optimized_weights.$$.txt") or die $!;
    print $final_weights_fh $final_weights_txt;
    close $final_weights_fh;
    
    print "$final_weights_txt\n";
    
    exit(0);
}


####
sub optimize_individual_evidence_weights {
    my ($train_entries_aref, $EVM_cmd_template, $ev_types_aref, $fixed_weights_href, $lambdas_href) = @_;

    my $average_fixed_weight = 1;
    my @fixed_weight_ev_types = keys %$fixed_weights_href;
    my $num_fixed_ev_types = scalar @fixed_weight_ev_types;
    if ($num_fixed_ev_types > 0) {
        my $sum_weights = 0;
        foreach my $val (values %$fixed_weights_href) {
            $sum_weights += $val;
        }
        $average_fixed_weight = $sum_weights / $num_fixed_ev_types;
    }

    my %best_individual_weights;

    ## find the lowest weight value that provides the optimal score:
    foreach my $ev_type (@$ev_types_aref) {
        
        print "-------------------------------------------------------------\n"
          . "MODE: INTRODUCE_$ev_type\n\n";
        
        my @weights_and_scores;

        # first score by introducing the evidence with zero weight:

        my ($BJHscore, $percent_gene_accuracy, $percent_exon_accuracy) = &score_weight_combination($train_entries_aref, $EVM_cmd_template, { %$fixed_weights_href, $ev_type => 0 } );
        
        push (@weights_and_scores, { score => $BJHscore,
                                     weights => { %$fixed_weights_href, $ev_type => 0 } });  # include zero weighting in score distribution.
        
        # now set to average weight:
        
        my %grad_desc_weights = &perform_gradient_descent($train_entries_aref, $EVM_cmd_template, 
                                                          { $ev_type => $average_fixed_weight} , 
                                                          $fixed_weights_href, $lambdas_href, \@weights_and_scores);
        
        @weights_and_scores = sort { $b->{score}<=>$a->{score} # desc by score 
                                     ||
                                    $a->{weights}->{$ev_type} <=> $b->{weights}->{$ev_type} # asc by weight
                                 }  @weights_and_scores;

        my $best_weight_score_combo = shift @weights_and_scores;
        
        print "Best weight and score combo is: " . Dumper ($best_weight_score_combo);

        $best_individual_weights{$ev_type} = $best_weight_score_combo->{weights}->{$ev_type};
        
    }


    ## archive them:
    {
        open (my $fh, ">best_individual_weights.$$.txt") or die $!;
        foreach my $ev_type (reverse sort {$best_individual_weights{$a}<=>$best_individual_weights{$b}} keys %best_individual_weights) {
            my $weight = $best_individual_weights{$ev_type};
            my $ev_class = $EVIDENCE_CLASS{$ev_type};
            print $fh "$ev_class\t$ev_type\t$weight\n";
        }
        close $fh;
    }
    
    
    return (%best_individual_weights);
}



####
sub get_train_listing {
    
    unless (-s $train_listing_file) {
        die "Error, cannot locate file $train_listing_file ";
    }
    
    my @entries;
    open (my $fh, $train_listing_file) or die "Error, cannot open $train_listing_file";
    while (<$fh>) {
        s/\s+//g;
        push (@entries, $_);
    }
    close $fh;
    return (@entries);
}

####
sub parse_ev_types {
    my ($gff3_filename, $ev_types_href) = @_;
    print "init: Parsing ev_types from $gff3_filename\n";
    open (my $fh, $gff3_filename) or die "Error, cannot open $gff3_filename";
    while (<$fh>) {
        if (/^\#/) { next;}
        unless (/\w/) { next;}
        my @x = split (/\t/);
        my $ev_type = $x[1];
        $ev_types_href->{$ev_type} = 1;
    }
    return;
}

####
sub ensure_disjoint {
    my @entries = @_;
    my %found;
    foreach my $entry (@entries) {
        if (exists $found{$entry}) {
            die "Error, found $entry as both an ev_type to train, and as a known weight.  ";
        }

    }
    return; # all's good.
}

####
sub get_disjoint_A {
    my ($hashref_A, $hashref_B) = @_;
    
    unless (ref $hashref_A eq "HASH" && ref $hashref_B eq "HASH") {
        confess "Error, gimme two hashrefs here!";
    }

    my %disj_A;
    
    foreach my $key_A (keys %$hashref_A) {
        unless (exists $hashref_B->{$key_A}) {
            $disj_A{$key_A} = $hashref_A->{$key_A};
        }
    }

    return (%disj_A);
}


####
sub train_weights {
    my ($train_entries_aref, $EVM_cmd_template, $ev_types_aref, $fixed_weights_href, $lambdas_href) = @_;
    
    ## fixed weights are only fixed relative to one another. Their values can change below.

    &ensure_disjoint(@$ev_types_aref, keys %$fixed_weights_href);
    
    print "Training for ev_types: " . join (",", @$ev_types_aref) . "\n";
    
    ## Do random iteration to find best weights:
    my %best_weights = &iterate_random_weights($train_entries_aref, $EVM_cmd_template, $ev_types_aref, $fixed_weights_href);
    print "Got best weights from random iterations: " . Dumper(\%best_weights);
    
    ## do gradient descent:
    ## break into the variable and relative constant weights (better term?) 
    my %relative_constants;
    my %variable_weights;
    
    foreach my $ev_type (keys %best_weights) {
        my $weight = $best_weights{$ev_type};
        
        if (exists ($fixed_weights_href->{$ev_type})) {
            $relative_constants{$ev_type} = $weight; 
        }
        else {
            $variable_weights{$ev_type} = $weight;
        }
    }
    
    my %final_weights = &perform_gradient_descent($train_entries_aref, $EVM_cmd_template, \%variable_weights, \%relative_constants, $lambdas_href);
    
    return (%final_weights);
    
}

####
sub iterate_random_weights {
    my ($train_entries_aref, $EVM_cmd_template, $ev_types_aref, $fixed_weights_href) = @_;
    
    print "Iterating random weights.\n";
    &ensure_disjoint(@$ev_types_aref, keys %$fixed_weights_href);
    
    my @ev_types = @$ev_types_aref;
    my $num_pieces_evidence = scalar (@ev_types);
    my $weight_train_filename = "weight_training.random_weights.${num_pieces_evidence}_ev_types." . time() . ".log";
    open (my $logfh, ">$weight_train_filename") or die "Error, cannot write to $weight_train_filename";
    
    my $best_score = 0;
    my %best_weights = ();
    
    ## all weights sum to 1 (even though they're not probabilities)
    my @random_weight_combos_employed; ## track combinations used to ensure good spatial sampling
    
  RANDOM_ITER:
    for my $i (1..$num_rand_iterations) {
        my %random_weights;
        eval {
            print "\n==========================\n";
            print "$i. random weight iteration:\n";
            %random_weights = &get_random_weights($ev_types_aref, $fixed_weights_href, \@random_weight_combos_employed);
            
        };

        if ($@) {
            print "random sampling was apparently sufficient.  Exiting random iterations.\n";
            last RANDOM_ITER;
        }

        
        my ($BJHscore, $percent_gene_accuracy, $percent_exon_accuracy) = &score_weight_combination($train_entries_aref, $EVM_cmd_template, \%random_weights);
        print "** BJHscore: $BJHscore\n";
        
        ## log the result:
        &log_result($logfh, "RAND_mode", 
                    { BJHscore => $BJHscore, 
                      gene_accuracy => $percent_gene_accuracy,
                      exon_accuracy => $percent_exon_accuracy,
                  },
                    \%random_weights);
        
        if ($BJHscore > $best_score) {
            $best_score = $BJHscore;
            %best_weights = %random_weights;
        }
        
    }
    
    close $logfh;
    
    return (%best_weights);
}
        
####
sub log_result {
    my ($fh, $prefix, $atts_href, $weights_href) = @_;
    
    my $logtext = "$prefix\t" . &href_to_text($weights_href) . "\t" . &href_to_text($atts_href) . "\n"; 
    
    ## for log
    print $fh $logtext;

    # for stdout...
    $logtext =~ s/\t/\n\t/g;
    print "SUMMARY:\n$logtext\n\n\n"; 
    
    return;
}

####
sub href_to_text {
    my $href = shift;

    my $text = "";

    foreach my $key (keys %$href) {
        my $val = $href->{$key};
        
        if ($val =~ /^\d+\.\d+$/) {
            ## float
            $val = sprintf "%.4f", $val;
        }
        
        $text .= "$key=>$val,";
    }

    chop $text; #remove trail comma
    return ($text);
}


####
sub adjust_random_weights_sum_to_1 {
    my ($var_weights_href, $fixed_weights_href) = @_;
    
    if ($SEE) {
        print "adjusting weights sum to 1:\n";
        print "Var weights: " . Dumper ($var_weights_href);
        print "Fixed weights: " . Dumper ($fixed_weights_href);
    }
    
    &ensure_disjoint(keys %$var_weights_href, keys %$fixed_weights_href);
    

    my $curr_weight_sum = 0;
    
    my %all_weights;
    ## take into account known weights, and keep their relative proportions the same:
    foreach my $ev_type (keys %$var_weights_href) {
        my $weight = $var_weights_href->{$ev_type};
        $all_weights{$ev_type} = $weight;
        $curr_weight_sum += $weight;
    }

    ## add in the others:
    my $num_fixed_weights = scalar (keys %$fixed_weights_href);
    my $rand_val_contribution = $num_fixed_weights * rand(1); # fixed weights contribute an average random contribution, then adjusted by existing weight value.
    
    foreach my $ev_type (keys %$fixed_weights_href) {
        my $weight = $fixed_weights_href->{$ev_type} * $rand_val_contribution;
        $curr_weight_sum += $weight;
        $all_weights{$ev_type} = $weight;
    }
    

    ## normalize so sum to 1
    foreach my $ev_type (keys %all_weights) {
        my $weight = $all_weights{$ev_type};
        $all_weights{$ev_type} = $weight / $curr_weight_sum;
    }
    
    return (%all_weights);
}


####
sub perform_gradient_descent {
    my ($train_entries_aref, $EVM_cmd_template, $variable_weights_href, $relative_constants_href, $lambdas_href, $weights_and_score_container_aref) = @_;
    
    print "Starting gradient descent.\n";
    #print Dumper ($variable_weights_href);
    #print Dumper ($relative_constants_href);
    
    my $gradient_lambda = $lambdas_href->{gradient_lambda};
    my $max_lambda = $lambdas_href->{max_lambda};
        
    my @ev_types = reverse sort {$variable_weights_href->{$a}<=>$variable_weights_href->{$b}} keys %$variable_weights_href;
    
    ## start with the highest scoring ev_type, and walk the gradient up or down.
    
    my %trained_weights = (%$variable_weights_href, %$relative_constants_href);
    
    my $num_variable_weights = scalar (keys %$variable_weights_href);

    my $hill_climbing_filename = "grad_descent.${num_variable_weights}_var_weights.log";
    open (my $logfh, ">$hill_climbing_filename") or die "Error, cannot write to $hill_climbing_filename";
    
    my ($BJHscore, $percent_gene_accuracy, $percent_exon_accuracy) = &score_weight_combination($train_entries_aref, $EVM_cmd_template, \%trained_weights, $weights_and_score_container_aref); 
    &log_result($logfh, "INIT:hill_climbing", 
                { BJHscore => $BJHscore, 
                  gene_accuracy => $percent_gene_accuracy,
                  exon_accuracy => $percent_exon_accuracy,
              },
                \%trained_weights);
    
    my $best_BJHscore = $BJHscore;
    my %best_weights = %trained_weights;

    my $gradient_iteration = 0;
    my $gradient_starting_score = $best_BJHscore;

    do {
        ## keep going if making global score improvements.
        $gradient_starting_score = $best_BJHscore; ## set to current best score.  See if make improvements this round.
        
        $gradient_iteration++;

        
        foreach my $ev_type (@ev_types) {
        
            my $score_climbing_flag = 1;
            my $direction = "STATIC"; # not up or down.
            
            my $BJHscore_w_lambda_up = $BJHscore; #init
            my $BJHscore_w_lambda_down = $BJHscore; #init
            
            my ($percent_gene_accuracy, $percent_exon_accuracy,
                $percent_gene_accuracy_up, $percent_exon_accuracy_up,
                $percent_gene_accuracy_down,$percent_exon_accuracy_down);
            
            my $lambda_employed = $gradient_lambda; # init
            
            my $already_scored_zero_weight = 0;

            
          HILL_CLIMBING:
            while ($score_climbing_flag) {
                
                
                print "\n\n\nHILL CLIMBING EXERCISE ($ev_type), lamda_employed=$lambda_employed\n";
                
                ###############################
                ## Adjust weights up and score
                
                my $best_weight = $best_weights{$ev_type};
                
                my %adjusted_weights_up = %best_weights; # &get_disjoint_A(\%best_weights, $relative_constants_href); ## pull out the variable portion
                #print "adjusted weights up before: " . Dumper (\%adjusted_weights_up);
                
                $adjusted_weights_up{$ev_type} = $best_weight + $best_weight * $lambda_employed; ## adjust w/ lambda
                #%adjusted_weights_up = &adjust_random_weights_sum_to_1(\%adjusted_weights_up, $relative_constants_href); ## join back to the constant portion and adjust
                #print "adjusted weights up after: " . Dumper (\%adjusted_weights_up);
               
                
                if ($direction ne "DOWN") {
                    my $adj_weight_up = $adjusted_weights_up{$ev_type};
                    print "*Scoring UP-weighted ($ev_type, lambda: $lambda_employed, weight: $adj_weight_up vs. best: $best_weight):\n";
                    ($BJHscore_w_lambda_up, $percent_gene_accuracy_up, $percent_exon_accuracy_up) = &score_weight_combination($train_entries_aref, $EVM_cmd_template, \%adjusted_weights_up, $weights_and_score_container_aref);
                    &log_result($logfh, "hill_climb($ev_type,$lambda_employed):UP", 
                                { BJHscore => $BJHscore_w_lambda_up, 
                                  gene_accuracy => $percent_gene_accuracy_up,
                                  exon_accuracy => $percent_exon_accuracy_up,
                                  gradient_iter => $gradient_iteration,
                              },
                                \%adjusted_weights_up);
                }
                
                
                #################################
                ## Adjust weights down and score:
 
                my %adjusted_weights_down = %best_weights; # &get_disjoint_A(\%best_weights, $relative_constants_href);
                #print "adjusted weights down before: " . Dumper (\%adjusted_weights_down);
                $adjusted_weights_down{$ev_type} = $best_weight - $best_weight * $lambda_employed;
                if ($adjusted_weights_down{$ev_type} < 0) {
                    $adjusted_weights_down{$ev_type} = 0; # not penalizing evidence with negative weights
                }
                #%adjusted_weights_down = &adjust_random_weights_sum_to_1(\%adjusted_weights_down, $relative_constants_href);
                #print "adjusted weights down after: " . Dumper (\%adjusted_weights_down);
                


                if ($direction ne "UP") { 
                    my $adj_weight_down = $adjusted_weights_down{$ev_type};
                    if ( ! ($already_scored_zero_weight && $adj_weight_down == 0) )  {
                        print "*Scoring DOWN-weighted ($ev_type, lambda: $lambda_employed, weight: $adj_weight_down vs. best: $best_weight)\n";
                        ($BJHscore_w_lambda_down, $percent_gene_accuracy_down, $percent_exon_accuracy_down) = &score_weight_combination($train_entries_aref, $EVM_cmd_template, \%adjusted_weights_down, $weights_and_score_container_aref);
                        &log_result($logfh, "hill_climb($ev_type,$lambda_employed):DOWN", 
                                    { BJHscore => $BJHscore_w_lambda_down, 
                                      gene_accuracy => $percent_gene_accuracy_down,
                                      exon_accuracy => $percent_exon_accuracy_down,
                                      gradient_iter => $gradient_iteration,
                                  },
                                    \%adjusted_weights_down);
                        
                        if ($adj_weight_down == 0) {
                            $already_scored_zero_weight = 1;  ## don't keep computing this value during this ev_type round.
                        }
                    }
                }
                
                print "\n*******************************************************************************************************************************************\n";
                print "** UPvsDOWN-dir:$direction ($ev_type(w($best_weight),Lambda=$lambda_employed):   "
                    . "UP:{ $BJHscore_w_lambda_up }    "
                  . "DOWN:{ $BJHscore_w_lambda_down }  "
                . "BEST_SO_FAR: { $best_BJHscore }\n";
                print "********************************************************************************************************************************************\n\n";
                

                my $SIG_FIGS_FACTOR = 1e2; # if score is 0.95, would need the score to change by at least 0.01 to notice a difference.
                
                if (  ## No change observed in score
                      ( $direction eq 'UP' && round($BJHscore_w_lambda_up*$SIG_FIGS_FACTOR) == round($best_BJHscore*$SIG_FIGS_FACTOR) ) 
                      ||
                      ( $direction eq 'DOWN' && round($BJHscore_w_lambda_down*$SIG_FIGS_FACTOR) == round($best_BJHscore*$SIG_FIGS_FACTOR) )
                      ||
                      ($direction eq 'STATIC' 
                       && (  ## make sure we look outside of this range to determine an up or down movement.
                             round($BJHscore_w_lambda_up*$SIG_FIGS_FACTOR) <= round($best_BJHscore*$SIG_FIGS_FACTOR)  #stay in static mode until we start to see some improvement or hit the max lambda val.
                             &&  round($BJHscore_w_lambda_down*$SIG_FIGS_FACTOR) <= round($best_BJHscore*$SIG_FIGS_FACTOR)
                           ) 
                       ) 
                      )
                      
                    
                {   ## compare interger equality at 2 figures. (not interested in sublte variations below this).
                    ## lambda isn't big enough to notice any change
                    ## Try again with a larger lambda if possible.
                    
                    print "-lambda($lambda_employed) was too small to have effect. ";
                    $lambda_employed *= $LAMBDA_AMPLIFICATION_FACTOR;
                    print "Increasing lambda to $lambda_employed.\n";
                    if ($lambda_employed > $max_lambda) {
                        print "-lambda exeeds the max of $max_lambda.  skipping further attempts to adjust $ev_type\n";
                        last HILL_CLIMBING;
                    }
                    elsif ($direction eq 'DOWN' && $adjusted_weights_down{$ev_type} == 0) {
                        print "-adjusted weight for $ev_type in down direction was zero and cannot be made smaller.  moving on.\n";
                        last HILL_CLIMBING;
                    }
                    else {
                        next HILL_CLIMBING;
                    }
                }
                
                my ($best_direction, $best_local_BJHscore, %best_weights_local);
                if ($BJHscore_w_lambda_up > $BJHscore_w_lambda_down) {
                    $best_direction = "UP";
                    $best_local_BJHscore = $BJHscore_w_lambda_up;
                    %best_weights_local = %adjusted_weights_up;
                    $percent_gene_accuracy = $percent_gene_accuracy_up;
                    $percent_exon_accuracy = $percent_exon_accuracy_up;
                    
                    ## going up next round, so set down to this up, avoid recomputing it.
                    $BJHscore_w_lambda_down = $best_local_BJHscore; 
                }
                else {
                    $best_direction = "DOWN";
                    $best_local_BJHscore = $BJHscore_w_lambda_down;
                    %best_weights_local = %adjusted_weights_down;
                    $percent_gene_accuracy = $percent_gene_accuracy_down;
                    $percent_exon_accuracy = $percent_exon_accuracy_down;
                    
                    ## going down next round, so set up to this down, avoid recomputing it.
                    $BJHscore_w_lambda_up = $best_local_BJHscore;
                }
                
                if ($direction ne 'STATIC' && $direction ne $best_direction) {
                    print "-best direction switched from $direction to $best_direction.  Stopping climbing/descent\n";
                    $score_climbing_flag = 0; # stop climbing if switched directions
                }
                
                my $prev_direction = $direction;
                $direction = $best_direction;
                
                if ($best_local_BJHscore > $best_BJHscore) {
                    $best_BJHscore = $best_local_BJHscore;
                    %best_weights = %best_weights_local;
                    
                    print "***** IMPROVING ON ACCURACY ($ev_type, $best_direction) *****\n";  
                    &log_result($logfh, "hill_climb_IMPROVEMENT($ev_type):$best_direction", 
                                { BJHscore => $best_BJHscore, 
                                  gene_accuracy => $percent_gene_accuracy,
                                  exon_accuracy => $percent_exon_accuracy,
                                  gradient_iter => $gradient_iteration,
                              },
                                \%best_weights);
                }
                else {
                
                    if ($prev_direction eq 'STATIC') {
                        print "** Switching from STATIC to direction $direction\n";
                        # not expecting a better score yet.  Only determined direction.
                        next HILL_CLIMBING;
                    }
                    else {
                        ## if haven't found a better score, then we're done looking for one here.
                        print "-haven't found a higher score by manipulating the $ev_type weight.  Moving on.\n";
                        $score_climbing_flag = 0;
                    }
                }
            }
        }
        
    } while ($best_BJHscore - $gradient_starting_score >= $required_gradient_improvement_points 
             && $num_variable_weights > 1 );
    
    
    return (%best_weights);
}

sub print_weight_values {
    my ($fh, $weights_href) = @_;
    
    my @ev_types = reverse sort {$weights_href->{$a}<=>$weights_href->{$b}} keys %$weights_href;
    
    my $sum = 0;
    foreach my $ev_type (@ev_types) {
        my $weight_value = $weights_href->{$ev_type};
        print $fh "\t$ev_type\t$weight_value\n";
        $sum += $weight_value;
    }
    return ($sum);
}




####
sub score_weight_combination {
    my ($train_entries_aref, $EVM_cmd_template, $weights_href, $weights_and_score_container_aref) = @_;
    
    $EVM_cmd_template .= " --INTERGENIC_SCORE_ADJUST_FACTOR $INTERGENIC_SCORE_ADJUST_FACTOR ";
    
    print "Scoring weights:  (intergenic factor: $INTERGENIC_SCORE_ADJUST_FACTOR)\n";
    my $sum = &print_weight_values(*STDOUT, $weights_href);
    print "\t----------------\n";
    print "\tSUM:\t$sum\n";
    
    my $score_begin_time = time();

    ## write the weights file:
    my $weights_file = "current_weights.txt";
    open (my $fh, ">$weights_file") or die $!;
    foreach my $ev_type (keys %$weights_href) {
        my $weight = $weights_href->{$ev_type};
        my $class = $EVIDENCE_CLASS{$ev_type};
        
        unless (defined $weight && defined $class) {
            confess "Error, no weight($weight) or class($class) for ev_type $ev_type ";
        }
        
        print $fh "$class\t$ev_type\t$weight\n";
    }
    close $fh;

    ## Stub, iterate over each, compute F-score
    
    $EVM_cmd_template .= " -w $weights_file ";

    $EVM_cmd_template .= " > __PATH__/evm.out 2>/dev/null";


    my $EVM_to_GFF3_cmd_template = "$bindir/../EvmUtils/EVM_to_GFF3.pl __PATH__/evm.out __GENOME_ACC__ > __PATH__/evm.gff3";
    
    my $snsp_analyzer = new SnSp::SnSp_analysis_manager("EVM");
    $snsp_analyzer->set_intergenic_included($FLANK_REGIONS_SnSp_BP);
    
    my $token = '__PATH__';

    my @cmds_to_run_on_grid;

    foreach my $train_entry (@$train_entries_aref) {
        
        #print "EVM exec ($train_entry)\n";
        
                
        my $evm_cmd = $EVM_cmd_template;
        
        my $full_path = $workdir . "/$train_entry";

        $evm_cmd =~ s/$token/$full_path/g;
        #&process_cmd($evm_cmd);
        
        # save cmd executed
        open (my $cmd_fh, ">$train_entry/cmd.exec") or die "Error, cannot open file for writing: ($train_entry/cmd.exec) \n";
        print $cmd_fh $evm_cmd;
        close $cmd_fh;

        push (@cmds_to_run_on_grid, $evm_cmd);

    }
    
    my $num_cmds = scalar (@cmds_to_run_on_grid);
    my $start_time = time();
    ## execute cmds on grid or locally:
    my $count = 0;
    
    if ($use_grid_computing_flag) {
        &run_commands_on_grid(@cmds_to_run_on_grid);
    }
    else {
        
        my $time_per_exec = 0;
        foreach my $cmd (@cmds_to_run_on_grid) { 
            $count++;
            if ($SEE) {
                print "CMD[$count] $cmd\n";
            }
            else {
                print "\rexec cmd $count of $num_cmds (prev took $time_per_exec s)    ";
            }
            my $time = time();
            &process_cmd($cmd);
            $time_per_exec = time() - $time;
        }
        print "done running cmds locally\n";
    }
    
    my $finish_time = time();
    
    my $num_seconds_runtime = $finish_time - $start_time;
    my $num_minutes_runtime = $num_seconds_runtime / 60;

    my $num_seconds_per_exec = $num_seconds_runtime / $num_cmds;

    printf("*$num_cmds EVM cmds took %.2f minutes = %.2f seconds per exec.\n", $num_minutes_runtime, $num_seconds_per_exec);
        
    open (my $gene_acc_fh, ">gene_accuracy_from_current_weights.summary") or die $!;

    foreach my $train_entry (@$train_entries_aref) {
        
        ## covert output to gff3 format first:
        my $evm_to_gff3_cmd = $EVM_to_GFF3_cmd_template;
        $evm_to_gff3_cmd =~ s/$token/$train_entry/g;
        
        ## get the accession:
        my $genome_fasta_file = $train_entry . "/$genomicSeqFile";
        my $header = `head -n1 $genome_fasta_file`;
        $header =~ />(\S+)/;
        my $genome_acc = $1 or die "Error, can't get accession from header of genome file $genome_fasta_file";
        $evm_to_gff3_cmd =~ s/__GENOME_ACC__/$genome_acc/;
        
        &process_cmd($evm_to_gff3_cmd);
        
        my %ev_type_to_genes = ( "EVM" => [] );
        
        ## append results to the SnSp analysis:
        
        my $tfpn_obj = $snsp_analyzer->get_TFPN("EVM");
        
        my $genes_correct_so_far = $tfpn_obj->{genes_predicted_correct};

        &do_SnSp_analysis($snsp_analyzer, "$train_entry/template.gff3", ["$train_entry/evm.gff3"], \%ev_type_to_genes);
        
        my $genes_correct_now = $tfpn_obj->{genes_predicted_correct};
        
        my $correct_flag = "";
        if ($genes_correct_now > $genes_correct_so_far) {
            $correct_flag = "CORRECT";
        }
        else {
            $correct_flag = "WRONG";
        }
                
        print $gene_acc_fh "$train_entry\t$correct_flag\n";

    }
    
    close $gene_acc_fh;

    my $TFPN_obj = $snsp_analyzer->get_TFPN("EVM");
    my $Fscore = $TFPN_obj->{Fscore};
    
    my $num_genes_analyzed = $snsp_analyzer->{number_genes_analyzed};
    my $num_exons_analyzed = $snsp_analyzer->{number_exons_analyzed};
    
    my $num_genes_correct = $TFPN_obj->{genes_predicted_correct};
    my $num_exons_correct = $TFPN_obj->{number_exons_correct};
    
    my $sensitivity = $TFPN_obj->{sensitivity};
    my $specificity = $TFPN_obj->{specificity};
    
    my $percent_gene_accuracy = $num_genes_correct / $num_genes_analyzed * 100;
    my $percent_exon_accuracy = $num_exons_correct / $num_exons_analyzed * 100;

    
    
    my $BJHscore = $Fscore + $percent_gene_accuracy/100 + $percent_exon_accuracy/100; ## modified scoring function (now the BJHscore instead of the Fscore alone) 
    
    print "### BJH: $BJHscore\n"
        . "### Fscore: $Fscore\n"
        . "### genes: ($num_genes_correct/$num_genes_analyzed) = $percent_gene_accuracy\n"
        . "### exons: ($num_exons_correct/$num_exons_analyzed) = $percent_exon_accuracy\n"
        . "### NSn: $sensitivity\n"
        . "### NSp: $specificity\n";
        
    
    my $score_end_time = time();
    my $full_scoring_time = ($score_end_time - $score_begin_time) / 60;
    
    printf("Scoring round took %.2f minutes\n", $full_scoring_time);
    

    if (ref $weights_and_score_container_aref eq "ARRAY") {
        push (@$weights_and_score_container_aref, { weights => { %$weights_href },
                                                    score => $BJHscore } );
    }
    
    
    return ($BJHscore, $percent_gene_accuracy, $percent_exon_accuracy);
}


####
sub process_cmd {
    my $cmd = shift;
    
    my $ret = system $cmd;
    if ($ret) {
        die "Error, cmd ($cmd) died with ret ($ret)\n";
    }
    return;
}

####
sub do_SnSp_analysis {
    my ($snsp_analyzer, $template_gff3_file, $other_gff3_files_aref, $ev_type_to_genes_href, $verbose_flag) = @_;

    
    my $gene_obj_indexer = {}; # using hash reference instead.
    &GFF3_utils::index_GFF3_gene_objs($template_gff3_file, $gene_obj_indexer);
    
    my ($template_gene_id) = keys %$gene_obj_indexer;
    my $template_gene = $gene_obj_indexer->{$template_gene_id};
    
    my ($template_lend, $template_rend) = sort {$a<=>$b} $template_gene->get_coords();
    
    # print "template span: $template_lend - $template_rend\n";
    
    ## get the others:
    foreach my $other_gff3_file (@$other_gff3_files_aref) {
        
        my $other_gene_obj_indexer = {}; # using hash reference instead.
        &GFF3_utils::index_GFF3_gene_objs($other_gff3_file, $other_gene_obj_indexer);
        
        my @gene_ids = keys %$other_gene_obj_indexer;
        
        foreach my $gene_id (@gene_ids) {
            my $gene_obj = $other_gene_obj_indexer->{$gene_id};
            
            my $source = $gene_obj->{source} || confess "Error, need the undocumented 'source' attribute for the gene_obj";
            
            #print "Source: $source\n";
            my $gene_list_aref = $ev_type_to_genes_href->{$source};
            unless (ref $gene_list_aref) { 
                #print "sorry, not interested in source $source\n";
                next;
            }
            
            
            ## only include it if it overlaps the coordinates of the template
            my ($lend, $rend) = sort {$a<=>$b} $gene_obj->get_coords();

            #print "\tother: $gene_id, $lend - $rend\n";
            unless ($lend < $template_rend && $rend > $template_lend) {
                #print "\tsorry, no overlap to template.\n";
                next; 
            }
            
            foreach my $gene_entry ($gene_obj, $gene_obj->get_additional_isoforms()) {
                $gene_entry->delete_isoforms(); ## unwrapping them for the SnSp computer.
                push (@$gene_list_aref, $gene_entry);
            }
        }
    }
        
    $snsp_analyzer->add_analysis_entry($template_gene, $ev_type_to_genes_href, $verbose_flag);
    
    return;
}


####
sub get_random_weights {
    my ($ev_types_aref, $fixed_weights_href, $random_weight_combos_employed_aref) = @_;


    my %random_weights;
    
    my $found_suitable_random_weights_flag = 0;

    my $try_count = 0;

    while (! $found_suitable_random_weights_flag) {
        
        $try_count++;
        
        ## try a random weight combination:
        foreach my $ev_type (@$ev_types_aref) {
            my $rand_weight = rand(1);
            $random_weights{$ev_type} = $rand_weight;
        }
        
        %random_weights = &adjust_random_weights_sum_to_1(\%random_weights, $fixed_weights_href);
        
        if (&lacks_significant_overlap_to_random_weights_tried(\%random_weights, $random_weight_combos_employed_aref)) {
            $found_suitable_random_weights_flag = 1;
        }
        else {
            if ($try_count > $MAX_RAND_ITERATION_TRY_COUNT) {
                die "tried max random iterations and none within min euclidean distance."; # break the eval (exception thrown here)
            }
        }
    }
        
    return (%random_weights);

}

####
sub lacks_significant_overlap_to_random_weights_tried {
    my ($random_weights_href, $random_weight_combos_employed_aref) = @_;

    unless (@$random_weight_combos_employed_aref) {
        ## nothing there, so take these random weights as first used data set.
        push (@$random_weight_combos_employed_aref, {%$random_weights_href});
        return (1);
    }

    ## check the Euclidean distance between each previously used entry and the current set.
    my $min_euclidian_dist = undef;
    

    foreach my $random_weights_used_href (@$random_weight_combos_employed_aref) {
        my $euclid_dist = &calc_Euclidian_distance($random_weights_used_href, $random_weights_href);
        
        # print "-euclid_dist: $euclid_dist\n";
        

        if (! defined $min_euclidian_dist) {
            $min_euclidian_dist = $euclid_dist;
        }
        else {
            if ($euclid_dist < $min_euclidian_dist) {
                $min_euclidian_dist = $euclid_dist;
            }
        }
    }
        
    print "min_euclid_dist: $min_euclidian_dist\n";
    
    if ($min_euclidian_dist >= $MIN_EUCLIDIAN_DISTANCE) {
        push (@$random_weight_combos_employed_aref, {%$random_weights_href});
        return (1); #OK
    }
    else {
        return (0); # too similar to a previously used random weight combo.
    }
}

####
sub calc_Euclidian_distance {
    my ($weights_A_href, $weights_B_href) = @_;
    
    ## should have same set of keys!
    if (&get_disjoint_A($weights_A_href, $weights_B_href) || &get_disjoint_A($weights_B_href, $weights_A_href)) {
        ## some keys in A and not in B or vice-versa!
        print "A-weights: " . Dumper ($weights_A_href);
        print "B-weights: " . Dumper ($weights_B_href);
        
        confess "Error, keys in weights hrefs are not identical.";
    }

    my @ev_types = keys %$weights_A_href;

    ## calc the RMS deviation as the Euclid Dist
    
    my $sum_squared_diff = 0;
    foreach my $ev_type (@ev_types) {
        my $weight_A = $weights_A_href->{$ev_type};
        my $weight_B = $weights_B_href->{$ev_type};

        my $delta = $weight_A - $weight_B;

        my $delta_squared = $delta ** 2;
        
        $sum_squared_diff += $delta_squared;
    }

    my $rms_dev = $sum_squared_diff ** 0.5;

    return ($rms_dev);
}
    
####
sub populate_weights_from_file {
    my ($filename, $weights_href) = @_;
    open (my $fh, $filename) or die "Error, cannot open file $filename";
    while (<$fh>) {
        unless (/\w/) { next; }
        if (/^\#/) { next; } #commented weights entry
        chomp;
        my ($ev_class, $ev_type, $weight) = split (/\s+/);
        $weights_href->{$ev_type} = $weight;
        $EVIDENCE_CLASS{$ev_type} = $ev_class;
    }
    close $fh;
}


####
sub find_maximum_weight_for_ev_type {
    my ($train_entries_aref, $EVM_cmd_template, $fixed_weights_href, $ev_type) = @_;
    
    open (my $logfh, ">maximizing_${ev_type}_weight.log");

    my $init_weight = 0;
    # set initial weight to the highest weight for any other ev_type so far:
    foreach my $ev (keys %$fixed_weights_href) {
        my $weight = $fixed_weights_href->{$ev};
        if ($weight > $init_weight) {
            $init_weight = $weight;
        }
    }
    
    my %applied_weights = %$fixed_weights_href;
    
    my %best_weights;

    my $best_score = 0;
    
    my $max_weight_value = 10 * $init_weight + $init_weight/2; # place an upper limit on the highest weight.
    
    for (my $weight = 0; $weight <= $max_weight_value; $weight += $init_weight) {
        
        $applied_weights{$ev_type} = $weight;
        
        my ($BJHscore, $percent_gene_accuracy, $percent_exon_accuracy)  = &score_weight_combination($train_entries_aref, $EVM_cmd_template, \%applied_weights);
            
        my $score = $BJHscore;
        if ($score >= $best_score) {
            $best_score = $score;
            %best_weights = %applied_weights;
        }
        
        
        &log_result($logfh, "[$ev_type: $weight]",
                    { BJHscore => $BJHscore, 
                      gene_accuracy => $percent_gene_accuracy,
                      exon_accuracy => $percent_exon_accuracy,
                    },
                    \%applied_weights,
            );
    }
    
    close $logfh;
    
    return (%best_weights);
    
    
}



####
sub measure_other_prediction_accuracies {
    my ($train_entries_aref, $prediction_prog_list_aref) = @_;

    my $snsp_analyzer = new SnSp::SnSp_analysis_manager(@$prediction_prog_list_aref);
    $snsp_analyzer->set_intergenic_included($FLANK_REGIONS_SnSp_BP);
    
    foreach my $train_entry (@$train_entries_aref) {
        print "\n\n----\nEntry: $train_entry\n";
        my %ev_type_to_genes;
        foreach my $ev_type (@$prediction_prog_list_aref) {
            $ev_type_to_genes{$ev_type} = []; #init to empty list
        }
        &do_SnSp_analysis($snsp_analyzer, "$train_entry/template.gff3", ["$train_entry/$genePredictionsFile"], \%ev_type_to_genes, 1);
    }
    
    return;
    
}

####
sub run_commands_on_grid {
    my @cmds = @_;

    unless ($READY_FOR_GRID_SUBMISSION) { 
        &prepare_for_grid_usage(scalar (@cmds));
    }
    &Run_Bsub::set_queue("broad");
    my $ret = &Run_Bsub::run(@cmds);

    if ($ret) {
        die "Error, at least one grid job failed. Ret: $ret\n";
    }
    else {
        print "Grid jobs completed successfully\n";
    }

    return;
    
}

####
sub write_weights_file {
    my ($filename, $weights_href) = @_;

    open (my $fh, ">$filename") or die "Error, cannot write file $filename";
    foreach my $ev_type (reverse sort {$weights_href->{$a}<=>$weights_href->{$b}} keys %$weights_href) {
        my $weight = $weights_href->{$ev_type};
        my $class = $EVIDENCE_CLASS{$ev_type};
        print $fh "$class\t$ev_type\t$weight\n";
    }
    close $fh;

    return;
}


sub round {
    my ($num) = shift;
    return ( int($num + 0.5) );
}

####
sub prepare_for_grid_usage {
    my ($num_cmds) = shift;

    unless ($cmds_per_node) {
        $cmds_per_node = int($num_cmds / 100);
        if ($cmds_per_node < 1) { 
            $cmds_per_node = 1;
        }
    }
    
    ## load the institute-specific modules for performing grid submissions.
    
    eval {
        push (@INC, $ENV{EUK_MODULES});
        require "$ENV{EUK_MODULES}/Run_Bsub.pm";
        my $packagename = "Run_Bsub";
        import $packagename;
        
        ## set up grid to run one cmd per node.
        $Run_Qsub::CMDS_PER_NODE = $cmds_per_node; 
        
    };
    if ($@) {
        die "Error, unable to load the Run_Bsub.pm module for grid submissions";
    }
    else {
        $READY_FOR_GRID_SUBMISSION = 1;
    }
    
    return;
}


####
sub compute_intergenic_score_adjust_factor {
    my ($train_entries_aref, $EVM_cmd_template, $prediction_progs_aref) = @_;

    my @scores_and_intergenic_factor;

    ## weight each prediction type equally for now.
    my %weights;
    foreach my $prog (@$prediction_progs_aref) {
        $weights{$prog} = 1;
    }
    
    for (my $i = 1; $i >= 0; $i -= 0.1) {
       
        $INTERGENIC_SCORE_ADJUST_FACTOR = $i;
 
        my ($BJHscore, $percent_gene_accuracy, $percent_exon_accuracy) = &score_weight_combination($train_entries_aref, $EVM_cmd_template, \%weights);
        push (@scores_and_intergenic_factor, { score => $BJHscore,
                                               factor => $i,
                                           } );
    }


    ## set the factor to the highest value that provides the optimum score:
    @scores_and_intergenic_factor = sort {$b->{score}<=>$a->{score}
                                          ||
                                              $b->{factor} <=> $a->{factor} } @scores_and_intergenic_factor;

    my $best = shift @scores_and_intergenic_factor;
    my $factor = $best->{factor};

    print "Best intergenic factor found to be: $factor\n";
    $INTERGENIC_SCORE_ADJUST_FACTOR = $factor;

    return;
}
