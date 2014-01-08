#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use Getopt::Long;
use Data::Dumper;

my $combo_dir = "combo_outs";

$|++;

my $params = <<__PARAMS__;

######################################################################################
#
#  Required:
# 
#  --genome               genome sequence in fasta format
#  --weights              weights for evidence types file
#  --gene_predictions       |  gene predictions gff3 file
#  --evaluationEntries             file containing the list of entries for evaulation
#
#  Optional but recommended:
#  --protein_alignments   protein alignments gff3 file
#  --transcript_alignments      transcript alignments gff3 file
#
#  Optional and miscellaneous
#
#  --repeats             gff3 file with repeats masked from genome file
#
#  
#  --terminalExons      supplementary file of additional terminal exons to consider (from PASA long-orfs)
#
#  --stop_codons             list of stop codons (default: TAA,TGA,TAG)
#                               *for Tetrahymena, set --stop_codons TGA
#  --min_intron_length       minimum length for an intron (default 20 bp)
#                                                                       
#  misc:
#   --use_grid_computing   (tigr-only for now!)
#
#   --NO_RECURSE              only do a single traceback along the sequence.
#          
#######################################################################################

__PARAMS__

    ;

my ($genomicSeqFile, $genePredictionsFile, $transcriptAlignmentsFile, $proteinAlignmentsFile,
    $repeatsFile, $weightsFile, $terminalExonsFile, $stop_codons_arg, $help_flag,
    $stitch_ends, $extend_to_terminal, $evaluationEntries, $min_intron_length,
    $use_grid_computing, $NO_RECURSE_FLAG,
    );


&GetOptions ("genome=s"=>\$genomicSeqFile,
             "gene_predictions=s"=>\$genePredictionsFile,
             "transcript_alignments=s"=>\$transcriptAlignmentsFile,
             "protein_alignments=s" => \$proteinAlignmentsFile,
             "repeats=s"=>\$repeatsFile,
             "weights=s"=>\$weightsFile,
             "terminalExonsFile=s"=>\$terminalExonsFile,
             "stop_codons=s" => \$stop_codons_arg,
             "help|h" => \$help_flag,
             
             "evaluationEntries=s" => \$evaluationEntries,
             "min_intron_length=i" => \$min_intron_length,
             "use_grid_computing" => \$use_grid_computing,
             "NO_RECURSE" => \$NO_RECURSE_FLAG,
             );                                         

unless ($genomicSeqFile && $genePredictionsFile && $weightsFile && $evaluationEntries) {
    die $params;
}

my $FLANK_REGIONS_SnSp_BP = 500; # default for regions flanking the reference gene structure to be used in SnSp calculations.



my @evaluation_entries_list = `cat $evaluationEntries`;
chomp @evaluation_entries_list;

my %weights;
my %ev_classes;

my (@genefinders, @proteins, @transcripts, @pasa, @genewise, @homology_predictions);

{ # parse weights file
    open (my $fh, $weightsFile) or die $!;
    while (<$fh>) {
        unless (/\w/) { next; }
        if (/^\#/) {
            next; 
        }
        chomp;
        my ($ev_class, $ev_type, $weight) = split (/\s+/);
        $weights{$ev_type} = $weight;
        $ev_classes{$ev_type} = $ev_class;
    }
    
    foreach my $ev_type (keys %ev_classes) {
        
        print "[$ev_type]\n";
        
        my $ev_class = $ev_classes{$ev_type};
        
        if ($ev_class eq "ABINITIO_PREDICTION") {
            push (@genefinders, $ev_type);
        }
        elsif ($ev_class eq 'OTHER_PREDICTION') {
            push (@homology_predictions, $ev_type);
        }
        
        elsif ($ev_class eq "TRANSCRIPT") {
            if ($ev_type =~ /^alignassembly/i) {
                push (@pasa, $ev_type);
            }
            else {
                push (@transcripts, $ev_type);
            }
        }
        elsif ($ev_class eq "PROTEIN") {
            if ($ev_type =~ /^genewise/i) {
                push (@genewise, $ev_type);
            }
            else {
                push (@proteins, $ev_type);
            }
        }

        else {
            die "Error, I don't know how to handle ev_type: $ev_type of class: $ev_class\n";
        }
    }
}


print "genefinders: " . Dumper (\@genefinders) . "\n"
    . "homology predictions: " . Dumper (\@homology_predictions) . "\n"
    . "proteins: " . Dumper (\@proteins) . "\n"
    . "transcripts: " . Dumper (\@transcripts) . "\n"
    . "pasa: " . Dumper (\@pasa) . "\n"
    . "genewise: " . Dumper (\@genewise) . "\n";


main: {
    
    &run_cmds();
    
    &summarize_results();

    exit(0);
    
}


####
sub run_cmds {
    if (-d $combo_dir) {
        system "rm -r $combo_dir";
    }
    
    mkdir $combo_dir or die "Error, cannotmkdir $combo_dir";

    my $bindir = $FindBin::Bin;
    
    my $core_cmd = "$bindir/train_weights.pl --genome $genomicSeqFile --gene_predictions $genePredictionsFile --training_entries $evaluationEntries "
        . " --just_score_existing_weights -w weights.test " ;
    if ($use_grid_computing) {
        $core_cmd .= " --use_grid_computing";
    }
    if ($NO_RECURSE_FLAG) {
        $core_cmd .= " --NO_RECURSE";
    }
    
    my @processes = (   
        # 1
        { 
            type => "just_genefinders",
            cmd => $core_cmd,
            weight_types => [@genefinders],
        },
                        
                        );

    if (@homology_predictions) {
        push (@processes, 
              {
                  type => "just_homology_preds",
                  cmd => $core_cmd,
                  weight_types => [@homology_predictions],
              }
              );
    

        push (@processes, 

              { type => "genefinders_plus_homology_preds",
                cmd => $core_cmd,
                weight_types => [@genefinders, @homology_predictions],
            }
              );

        print "*** Note, from here on bundling all homology predictions with the genefinders ***\n";
        push (@genefinders, @homology_predictions);
        
    }
    
    
    if ($proteinAlignmentsFile) {
        push (@processes, 
              # 2
              { 
                  type => "proteinsAAT",
                  cmd => $core_cmd . " --protein_alignments $proteinAlignmentsFile",
                  weight_types => [@genefinders, @proteins],
              },
              );
        
        ## try extending protein alignment termini:
        #push (@processes, 
        #      # 2
        #      { 
        #          type => "proteinsAAT--extend",
        #          cmd => $core_cmd . " --protein_alignments $proteinAlignmentsFile",
        #          weight_types => [@genefinders, @proteins],
        #          extend_list => [@proteins],
        #      },
        #      );
        
        

    }
    
    
    # 3
    if ($transcriptAlignmentsFile) {
        push (@processes, 
              { 
                  type => "estsAAT",
                  cmd => $core_cmd . " --transcript_alignments $transcriptAlignmentsFile",
                  weight_types => [@genefinders, @transcripts],
              },
              );
        
        ## try extending the transcripts:
        #push (@processes, 
        #      { 
        #          type => "estsAAT--extend",
        #          cmd => $core_cmd . " --transcript_alignments $transcriptAlignmentsFile",
        #          weight_types => [@genefinders, @transcripts],
        #          extend_list => [@transcripts],
        #      },
        #      );
    }
    
    if ($proteinAlignmentsFile && $transcriptAlignmentsFile) {
        # 4
        
        push (@processes, 
              { 
                  type => "proteinsAAT_estsAAT",
                  cmd => $core_cmd . " --transcript_alignments $transcriptAlignmentsFile --protein_alignments $proteinAlignmentsFile",
                  weight_types => [@genefinders, @proteins, @transcripts],
              },
              
              #{ 
              #    type => "proteinsAAT_estsAAT-extend",
              #    cmd => $core_cmd . " --transcript_alignments $transcriptAlignmentsFile --protein_alignments $proteinAlignmentsFile",
              #    weight_types => [@genefinders, @proteins, @transcripts],
              #    extend_list => [@proteins, @transcripts]
              #},
              
              );
        


        
        
        if (@genewise) {
            push (@processes, 
                  { 
                      type => "proteinsAAT_estsAAT_genewise",
                      cmd => $core_cmd . " --transcript_alignments $transcriptAlignmentsFile --protein_alignments $proteinAlignmentsFile",
                      weight_types => [@genefinders, @proteins, @transcripts, @genewise],
                  },

               #   { 
               #       type => "proteinsAAT_estsAAT_genewise-extend",
               #       cmd => $core_cmd . " --transcript_alignments $transcriptAlignmentsFile --protein_alignments $proteinAlignmentsFile",
               #       weight_types => [@genefinders, @proteins, @transcripts, @genewise],
               #       extend_list => [@proteins, @transcripts, @genewise],
               #   },

                  
               #   { 
               #       type => "proteinsAAT_estsAAT_genewise-ESTs-extend-only",
               #       cmd => $core_cmd . " --transcript_alignments $transcriptAlignmentsFile --protein_alignments $proteinAlignmentsFile",
               #       weight_types => [@genefinders, @proteins, @transcripts, @genewise],
               #       extend_list => [@transcripts,],
               #   },
                  
                  
                  );
        }
    
    }
    
    
    #5
    if ($proteinAlignmentsFile && @genewise) {
        
        push (@processes, 

              {
                  type => "proteinsAAT_genewise",
                  cmd => $core_cmd . " --protein_alignments $proteinAlignmentsFile",
                  weight_types => [@genefinders, @proteins, @genewise],
              },
              

              # 6 
              
              #{ 
              #    type => "proteinsAAT_genewise-extend",
              #    cmd => $core_cmd . " --protein_alignments $proteinAlignmentsFile ",
              #    weight_types => [@genefinders, @proteins, @genewise],
              #    extend_list => [@genewise],
              #},
              
              {
                  type => "genewise",
                  cmd => $core_cmd . " --protein_alignments $proteinAlignmentsFile ",
                  weight_types => [@genefinders, @genewise],
              },
              
              #{ 
              #    type => "genewise-extend",
              #    cmd => $core_cmd . " --protein_alignments $proteinAlignmentsFile ",
              #    weight_types => [@genefinders, @genewise],
              #    extend_list => [@genewise],
              #},
              
              
              );
    }
    
    if ($transcriptAlignmentsFile && @pasa) {
        
        #push (@processes, 
        #      # 9
        #      #{ 
        #      #    type => "pasa",
        #      #    cmd => $core_cmd . " --transcript_alignments $transcriptAlignmentsFile",
        #      #    weight_types => [@genefinders, @pasa],
        #      #},
        #      #
        #
        #          { 
        #          type => "pasa-extend",
        #          cmd => $core_cmd . " --transcript_alignments $transcriptAlignmentsFile",
        #          weight_types => [@genefinders, @pasa],
        #          extend_list => [@pasa],
        #      },
        #      
        #      );
        
        if ($terminalExonsFile) {
            push (@processes, 
        
                  # 11
                  {
                      type => "pasa_term",
                      cmd => $core_cmd . " --transcript_alignments $transcriptAlignmentsFile --terminalExons $terminalExonsFile",
                      weight_types => [@genefinders, @pasa],
                  },
                  
                  
                  );
        
            if ($proteinAlignmentsFile) {
                
                push (@processes, 
                      # 13

                       {
                          type => "all_evidence_pasa_term",
                          cmd => $core_cmd . " --protein_alignments $proteinAlignmentsFile --transcript_alignments $transcriptAlignmentsFile --terminalExons $terminalExonsFile ",
                          weight_types => [@genefinders, @transcripts, @pasa, @genewise, @proteins],
                      },
    
                      

                      #{
                      #    type => "all_evidence_pasa_term-extend",
                      #    cmd => $core_cmd . " --protein_alignments $proteinAlignmentsFile --transcript_alignments $transcriptAlignmentsFile --terminalExons $terminalExonsFile ",
                      #    weight_types => [@genefinders, @transcripts, @pasa, @genewise, @proteins],
                      #    extend_list => [@genewise, @transcripts, @genewise, @proteins],
                      #},
                      #
                      
                      
                      );
            }
            
        }
    }
    
    ####
    foreach my $cmdset (@processes) {
        my ($cmd, $type, $weight_types_aref)  = ($cmdset->{cmd}, $cmdset->{type}, $cmdset->{weight_types});

        print "Cmdset: " . Dumper ($cmdset);
   
        ## create a weights file:
        my $weights_file = "weights.test";
        open (my $fh, ">$weights_file") or die $!;
        foreach my $type (@$weight_types_aref) {
            my $weight = $weights{$type};
            my $class = $ev_classes{$type};
            unless (defined $weight) { 
                die "Error, no weight for $type\n";
            }
            print $fh "$class\t$type\t$weight\n";
        }
        close $fh;
        
        if (my $stitch_list = $cmdset->{stitch_list}) {
            open (my $fh, ">stitch_list_file") or die $!;
            foreach my $stitch_type (@$stitch_list) {
                print $fh $ev_classes{$stitch_type} . "\t$stitch_type\n";
            }
            close $fh;
            $cmd .= " --stitch_ends stitch_list_file ";
        }
        
        if (my $extend_list = $cmdset->{extend_list}) {
            open (my $fh, ">extend_list_file") or die $!;
            foreach my $extend_type (@$extend_list) {
                print $fh $ev_classes{$extend_type} . "\t$extend_type\n";
            }
            close $fh;
            $cmd .= " --extend_to_terminal extend_list_file ";
        }
        
        ## run it:
        print "\n\n** Running: $type\nCMD: $cmd\n";
        my $result = `$cmd`;
        if ($?) {
            die "Error, cmd $cmd died with ret ($?)\n";
        }
        
        open ($fh,  ">$combo_dir/type.$type") or die $!;
        print $fh $result;

        rename("gene_accuracy_from_current_weights.summary", "$combo_dir/evaluation.$type") or die "could not rename file $!";
        

        ## archive the outputs!
        if (-e "$combo_dir/$type.evm.gff3") { die "Error, archive for evm.gff3 already exists at $combo_dir/$type.evm.gff3";}
        if (-e "$combo_dir/$type.evm.out") { die "Error, archive for evm.out already exists at $combo_dir/$type.evm.out"; }
        
        foreach my $eval_entry (@evaluation_entries_list) {
            
            my $evm_gff3_file = "$eval_entry/evm.gff3";
            unless (-e $evm_gff3_file) {
                die "Error, cannot find file $evm_gff3_file";
            }
            my $ret = system ("cat $evm_gff3_file >> $combo_dir/$type.evm.gff3");
            
            if ($ret) {
                die "gff3 aggregation failed for $type, file $evm_gff3_file";
            }
            
            my $evm_out_file = "$eval_entry/evm.out";
            unless (-e $evm_out_file) {
                die "Error, cannot locate file $evm_out_file";
            }
            $ret = system "echo !! Entry: $eval_entry, type $type >> $combo_dir/$type.evm.out";
            if ($ret) { 
                die "Error echoing comment to $combo_dir/$type.evm.out";
            }
            $ret = system "cat $evm_out_file >> $combo_dir/$type.evm.out";
            if ($ret) {
                die "Error archiving $evm_out_file ";
            }
        }
        
        ## get the EVM pred exon summary:
        my $exon_cmd = "$bindir/analyze_predicted_exons.pl $evaluationEntries evm.gff3 ";
        my $ret = system $exon_cmd;
        if ($ret) {
            die "Error, cmd $cmd died with ret $ret";
        }
        
        foreach my $file qw (exon_coords_to_type.summary exon_type_counts.summary genes_to_type.summary) {
            rename($file, "$combo_dir/$type.$file") or die "Error, could not rename $file";
        }
        
        print "Done with cmd: $cmd\n$result\n";
        
        ######################
        ## Run Eval:
        # -prepare inputs for eval:
        print "Running Eval for $type\n";
        
        my $prepare_cmd = "$bindir/prepare_Eval_inputs.pl --dir_listing $evaluationEntries --template_ev_type . --template_filename template.gff3 --genome_filename $genomicSeqFile --other_filename evm.gff3 --other_ev_type EVM --intergenic_flank $FLANK_REGIONS_SnSp_BP";
        $ret = system $prepare_cmd;
        if ($ret) {
            die "Error, cmd $cmd died with ret $ret";
        }
        
        # -now run eval
        $cmd = "evaluate_gtf.pl -g reference.GTF EVM.GTF > wu-eval.out";
        $ret = system $cmd;
        if ($ret) {
            die "Error, cmd $cmd died with ret $ret";
        }
        # relocate results:
        foreach my $file qw (reference.GTF EVM.GTF wu-eval.out) {
            rename($file, "$combo_dir/$type.$file") or die "Error, could not rename file $file";
        }
        
    }
    
    print "Finished.\n\n";
    return;
}



####
sub summarize_results {

    ## Summarize the results:
    my @structs;
    opendir (my $dir, $combo_dir) or die $!;
    while (my $file = readdir($dir)) {
        unless ($file =~ /^type\./) { # hidden file
            next;
        }
        my $struct = { type => $file };
        open (my $fh, "$combo_dir/$file") or die $!;
        while (<$fh>) {
            if (/^\#\#\#/) {
                chomp;
                my @x = split (/\s+/);
                my $param = $x[1];
                my $value = pop @x;
                $struct->{$param} = $value;
            }
        }
        close $fh;
        push (@structs, $struct);
    }
    
    @structs = reverse sort {$a->{'genes:'}<=>$b->{'genes:'}} @structs;
    foreach my $struct (@structs) {
        my $type = $struct->{type};
        my ($bjh, $fscore, $exons, $genes) = ($struct->{'BJH:'}, $struct->{'Fscore:'}, $struct->{'exons:'}, $struct->{'genes:'});
        
        print "$type\t$bjh\t$fscore\t$exons\t$genes\n";
    }
    return;
}


                        
                    
                        

