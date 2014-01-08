#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use Getopt::Long;
use Data::Dumper;


$|++;

my $params = <<__PARAMS__;

######################################################################################
#
#  Required:
# 
#  --genome               genome sequence in fasta format
#  --weights              weights for evidence types file (note weight values are not extracted here!)
#  --gene_predictions       |  gene predictions gff3 file
#  --evaluationEntries             file containing the list of entries for evaulation
#
#  Optional but recommended:
#  --protein_alignments   protein alignments gff3 file
#  --transcript_alignments      transcript alignments gff3 file
#
#  --best_indiv_weights     file containing the ev_type(space)best_weight_value
#  --genefinder_weights     weights for ab-initio gene predictions (homology_prediction types are not included).
#
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

    $best_indiv_weights, $genefinder_weights, 
    
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
             
             "best_indiv_weights=s" => \$best_indiv_weights,
             "genefinder_weights=s" => \$genefinder_weights,
             
             );                                         

unless ($genomicSeqFile && $genePredictionsFile && $weightsFile && $evaluationEntries && $best_indiv_weights && $genefinder_weights) {
    die $params;
}

my $FLANK_REGIONS_SnSp_BP = 500; # default for regions flanking the reference gene structure to be used in SnSp calculations.



my @evaluation_entries_list = `cat $evaluationEntries`;
chomp @evaluation_entries_list;


my $combo_dir = "combos.BestIndiv";
unless (-d $combo_dir) {
    mkdir $combo_dir or die "Error, cannot mkdir $combo_dir";
    chmod 0777, $combo_dir;
}


my %ev_classes;
my $genefinder_weights_text = "";

my $bindir = $FindBin::Bin;

{ 
    
    my $core_cmd = "$bindir/train_weights.pl --genome $genomicSeqFile --gene_predictions $genePredictionsFile --training_entries $evaluationEntries "
        . " --just_score_existing_weights -w weights.test " ;
    if ($use_grid_computing) {
        $core_cmd .= " --use_grid_computing";
    }
    if ($NO_RECURSE_FLAG) {
        $core_cmd .= " --NO_RECURSE";
    }

    my %final_weights = &parse_weights_file($weightsFile);
    my %genefinder_weights = &parse_weights_file($genefinder_weights);
    my %best_individual_weights = &parse_weights_file($best_indiv_weights);
        
    {
        ## run just genefinders:
        my $process = { 
            type => "just_genefinders",
            cmd => $core_cmd,
        };
        &process_cmd($process, \%genefinder_weights);
    }
    
    {
        ## run all minus pasa:
        my %final_weights_minus_pasa;
        foreach my $ev_type (keys %final_weights) {
            unless ($ev_type =~ /alignassembly/i) {
                $final_weights_minus_pasa{$ev_type} = $final_weights{$ev_type};
            }
        }
        
        
        my $process = {
            type => "final_weights_minus_pasa",
            cmd => $core_cmd . " --transcript_alignments $transcriptAlignmentsFile --protein_alignments $proteinAlignmentsFile",
        };
        &process_cmd($process, \%final_weights_minus_pasa);
    }
        

    { ## all weights including pasa
        
        my $process = {
            type => "final_weights_plus_pasa",
            cmd => $core_cmd . " --transcript_alignments $transcriptAlignmentsFile "
                . " --protein_alignments $proteinAlignmentsFile"
                . " --terminalExons $terminalExonsFile "
                ,
        };
        &process_cmd($process, \%final_weights);
    }
    

    ## run through the best individual weights:
    
    foreach my $ev_type (keys %best_individual_weights) {
        
        my $ev_class = $ev_classes{$ev_type};
        
        my %weights_to_use = (%genefinder_weights, $ev_type => $best_individual_weights{$ev_type});
        
        my $evm_cmd = $core_cmd;
        
        if ($ev_class eq 'TRANSCRIPT') {
            $evm_cmd .= " --transcript_alignments $transcriptAlignmentsFile ";
            
            if ($ev_type =~ /alignAssembly/ && $terminalExonsFile) {
                $evm_cmd .= " --terminalExons $terminalExonsFile ";
            }
        }
        elsif ($ev_class eq 'PROTEIN') {
            $evm_cmd .= " --protein_alignments $proteinAlignmentsFile ";
        }
        
        
        my $process = { 
            type => "$ev_type.best_individual",
            cmd => $evm_cmd,
        };
        
        &process_cmd($process, \%weights_to_use);
        
    }
}

exit(0);

sub write_weights_file {
    my ($weights_href) = @_;
    
    open (my $fh, ">weights.test") or die "Error, cannot write weights";
    foreach my $ev_type (keys %$weights_href) {
        # write weights file.
        my $weight = $weights_href->{$ev_type};
        my $ev_class = $ev_classes{$ev_type};
        
        print $fh "$ev_class\t$ev_type\t$weight\n";
    }

    
    close $fh;
    return;
}


####
sub process_cmd {
    my $struct = shift;
    my $weights_href = shift;
    unless ($weights_href) {
        die "Error, need weights href";
    }
    &write_weights_file($weights_href);
    
    my $cmd = $struct->{cmd};
    my $type = $struct->{type};

    unless ($type) {
        die "Error, no type for struct: " . Dumper ($struct);
    }
    
    ## run it:
    print "\n\n** Running: $type\nCMD: $cmd\n";
    my $result = `$cmd`;
    if ($?) {
        die "Error, cmd $cmd died with ret ($?)\n";
    }
    
    open (my $fh,  ">$combo_dir/type.$type") or die $!;
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


####
sub parse_weights_file {
    my ($filename) = @_;
    

    my %weights;
    open (my $fh, $filename) or die $!;
    while (<$fh>) {
        unless (/\w/) { next; }
        if (/^\#/) {
            next; 
        }
        chomp;
        my ($ev_class, $ev_type, $weight) = split (/\s+/);
        $ev_classes{$ev_type} = $ev_class;
        
        unless ($weight =~ /\d/) { 
            die "Error, $_ weight value lacks number";
        }
        

        $weights{$ev_type} = $weight;
    }
    close $fh;
    
    return (%weights);
}
                        

