#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case bundling);
use FindBin;
use File::Basename;


my $usage = <<_EOUSAGE_;

###############################################################################################
#
#  --weights_file          :standard weights file indicating gene prediction tyes.
#
#  --gene_predictions      :base filename for the gff3 file containing gene predictions.
#  --genome                :base filename for the genome sequence fasta file
#  --evaluationEntries     :list of directories housing the relevant gene predictions files.
#
###############################################################################################


_EOUSAGE_

    ;


my ($weights_file, $gene_predictions, $evaluationEntries, $help, $genomicSeqFile);

&GetOptions ('weights_file=s' => \$weights_file,
             'gene_predictions=s' => \$gene_predictions,
             'evaluationEntries=s' => \$evaluationEntries,
             'h'=> \$help,
             'genome=s' => \$genomicSeqFile,
             );

if ($help || ! ($weights_file && $gene_predictions && $evaluationEntries)) {
    die $usage;
}


my $FLANK_REGIONS_SnSp_BP = 500; # flanking intergenic region sizes to include in evaluation.
my $bindir = $FindBin::Bin;


main: {
    
    ## prepare output directory
    my $baseEvalEntriesName = basename($evaluationEntries);
    my $outputdir = "Eval.$baseEvalEntriesName";
    if (-d $outputdir) {
        die "Error, intended output directory $outputdir already exists";
    }
    mkdir ($outputdir) or die "Error, cannot mkdir $outputdir";
    chmod (0777, $outputdir);
    
    # get list of prediction types:
    my @pred_types = &get_prediction_types($weights_file);
    
    foreach my $pred_type (@pred_types) {

        print "Processing prediction type: $pred_type\n";
        
        ## prepare the data types and run eval
        my $cmd = "$bindir/prepare_Eval_inputs.pl --dir_listing $evaluationEntries --template_ev_type . "
            . "--template_filename template.gff3 --genome_filename $genomicSeqFile "
            . "--other_filename $gene_predictions --other_ev_type $pred_type "
            . "--intergenic_flank $FLANK_REGIONS_SnSp_BP";
        
        my $ret = system ($cmd);
        if ($ret) {
            die "Error, cmd $cmd died with ret $ret";
        }

        ## run Eval:
        $cmd = "evaluate_gtf.pl -g reference.GTF $pred_type.GTF > wu-eval.$pred_type.out";
        $ret = system $cmd;
        if ($ret) {
            die "Error, cmd $cmd died with ret $ret";
        }

        ## rename the output files:
        foreach my $file ("reference.GTF", "$pred_type.GTF", "wu-eval.$pred_type.out") {
            rename ($file, "$outputdir/$file") or die "Error, could not rename $file to $outputdir/$file";
        }
        
        print "\nDone processing $pred_type\n\n";

    }

    print "Finished.\n";

    exit(0);
    
}
    
             
####
sub get_prediction_types {
    my $weights_file = shift;
    
    my @pred_types;
    
    open (my $fh, $weights_file) or die "Error, cannot open $weights_file";
    while (<$fh>) {
		if (/^\#/) { next; }
        chomp;
        my ($class, $type, $weight) = split (/\s+/);
        if ($class =~ /PREDICTION/) {
            push (@pred_types, $type);
        }
    }
    close $fh;
    
    return (@pred_types);
}


    
