#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use CdbTools;
use Gene_obj_indexer;
use GFF3_utils;
use Getopt::Long qw(:config no_ignore_case bundling);

umask(0000);

my $usage = <<_EOUSAGE_;

####################################################################################################
#  
#  --template_accessions     the accessions (gene IDs) for those genes in the template gff3 file
#  --template_gff3_file      gold standard gene structures in gff3 format
#  --genome                  genome sequence in multi-fasta format
#  --evidence_gff3_files     comma-separated list of gff3 files containing evidence 
#                               ie.  gene_predictions.gff3,transcript_alignments.gff3,protein_alignments.gff3
#
#  --flanking_region         length of sequence (and assoc data) to flank the target gene [default: 30000]
#
# -v   verbose setting.
#
#####################################################################################################

_EOUSAGE_

    ;


my ($evidence_gff3_files, $genome_fasta_file, $gold_standard_list_file, $gold_standard_gff3_superset_file, $help, $VERBOSE);

my $OFFSET_INCLUDE = 30000; # 30kb on each end to include

&GetOptions ( 
              ## evm opts
              "template_accessions=s" => \$gold_standard_list_file,
              "template_gff3_file=s" => \$gold_standard_gff3_superset_file,
              "genome=s" => \$genome_fasta_file,
              "evidence_gff3_files=s" => \$evidence_gff3_files,
              "help|h" => \$help,
              "flanking_region=i" => \$OFFSET_INCLUDE,
              "v" => \$VERBOSE,
              
              );

if ($help) { die $usage; }

unless ($evidence_gff3_files && $genome_fasta_file && $gold_standard_list_file && $gold_standard_gff3_superset_file) {
    die $usage;
}


my @evidence_files = split (/,/, $evidence_gff3_files);
foreach my $ev_file (@evidence_files) { 
    $ev_file =~ s/\s//g;
}


my $util_dir = $FindBin::Bin;

## globals:
my %acc_to_genome_region;
my $fh;


## get gold standard entries
open ($fh, $gold_standard_list_file) or die "Error, cannot open $gold_standard_list_file";
while (<$fh>) {
    unless (/\w/) { next; }
    my @x = split (/\s+/);
    my $acc = shift @x;
    my $dirname = $acc;
    unless ($dirname) {
        die "Error, no accession retrieved from $gold_standard_list_file ";
    }

    $dirname =~ s/\W/_/g; 
    $acc_to_genome_region{$acc} = { dirname => $dirname }; #init to hashref
}
close $fh;

## init data repository
my $train_dir = "evm_train_dir";
unless (-d $train_dir) {
    mkdir ($train_dir) or die "Error, cannot mkdir $train_dir, $!";
}

## track entries by a train.entries file
my $train_listing = "all_train.entries";
open (my $trainfh, ">$train_listing") or die $!;
foreach my $acc (keys %acc_to_genome_region) {
    my $dirname = $acc_to_genome_region{$acc}->{dirname};
    my $dirpath = "$train_dir/$dirname";
    print $trainfh "$dirpath\n";
    unless (-d $dirpath) {
        mkdir ($dirpath) or die "Error, cannot mkdir $dirpath, $!";
        print "Creating directory: $dirpath\n" if $VERBOSE;
    }
}
close $trainfh;


## partition evidence data based on contig
print "-Parititioning all evidence based on contig identifier: --> ./evidence_partitioned/\*\n\n";

my %partitioned_files;
my $partitions_dir = "evidence_partitioned";
unless (-d $partitions_dir) {
    mkdir ($partitions_dir) or die "Error, cannot mkdir $partitions_dir $!";
    foreach my $evidence_file (@evidence_files) {
        print "-partitioning $evidence_file based on contig.\n";
        open (my $infh, $evidence_file) or die "Error, cannot open file $evidence_file";
        my $outfh;
        my $prev_contig = "";
        while (<$infh>) {
            if (/^\#/) { next; }
            unless (/\w/) { next;}
            my $line = $_;
            my @x = split (/\t/);
            my $contig_id = $x[0];
            if ($contig_id ne $prev_contig) {
                my $partitioned_file = "$partitions_dir/$evidence_file.$contig_id";
                $partitioned_files{$partitioned_file} = 1;
                if ($outfh) { 
                    close $outfh;
                }
                open ($outfh, ">>$partitioned_file") or die $!;
                print "\t-writing to $partitioned_file\n" if $VERBOSE;
            }
            print $outfh $line;
            $prev_contig = $contig_id;
        }
        close $infh;
    }
}

## index the gff3 files:
foreach my $gff3_file (keys %partitioned_files) {
    my $cmd = "$util_dir/iit_store -G -o $gff3_file.iit $gff3_file";
    my $ret = system $cmd;
    die "ERROR, cmd ($cmd) died with ret($ret)\n" if $ret;
    
    print "\n-Indexing $gff3_file, cmd: $cmd\n\n" if $VERBOSE;

}


## determine contig assignment and sequence ranges
my %contig_id_to_gene_list;
my $tmp_gene_obj_indexer = new Gene_obj_indexer( {create => "$gold_standard_gff3_superset_file.tmp.inx"} );
&GFF3_utils::index_GFF3_gene_objs($gold_standard_gff3_superset_file, $tmp_gene_obj_indexer);

my $gene_obj_indexer = new Gene_obj_indexer( { create => "$gold_standard_gff3_superset_file.inx" } );
foreach my $gene_key ($tmp_gene_obj_indexer->get_keys()) {
    my $gene_obj = $tmp_gene_obj_indexer->get_gene($gene_key);
    
    # store by the real gene ID.  (TODO: get rid of this!  Reconfigure the GFF3_utils)
    $gene_obj_indexer->store_gene($gene_obj->{TU_feat_name}, $gene_obj);
}

$tmp_gene_obj_indexer = undef;
unlink ("$gold_standard_gff3_superset_file.tmp.inx");


# group together those genes found on the same contig:
foreach my $acc (keys %acc_to_genome_region) {
    
    my $struct = $acc_to_genome_region{$acc};
    
    my $gene_obj = $gene_obj_indexer->get_gene($acc);
    my $asmbl_id = $gene_obj->{asmbl_id} or die "Error, no asmbl_id for " . $gene_obj->toString();
    my ($gene_lend, $gene_rend) = sort {$a<=>$b} $gene_obj->get_gene_span();

    $struct->{contig_id} = $asmbl_id;
    $struct->{lend} = $gene_lend;
    $struct->{rend} = $gene_rend;
    
    my $gene_list_aref = $contig_id_to_gene_list{$asmbl_id};
    unless (ref $gene_list_aref) {
        $gene_list_aref = $contig_id_to_gene_list{$asmbl_id} = [];
    }
    push (@$gene_list_aref, $acc);
        
}

## prepare genome sequence and gold standard structure, with adjusted coordinates:
# process via loop through contigs and corresponding gene lists:

foreach my $contig_id (keys %contig_id_to_gene_list) {
    my $gene_list_aref = $contig_id_to_gene_list{$contig_id};

    my $genome_seq = cdbyank_linear($contig_id, $genome_fasta_file, 'n');
    my $genome_seq_length = length($genome_seq);

    foreach my $gene_id (@$gene_list_aref) {
        print "processing template structure and genome seq: $contig_id, $gene_id\n";
        my $struct = $acc_to_genome_region{$gene_id};
        my $lend = $struct->{lend};
        my $rend = $struct->{rend};
        
        my $genome_start_pt = $lend - $OFFSET_INCLUDE;
        if ($genome_start_pt < 1) {
            $genome_start_pt = 1;
        }

        my $genome_end_pt = $rend + $OFFSET_INCLUDE;
        if ($genome_end_pt > $genome_seq_length) {
            $genome_end_pt = $genome_seq_length;
        }

        my $extract_length = $genome_end_pt - $genome_start_pt + 1;
        
        my $genome_subseq = substr($genome_seq, $genome_start_pt - 1, $extract_length);
        
        my $gene_obj = $gene_obj_indexer->get_gene($gene_id);
        
        my $gene_orientation = $gene_obj->get_orientation();

        $gene_obj->adjust_gene_coordinates(-1 * ($genome_start_pt - 1)); # want gene coords to map to genome substring
        
        # $gene_obj->create_all_sequence_types(\$genome_subseq);
        # print "$gene_id\t" . $gene_obj->get_CDS_sequence() . "\n\n";

        ## write the inputs:
        my $dirname = $struct->{dirname};
        my $dirpath = "$train_dir/$dirname";
        
        my $template_structure_file = "$dirpath/template.gff3";
        open (my $fh, ">$template_structure_file") or die $!;
        print $fh $gene_obj->to_GFF3_format();
        close $fh;
        print "\t-wrote $template_structure_file\n" if $VERBOSE;


        $genome_subseq =~ s/(\S{60})/$1\n/g;
        chomp $genome_subseq;
        open ($fh, ">$dirpath/$genome_fasta_file") or die $!;
        print $fh ">$contig_id [$genome_start_pt-$genome_end_pt]\n$genome_subseq\n";
        close $fh;
        print "\t-wrote $dirpath/$genome_fasta_file\n" if $VERBOSE;
        
        my $orient_token = ($gene_orientation eq '+') ? 'fwd' : 'rev';
        
        ## Gather the evidence within range: (there must be a faster way than this!!!)
      EVIDENCE_FILE:
        foreach my $evidence_file (@evidence_files) {
            print "\tprocessing $evidence_file, $gene_id\n";
            open (my $outfh, ">$dirpath/$evidence_file") or die $!; 
            my $evidence_partition_file = "$partitions_dir/$evidence_file.$contig_id";
            unless (-e $evidence_partition_file) {
                print STDERR "Warning, cannot find $partitions_dir/$evidence_file.$contig_id, maybe no evidence really exists here...\n";
                close $outfh;
                next EVIDENCE_FILE;
            }
            
            ## get the rows of the gff3 file (indexed):
            my @gff3_entries;
            my $cmd_template = "$util_dir/iit_get -A $evidence_partition_file.iit $genome_start_pt $genome_end_pt " . ${contig_id}; #${contig_id}_${orient_token} ";
            # retrieve hits from both strands within that range:
            foreach my $orient qw (+ -) {
                my $cmd = $cmd_template . "$orient";
                print "CMD: $cmd\n" if $VERBOSE;
                my @gff = `$cmd`;
                if ($?) {
                    die "Error, cmd: $cmd\ndied with ret $? ";
                }
                if (@gff) { 
                    push (@gff3_entries, @gff);
                }
            }
            
            foreach (@gff3_entries) {
                my @x = split (/\t/);
                my ($feat_lend, $feat_rend) = ($x[3], $x[4]);
                if ($feat_lend < $genome_end_pt && $feat_rend > $genome_start_pt) { # overlap
                    $feat_lend -= ($genome_start_pt - 1);
                    $feat_rend -= ($genome_start_pt - 1);
                    # keep coords w/in bounds:
                    if ($feat_lend < 1) {
                        $feat_lend = 1;
                    }
                    if ($feat_rend > $extract_length) {
                        $feat_rend = $extract_length;
                    }
                    
                    $x[3] = $feat_lend;
                    $x[4] = $feat_rend;
                    
                    print $outfh join ("\t", @x);
                }
            }
            close $outfh;
        }
    }
}


