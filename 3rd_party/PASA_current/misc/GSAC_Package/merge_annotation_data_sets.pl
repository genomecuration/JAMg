#!/usr/local/bin/perl

use strict;
use warnings;
use lib $ENV{EUK_MODULES};
use Fasta_reader;

use Getopt::Long qw(:config no_ignore_case bundling);


my $usage = <<_EOUSAGE_;

###############################################################################
#
#  --comparison     the GSAC comparison output file
#  --protein_files   list of the protein files (ie. "file_A.pep,fileB.pep"
#  --top_hit_btab_files  the list of top hit btab files
#
###############################################################################


_EOUSAGE_

    ;

my ($comparison_file, $protein_files_list, $btab_files_list, $help);

&GetOptions( "comparison=s" => \$comparison_file,
             "protein_files=s" => \$protein_files_list,
             "top_hit_btab_files=s" => \$btab_files_list,
             "help|h" => \$help);


if ($help) { die $usage; }

unless ($comparison_file && $protein_files_list && $btab_files_list) {
    die $usage;
}

my @protein_fasta_files = split (/,/, $protein_files_list);
unless (scalar @protein_fasta_files == 2) { 
    die "Error, need exactly two protein fasta files: $protein_files_list";
}
my @btab_files = split (/,/, $btab_files_list);
unless (scalar @btab_files == 2) {
    die "Error, need exactly two btab files: $btab_files_list";
}


my %seqlengths; # length of protein sequence
my %alignment_lengths; # length of alignment to best database hit

open (my $logfh, ">comparison.$$.log") or die $!;


main: {
    foreach my $fasta_file (@protein_fasta_files) {
        &parse_sequence_lengths($fasta_file);
    }

    foreach my $btab_file (@btab_files) {
        &parse_btab_file($btab_file);
    }

    
    open (my $fh, $comparison_file) or die $!;
    while (<$fh>) {
        if (/^\#/) { next; }
        chomp;
        my $input_line = $_;
        
        print $logfh "// $input_line\n";

        my ($type, $contig, $model_list) = split (/\t/);
        my $models_selected = &resolve_conflicts($type, $model_list);
        
        unless ($models_selected) { die "Error, no models selected"; }
        
        print "$input_line\t$models_selected\n";
        
        print $logfh "*result: $models_selected\n\n";

    }
    close $fh;
}

exit(0);



####
sub resolve_conflicts {
    my ($type, $model_list) = @_;
    my %fully_qualified_name;
        
    my @models = split (/,/, $model_list);

    my %org_to_model_list;
    
    foreach my $model (@models) {
        my ($tu_feat, $model_feat, $org) = split (/:+/, $model);
        
        $fully_qualified_name{$model_feat} = $model;
        
        my $list_ref = $org_to_model_list{$org};
        unless (ref $list_ref) {
            $list_ref = $org_to_model_list{$org} = [];
        }

        push (@$list_ref, $model_feat);
    }
    
    my @org_types = keys %org_to_model_list;
    
    my @keep;
    if (scalar @org_types == 2) {
        my ($typeA, $typeB) = @org_types;
        my @modelsA = @{$org_to_model_list{$typeA}};
        my @modelsB = @{$org_to_model_list{$typeB}};
        
        @keep = &find_best_model_set(\@modelsA, \@modelsB);
    }
    elsif (scalar @org_types == 1) {
        @keep = @{$org_to_model_list{$org_types[0]}};
    }
    else {
        die "Error, do not have 1 or two sources for gene models: @org_types ";
    }
    
    my $keep_string = "";
    
    foreach my $model_feat (@keep) {
        my $full_name = $fully_qualified_name{$model_feat};
        $keep_string .= "$full_name,";
    }
    chop $keep_string; #strip trailing comma
    
    return ($keep_string);
    
}


####
sub find_best_model_set {
    my ($listA_aref, $listB_aref) = @_;
    
    my $listA_sum_align_length = &get_sum_align_length (@$listA_aref);
    my $listB_sum_align_length = &get_sum_align_length (@$listB_aref);
    print $logfh "alignment_coverage (@$listA_aref, sum:$listA_sum_align_length, @$listB_aref, sum:$listB_sum_align_length)\n";
    
    if ($listA_sum_align_length != $listB_sum_align_length) {
        
        if ($listA_sum_align_length > $listB_sum_align_length) {
            print $logfh "Case A1\n";
            return (@$listA_aref);
        }
        else {
            print $logfh "Case A2\n";
            return (@$listB_aref);
        }
        
    }
    
    else {
        ## Do CDS check
        my $listA_sum_prot_length = &get_sum_prot_lengths(@$listA_aref);
        my $listB_sum_prot_length = &get_sum_prot_lengths(@$listB_aref);
        
        print $logfh "prot_length_sums (@$listA_aref, sum:$listA_sum_prot_length, @$listB_aref, sum:$listB_sum_prot_length)\n";

        if ($listA_sum_prot_length >= $listB_sum_prot_length) {
            print $logfh "Case B1\n";
            return (@$listA_aref);
        }
        else {
            print $logfh "Case B2\n";
            return (@$listB_aref);
        }
    }
    
    
}


####
sub get_sum_align_length {
    my @feat_name_list = @_;
    my $sum = 0;
    foreach my $feat_name (@feat_name_list) {
        $sum += $alignment_lengths{$feat_name} || 0;
    }
    return ($sum);
}

####
sub get_sum_prot_lengths {
    my @feat_name_list = @_;
    my $sum = 0;
    foreach my $feat_name (@feat_name_list) {
        my $seqlen = $seqlengths{$feat_name} or die "Error, no sequence length for $feat_name\n";
        $sum += $seqlen;
    }

    return ($sum);
}


####
sub parse_sequence_lengths {
    my $fasta_file = shift;

    my $fasta_reader = new Fasta_reader($fasta_file);
    while (my $seq_obj = $fasta_reader->next()) {
        my $accession = $seq_obj->get_accession();
        my $sequence = $seq_obj->get_sequence();
        
        my $seq_length = length($sequence);
        $seqlengths{$accession} = $seq_length;
    }

    return;
}

####
sub parse_btab_file {
    my $btab_file = shift;

    open (my $fh, $btab_file) or die "Error, cannot open file $btab_file";
    while (<$fh>) {
        chomp;
        my @x = split (/\t/);
        my ($accession, $per_id, $end5, $end3) = ($x[0], $x[10], $x[8], $x[9]);
        
        if (defined($per_id) && $per_id > 1) {
            # got a hit
            my $alignment_length = abs ($end3 - $end5) + 1;
            $alignment_lengths{$accession} = $alignment_length;
        }
    }
    close $fh;

    return;

}
