#!/usr/local/bin/perl

use lib ($ENV{EGC_SCRIPTS}, $ENV{EUK_MODULES});

use strict;
use warnings;
use CdbTools;
use Egc_library;

my $exon_bp = 10;
my $intron_bp = 20;
my $pictogram_length = $exon_bp + $intron_bp + 2;

my %data;
open (my $fh, "all_alt_splices.dat") or die $!;
while (<$fh>) {
    chomp;
    my ($type, $asmbl_id, $orient, 
        $pasa_acc_A, $end5_A, $end3_A, 
        $pasa_acc_B, $end5_B, $end3_B,
        $splice_length_diff) = split (/\t/);
    
    unless ($splice_length_diff <= 10) { next;}

    ## want 3'-most AG acceptors and 5'-most GT donors
    my $want_A_flag = 0;
    
    if ($orient eq "+") {
        
        if ($type eq "ACC") {
            
            # -----AG++++++  [A]
            # -------AG++++  [B]

            if ($end5_A > $end5_B) {
                $want_A_flag = 1;
            }
       
        } else {
            # donor
            
            # ++++GT-------------- [A]
            # ++++++++GT---------- [B]


            if ($end5_A < $end5_B) {
                $want_A_flag = 1;
            }
        }
    }
    else {
        ## minus strand

        if ($type eq 'ACC') {
            
            # ++++++++AG----------- [A]
            # ++++++++++++AG------- [B]

            if ($end5_A < $end5_B) {
                $want_A_flag = 1;
            }
            
        }
        else {
            # donor
            # --------------GT++++++ [A]
            # -----------GT+++++++++ [B]
            
            if ($end5_A > $end5_B) {
                $want_A_flag = 1;
            }
        }


    }

    if ($want_A_flag) {
        my $splice_A_key = "$end5_A,$end3_A,$orient";
        $data{$type}->{$asmbl_id}->{$splice_length_diff}->{$splice_A_key} = 1;
    }
    else {
        my $splice_B_key = "$end5_B,$end3_B,$orient";
        $data{$type}->{$asmbl_id}->{$splice_length_diff}->{$splice_B_key} = 1;
    }
}
close $fh;




foreach my $splice_type (keys %data) {
    
    my $asmbl_href = $data{$splice_type};
    foreach my $asmbl_id (keys %$asmbl_href) {

        my $genome_seq = cdbyank_linear ($asmbl_id, "osa1.genome");
                
        my $splice_length_diff_href = $asmbl_href->{$asmbl_id};
        
        foreach my $splice_length_diff (keys %$splice_length_diff_href) {
            open (my $fh, ">>$splice_type.delta_$splice_length_diff") or die $!;
            
            my $splice_coords_href = $splice_length_diff_href->{$splice_length_diff};

            
            foreach my $splice_coords_entry (keys %$splice_coords_href) {
                my ($lend, $rend, $orient) = split (/,/, $splice_coords_entry);
                
                my $subseq = "";
               
                if ($splice_type eq "ACC") {
                    if ($orient eq '+') {
                        my $position = $lend - $intron_bp;
                        $subseq = substr($genome_seq, $position -1, $pictogram_length);
                    }
                    else {
                        # minus strand
                        my $position = $lend - $exon_bp;
                        $subseq = substr($genome_seq, $position -1, $pictogram_length);
                        $subseq = &reverse_complement($subseq);
                    }
                    
                    ## make chars preceding exon lowercase
                    my @chars = split (//, uc $subseq);
                    for (my $i = 0; $i < $intron_bp+2; $i++) {
                        $chars[$i] = lc $chars[$i];
                    }
                    $subseq = join ("", @chars);
                    
                }
                else {
                    
                    ## donor site
                    if ($orient eq '+') {
                        my $position = $lend - $exon_bp;
                        $subseq = substr($genome_seq, $position -1, $pictogram_length);
                    }
                    else {
                        my $position = $lend - $intron_bp;
                        $subseq = substr($genome_seq, $position -1, $pictogram_length);
                        $subseq = &reverse_complement($subseq);
                    }
                    
                    ## make chars succeeding the exon lowercase
                    my @chars = split (//, uc $subseq);
                    for (my $i = $exon_bp; $i < $pictogram_length; $i++) {
                        $chars[$i] = lc $chars[$i];
                    }
                    $subseq = join ("", @chars);
                    
                }
                
                my $report = ">${asmbl_id}-$lend-$rend-[$orient,delta:$splice_length_diff]-$splice_type\n$subseq\n";
                print $report;
                print $fh $report;
            }
        }
    }
}





