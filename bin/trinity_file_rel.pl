#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;

#      --samples_file <string>         tab-delimited text file indicating biological replicate relationships.
#                                   ex.
#                                        cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
#                                        cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
#                                        cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
#                                        cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
#



=pod

=head1 USAGE

        -naming     :s  => Tell me (case-insensitive) if you data come from sra seqcenter1 or seqcenter2. Determines format of files:
                        SRA             "*_1.fastq*"    DEFAULT
                        seqcenter1      "*_1_sequence.fastq*"
                        seqcenter2      "*R1_*.fastq*"


=cut

my $naming;
my ($is_sra,$is_seqcenter1,$is_seqcenter2) = (1,0,0);


$is_sra = 1 if $naming && $naming=~/sra/i;
$is_seqcenter1 = 1 if $naming && $naming=~/seqcenter1/i;
$is_seqcenter2 = 1 if $naming && ($naming=~/seqcenter2/i || $naming=~/rama/i  ) ;

# change this to find all the files from left pairs
my @files;
@files = glob("*_1.fastq*trimmomatic") if $is_sra; #SRA
@files = glob("*_1_sequence.fastq*trimmomatic") if $is_seqcenter1;
@files = glob("*R1_*.fastq*trimmomatic") if $is_seqcenter2;


print "\nWill process these files:\n\t".join("\n\t",@files)."\n\n";

sleep(1);
open (OUT,">filelist");

foreach my $f (sort @files){
        my $pair = $f;
	my $basename = $f;
	if ($basename=~/(.+)\./){$basename = $1;}

        # change this to grab the pair's filename by substituting something
	if ($is_sra){
	        $pair=~s/_1.fastq/_2.fastq/;
		$basename=~s/_1.fastq//;
	}elsif ($is_seqcenter1){
	        $pair=~s/_1_sequence/_2_sequence/;
		$basename=~s/_1_sequence//;
	}elsif ($is_seqcenter2){
	        $pair=~s/R1_/R2_/ if $is_seqcenter2;
		$basename=~s/R1_//;
	}
	print OUT $basename."\t".$basename."\t".$f;
	print OUT "\t".$pair if $pair && -s $pair;
	print OUT "\n";
}


close OUT;
