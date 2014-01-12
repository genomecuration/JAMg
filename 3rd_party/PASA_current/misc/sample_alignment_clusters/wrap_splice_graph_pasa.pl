#!/usr/local/bin/perl

use FindBin;
use lib ($FindBin::Bin);
use Pasa_init;
use lib ("../PerlLib");
use strict;
use CDNA::Splice_graph_assembler;
use CDNA::CDNA_alignment;

open (STDERR, "&>STDOUT");

my $usage =  <<_EOH_;

usage: $0 alignment.textfile [archive_failures] 

The alignment.textfile should have the following format:

cdna_acc,transcribed_orient,segment_coords,...

The transcribed_orient should be '+|-|?'  
Use '?' in cases where you have single exon alignments of ambiguous transcribed orientation.

ie. 
// cluster: 3111
gi|86081823|gb|DR377582.1|DR377582,-,50113-50337,50419-50631,50883-50924
gi|86081670|gb|DR377427.1|DR377427,-,50109-50337,50419-50631,50883-50963,51126-51187
gi|86081223|gb|DR376980.1|DR376980,-,50109-50337,50419-50631,50883-50926
gi|86078045|gb|DR373802.1|DR373802,-,50113-50337,50419-50631,50883-50963,51064-51108
gi|86077457|gb|DR373214.1|DR373214,-,50113-50337,50419-50631,50883-50963,51126-51155
gi|86076745|gb|DR372502.1|DR372502,-,50177-50337,50419-50631,50883-50963,51126-51194
gi|86074966|gb|DR370723.1|DR370723,-,50444-50631,50883-50963,51126-51183
gi|86049609|gb|DR345364.1|DR345364,-,50236-50337,50419-50631,50883-50963,51126-51191
gi|86049607|gb|DR345362.1|DR345362,-,50920-50963,51126-51187
gi|86049606|gb|DR345361.1|DR345361,-,50920-50963,51126-51187    

// cluster ... etc, etc...

Namely,
accession,orientation,coordinates

_EOH_

    ;

my $input_file = $ARGV[0] or die $usage;
my $archive_failures_flag = $ARGV[1];


$/ = "\n//";

my $bindir = $FindBin::Bin;

my $failure_counter = 0;

open (my $fh, $input_file) or die "Error, cannot open $input_file";

my $tmpfile = "$ENV{HOSTNAME}.$$.tmp_cluster.txt";

while (my $input = <$fh>) {
    
    open (my $outfh, ">$tmpfile") or die $!;
    print $outfh $input;
    close $outfh;

    my $cmd = "$bindir/splice_graph_pasa.pl < $tmpfile ";
    my $ret = system $cmd;
    if ($ret) {
        if ($archive_failures_flag) {
            $failure_counter++;
            my $failure_dir = "$input_file.failures";
            unless (-d $failure_dir) {
                mkdir $failure_dir;
            }
            rename ("tmp_cluster.txt", "$failure_dir/failed_input.$$.$failure_counter.txt") or die "Error, cannot relocate failed cluster inputs";
        }
        else {
            
            die "Error, cmd ($cmd) died with ret ($ret)\n";
        }
    }
    
}

exit(0);
