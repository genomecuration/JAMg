#!/usr/bin/env perl

use strict;
use warnings;

use FindBin qw/$RealBin/;
use Getopt::Long;

my $trinity_exec = "$RealBin/../3rd_party/trinityrnaseq/Trinity";

my $usage = <<__EOUSAGE__;

####################################################################################
#
# Required:
#
# --reads_list_file <string>      file containing list of filenames corresponding 
#                                  to the reads.fasta
#
# Optional:
#
# --paired                        reads are paired (default: not paired)
#
# --SS                            strand-specific  (reads are already oriented
#
# --jaccard_clip                  run jaccard clip
#
# --bfly_opts <string>            options to pass on to butterfly
#
#####################################################################################


__EOUSAGE__

    ;


my $reads_file;
my ($paired_flag,$SS_flag,$jaccard_clip,$normalize,$bfly_opts,$help_flag);

&GetOptions (
             'reads_list_file:s' => \$reads_file,
             'paired' => \$paired_flag,
             'SS' => \$SS_flag,
             'jaccard_clip' => \$jaccard_clip,
             'bfly_opts:s' => \$bfly_opts,
	     'normalize' => \$normalize,             
             'help' => \$help_flag,
             );
die $usage if ($help_flag);
die $usage unless ($reads_file && -s $reads_file);


open (my $fh, $reads_file) or die "Error, cannot open file $reads_file";
while (<$fh>) {
    chomp;
    my @x = split(/\s+/);
    
    my $file = pop @x;
    
    my $cmd = "$trinity_exec --seqType fa --single \"$file\" --JM 2G --CPU 4 --output \"$file.out\" --full_cleanup_ET ";
    $cmd .= " --run_as_paired " if $paired_flag;
    $cmd .= " --SS_lib_type F " if $$_flag;
    $cmd .= " --jaccard_clip " if $jaccard_clip;
    $cmd .= " --normalize_reads --normalize_max_read_cov 100 " if ($normalize);
    $cmd .= " --bfly_opts \"$bfly_opts\" " if $bfly_opts;
    print "$cmd\n";
}

exit(0);




		
