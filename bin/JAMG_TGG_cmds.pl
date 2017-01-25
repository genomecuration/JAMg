#!/usr/bin/env perl

=pod

=head1 USAGE 


Required:

 --reads_list_file <string>      file containing list of filenames corresponding 
                                  to the reads.fasta

Optional:

 --paired                        reads are paired (default: not paired)

 --SS                            strand-specific  (reads are already oriented

 --jaccard_clip                  run jaccard clip

 --bfly_opts <string>            options to pass on to butterfly

=cut

use strict;
use warnings;
use Pod::Usage;
use FindBin qw/$RealBin/;
use Getopt::Long;

my $trinity_exec = "$RealBin/../3rd_party/trinityrnaseq/Trinity";


my $reads_file;
my ($paired_flag,$SS_flag,$jaccard_clip,$normalize,$bfly_opts,$help_flag);
my $threads = 4;
my $inch_cpus = 4;
my $memory = '2G';


pod2usage $! unless &GetOptions (
             'reads_list_file:s' => \$reads_file,
             'paired' => \$paired_flag,
             'SS' => \$SS_flag,
             'jaccard_clip' => \$jaccard_clip,
             'bfly_opts:s' => \$bfly_opts,
	     'normalize' => \$normalize,
	     'threads:i' => \$threads,
	     'inch_cpus:i' => \$inch_cpus,
	     'memoryg:s' => \$memory,
             'help' => \$help_flag,
             );
pod2usage if ($help_flag);
pod2usage unless ($reads_file && -s $reads_file);


open (my $fh, $reads_file) or die "Error, cannot open file $reads_file";
while (<$fh>) {
    chomp;
    my @x = split(/\s+/);
    
    my $file = pop @x;
    
    my $cmd = "$trinity_exec --seqType fa --single \"$file\" --max_memory $memory --CPU $threads --inchworm_cpu $inch_cpus --no_bowtie --output \"$file.out\" --full_cleanup ";
    $cmd .= " --run_as_paired " if $paired_flag;
    $cmd .= " --SS_lib_type F " if $SS_flag;
    $cmd .= " --jaccard_clip " if $jaccard_clip;
    $cmd .= " --normalize_reads --normalize_max_read_cov 100 " if $normalize;
    $cmd .= " --bfly_opts \"$bfly_opts\" " if $bfly_opts;
    print "$cmd  >/dev/null\n";
}

exit(0);




		
