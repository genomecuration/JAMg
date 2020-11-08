#!/usr/bin/perl

=pod
=head1 NAME

prepare_trinity_normalisation.pl

=head1 USAGE

Show command options for trinityrnaseq normalisation (v2.11)

Optional:

 -pattern1           Pattern for automatching left pair files with *'pattern1'*.fastq (defaults to '_1_')
 -pattern2           Pattern for automatching right pair (defaults to '_2_')
 -filetype       :s  Only process files ending with this text. Do NOT use a wildcard (e.g no *fastq, just fastq)
 -cpus           :i  Number of CPUs/threads, defaults to 20
 -memory             Memory, use suffix G M b (def '100G')
 -input_dir	     files are present in this directory (defaults to current directory)

=cut

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;


my $cwd = `pwd`;
chomp($cwd);
my $input_dir = $cwd;
my $memory = '100G';
my $cpus = 20;
my $help;

my ($pattern1,$pattern2,$filetype) = ('_1_','_2_','');

&GetOptions(
             'pattern1:s'      => \$pattern1,
             'pattern2:s'      => \$pattern2,
             'input_dir:s'     => \$input_dir,
             'memory:s'        => \$memory,
             'cpus|threads:i'  => \$cpus,
             'filetype:s'      => \$filetype,
             'help'            => \$help,


);
pod2usage if $help;

my @files1 = glob("$input_dir/*$pattern1*$filetype");
my @files2 = glob("$input_dir/*$pattern2*$filetype");

pod2usage "No data found\n" unless scalar(@files1)>0 || scalar(@files2)>0;

my $cmd = " --seqType fq --max_memory $memory --CPU $cpus ";
$cmd .= " --left " . join(",", @files1) if scalar(@files1)>0;
$cmd .= " --right " . join(",", @files2) if scalar(@files2)>0;
$cmd .= " --no_bowtie --normalize_max_read_cov 200 --workdir /dev/shm/".$ENV{'USER'}."/trinity --output trinity";

print "\n\n".$cmd."\n\n";
