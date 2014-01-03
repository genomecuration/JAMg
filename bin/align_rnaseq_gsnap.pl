#!/usr/bin/env perl

=pod
=head1 NAME

align_rnaseq_gsnap.pl

=head1 USAGE

Run GSNAP on lots of RNASeq data.

Will automatically run against all files that match a pattern, eg if -pattern is '_1_' then "*_[12]_*" and also any files given in the command line.
Pairs are matched using the _[12]_ Files must only have one of _1_ or _2_ in their filename. The pattern can be changed with -pattern1 and -pattern2.

Mandatory:

 -fasta :s           FASTA of genome
 -dbname :s          Name of database for GMAP. Will create if it doesn't exist.
 
 
Optional:

 -input_dir :s       Directory with read files (defaults to current working directory)
 -intron_db :s       GMAP intron splice database
 -gmap_dir :s        Where the GMAP databases are meant to live (def. ~/databases/gmap)
 -intron_size :i     Maximum intron length (def. 70,000)
 -cpus :i            Number of CPUs/threads (def. 10)
 -help
 -pattern1            Pattern for automatching left pair files with *$pattern*.fastq (defaults to '_1_')
 -pattern2            Pattern for automatching right pair (defaults to '_2_')
 -nofail             Don't print out failures (I/O friendlyness if sparse hits expected). Otherwise captured as FASTQ
 -suffix             Build/use suffix array (fast, downweights SNPs, use for non-polymorphic genomes)
 -path_number        Maximum number of hits for the read pair. If more that these many hits, then nothing is returned (defaults to 50)

=head1 AUTHORS

 Alexie Papanicolaou

        CSIRO Ecosystem Sciences
        alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

Copyright 2012-2014 the Commonwealth Scientific and Industrial Research Organization. 
This software is released under the Mozilla Public License v.2.

It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.mozilla.org/MPL/2.0.

=cut

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;
use Time::localtime;
use File::Basename;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";

my ( $gmap_build_exec, $gsnap_exec, $samtools_exec ) =
  &check_program( "gmap_build", "gsnap", "samtools" );
&samtools_version_check($samtools_exec);
my ($input_dir,$pattern2,$debug,$genome,$genome_dbname,$nofails,$suffix,$help, $intron_splice_db);
my $cwd   = `pwd`;
chomp($cwd);
my $gmap_dir = $ENV{'HOME'}.'/databases/gmap/';
my $repeat_path_number = 50;
my $intron_length = 70000;
my $cpus          = 10;
my $memory        = '35G';
my $pattern = '_1_';
&GetOptions(
	     'debug'         => \$debug,
             'fasta:s'        => \$genome,
             'dbname:s'       => \$genome_dbname,
             'intron_db:s'    => \$intron_splice_db,
             'gmap_dir:s'     => \$gmap_dir,
             'intron_size:i'  => \$intron_length,
             'cpus|threads:i' => \$cpus,
             'memory:s'       => \$memory,
             'help'           => \$help,
             'pattern1:s'      => \$pattern,
             'pattern2:s'      => \$pattern2,
             'nofail'         => \$nofails,
             'suffix'        => \$suffix,
	     'input_dir:s'   => \$input_dir,
	     'path_number:i' => \$repeat_path_number,
);

pod2usage if $help;
pod2usage "No genome FASTA\n" unless $genome && -s $genome;
pod2usage "No GMAP genome database name\n" unless $genome_dbname;
pod2usage "GMAP database does not exist: $gmap_dir\n" unless -d $gmap_dir;

$input_dir = $cwd unless $input_dir;
my $samtools_sort_CPUs = int( $cpus / 2 ) > 2 ? int( $cpus / 2 ) : 2;
my $suff = "";
if ($memory =~s/([A-Z])$//){
 $suff = $1;
}

$memory = int(( $memory / $samtools_sort_CPUs ) ) < 1 ? '1G' :
  int(  ( $memory / $samtools_sort_CPUs ) )
  . $suff;    # samtools sort uses -memory per CPU

unless ($pattern2){
	$pattern2 = $pattern;
	$pattern2 =~ s/1/2/;
}

my @files = glob("$input_dir/*$pattern*");
push( @files, @ARGV );
my @verified_files;
for ( my $i = 0 ; $i < @files ; $i++ ) {
 next if $files[$i] =~/\.bz2$/ || $files[$i] =~/\.gz$/;
 if ( -s $files[$i] ) {
  push( @verified_files, $files[$i] );
 }
 else {
  warn "Skipping: did not find file " . $files[$i] . "\n";
 }
}
@files = @verified_files;
die "No files found!\n" unless @files;


my ($build_cmd,$align_cmd);
if ($suffix){
 $build_cmd = "$gmap_build_exec -D $gmap_dir -d $genome_dbname -T /tmp/\$USER -k 13 -b 10 -q 1 -e 0 $genome >/dev/null";
 $align_cmd ="$gsnap_exec -B 5 -D $gmap_dir -d $genome_dbname --nthreads=$cpus  --pairmax-rna=$intron_length  --localsplicedist=$intron_length -N 1 -Q --npaths=$repeat_path_number --format=sam --sam-use-0M --no-sam-headers ";
}else{
 $build_cmd = "$gmap_build_exec -D $gmap_dir -d $genome_dbname -T /tmp/\$USER -k 13 -b 10 -q 1 -e 0 --no-sarray $genome >/dev/null";
 $align_cmd ="$gsnap_exec --use-sarray=0 -B 5 -D $gmap_dir -d $genome_dbname --nthreads=$cpus  --pairmax-rna=$intron_length  --localsplicedist=$intron_length -N 1 -Q --npaths=$repeat_path_number --format=sam --sam-use-0M --no-sam-headers ";
}

$align_cmd .= " --nofails "            if $nofails;
$align_cmd .= " --fails-as-input "     if !$nofails;
$align_cmd .= " -s $intron_splice_db " if $intron_splice_db;

foreach my $file ( sort @files ) {
 my $pair = $file;
 $pair =~ s/$pattern/$pattern2/;
 next if $pair eq $file;
 unless ( -s $pair ) {
  warn "Didn't find pair of $file. Skipping\n";
  next;
 }
 my $base = basename($file);
 $base =~ s/$pattern.+//;
 my $group_id = $base;
 $base .= "_vs_$genome_dbname";
 if ( -s "gsnap.$base.log" ) {
  open( LOG, "gsnap.$base.log" );
  my @log = <LOG>;
  close LOG;
  next if $log[-1] && $log[-1] =~ /^GSNAP Completed/;
#  my @del = glob("gsnap.$base.*");
#  foreach (@del) { unlink($_); }
 }
 open( LOG, ">gsnap.$base.log" );
 &process_cmd($build_cmd) unless -d $gmap_dir . '/' . $genome_dbname;
 my $file_align_cmd = $align_cmd;
 $file_align_cmd .=
   " --split-output=gsnap.$base --read-group-id=$group_id $file $pair ";
 &process_cmd( $file_align_cmd, '.', "gsnap.$base*" )
   unless (    -s "gsnap.$base.concordant_uniq"
            || -s "gsnap.$base.concordant_uniq.bam" );
 unless ( -s "gsnap.$base.concordant_uniq.bam" ) {
  &process_cmd(
"$samtools_exec view -h -u -T $genome gsnap.$base.concordant_uniq | $samtools_exec sort -@ $samtools_sort_CPUs -l 9 -m $memory - gsnap.$base.concordant_uniq"
  );
  &process_cmd("$samtools_exec index gsnap.$base.concordant_uniq.bam");
  print LOG "\ngsnap.$base.concordant_uniq.bam:\n";
  &process_cmd(
    "$samtools_exec flagstat gsnap.$base.concordant_uniq.bam >> gsnap.$base.log"
  );
  unlink("gsnap.$base.concordant_uniq");
 }
 unless ( -s "gsnap.$base.concordant_mult.bam" ) {
  &process_cmd(
"$samtools_exec view -h -u -T $genome gsnap.$base.concordant_mult | $samtools_exec sort -@ $samtools_sort_CPUs -l 9 -m $memory - gsnap.$base.concordant_mult"
  );
  &process_cmd("$samtools_exec index gsnap.$base.concordant_mult.bam");
  print LOG "\ngsnap.$base.concordant_mult.bam:\n";
  &process_cmd(
    "$samtools_exec flagstat gsnap.$base.concordant_mult.bam >> gsnap.$base.log"
  );
  unlink("gsnap.$base.concordant_mult");
 }
 unless ( -s "gsnap.$base.concordant_uniq_mult.bam" ) {
  &process_cmd(
"$samtools_exec merge -@ $cpus  -l 9 gsnap.$base.concordant_uniq_mult.bam gsnap.$base.concordant_uniq.bam gsnap.$base.concordant_mult.bam"
  );
  &process_cmd("$samtools_exec index gsnap.$base.concordant_uniq_mult.bam");
  print LOG "\ngsnap.$base.concordant_uniq_mult.bam:\n";
  &process_cmd(
"$samtools_exec flagstat gsnap.$base.concordant_uniq_mult.bam >> gsnap.$base.log"
  );
 }

 print LOG "\nGSNAP Completed!\n";
 close LOG;
}

########################################
sub process_cmd {
 my ( $cmd, $dir, $delete_pattern ) = @_;
 print &mytime . "CMD: $cmd\n"     if $debug;
 chdir($dir)                       if $dir;
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  chdir($cwd) if $dir;
  &process_cmd_delete_fails( $dir, $delete_pattern ) if $delete_pattern;
  die "Error, cmd died with ret $ret\n";
 }
 chdir($cwd) if $dir;
 return;
}

sub mytime() {
 my @mabbr =
   qw(January February March April May June July August September October November December);
 my @wabbr = qw(Sunday Monday Tuesday Wednesday Thursday Friday Saturday);
 my $sec   = localtime->sec() < 10 ? '0' . localtime->sec() : localtime->sec();
 my $min   = localtime->min() < 10 ? '0' . localtime->min() : localtime->min();
 my $hour =
   localtime->hour() < 10 ? '0' . localtime->hour() : localtime->hour();
 my $wday = $wabbr[ localtime->wday ];
 my $mday = localtime->mday;
 my $mon  = $mabbr[ localtime->mon ];
 my $year = localtime->year() + 1900;
 return "$wday, $mon $mday, $year: $hour:$min:$sec\t";
}

###
sub check_program() {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  pod2usage "Error, path to a required program ($prog) cannot be found\n\n"
    unless $path =~ /^\//;
  chomp($path);
  $path = readlink($path) if -l $path;
  push( @paths, $path );
 }
 return @paths;
}

sub samtools_version_check() {
 my $samtools_exec = shift;
 my @version_lines = `$samtools_exec 2>&1`;
 foreach my $ln (@version_lines) {
  if ( $ln =~ /^Version:\s+\d+\.(\d+).(\d+)/i ) {
   die "Samtools version 1.19+ is needed\n" unless $1 >= 1 && $2 >= 19;

   #print "Good: Samtools version 1.19+ found\n";
  }
 }
}

sub process_cmd_delete_fails {
 my $dir     = shift;
 my $pattern = shift;
 my @delete  = glob( $dir . "/" . $pattern );
 foreach (@delete) { unlink $_; }
 die
"\nBreak requested while $dir/$pattern was being processed. Unfinished output deleted.\n";
}
