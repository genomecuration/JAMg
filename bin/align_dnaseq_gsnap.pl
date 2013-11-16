#!/usr/bin/env perl

=pod
=head1 NAME

align_rnaseq_gsnap.pl

=head1 USAGE

Run GSNAP on lots of RNASeq data.

Will automatically run against all files that match "*_[12]_*fastq" and any files given in the command line.
Pairs are matched using the _[12]_ Files must only have one _1_ or _2_ pattern in their filename.

Mandatory:

 -fasta :s           FASTA of genome
 -dbname :s          Name of database for GMAP. Will create if it doesn't exist.
 
 
Optional:

 -gmap_dir :s        Where the GMAP databases are meant to live (def. /databases/gmap)
 -cpus :i            Number of CPUs/threads (def. 10)
 -help
 -pattern            Pattern for automatching left pair files with *$pattern*.fastq (defaults to _1_)
 -nofail             Don't print out failures (I/O friendlyness if sparse hits expected). Otherwise captured as FASTQ

=cut

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;
use Time::localtime;
use File::Basename;

my ( $gmap_build_exec, $gsnap_exec, $samtools_exec ) =
  &check_program( "gmap_build", "gsnap", "samtools" );
&samtools_version_check($samtools_exec);

my $debug = 1;
my $cwd   = `pwd`;
chomp($cwd);
my $gmap_dir = '/databases/gmap/';
my $genome;

# '/archive/pap056/work/assembly_qc/harm_csiro4b/assembly_csiro4b.scaffolds.public.fsa';
my $genome_dbname;

# 'harm_csiro4b_unmasked';
my $cpus   = 10;
my $memory = '35G';
my ($help);
my $pattern = '_1_';
my $nofails;

&GetOptions(
             'fasta:s'        => \$genome,
             'dbname:s'       => \$genome_dbname,
             'gmap_dir:s'     => \$gmap_dir,
             'cpus|threads:i' => \$cpus,
             'memory:s'       => \$memory,
             'help'           => \$help,
             'pattern:s'      => \$pattern,
             'nofail'         => \$nofails
);

pod2usage if $help;
pod2usage "No genome FASTA\n" unless $genome && -s $genome;
pod2usage "No GMAP genome database name\n" unless $genome_dbname;
pod2usage "GMAP database does not exist: $gmap_dir\n" unless -d $gmap_dir;

my $samtools_sort_CPUs = int( $cpus / 2 ) > 2 ? int( $cpus / 2 ) : 2;
my $suff = "";
if ($memory =~s/([A-Z])$//){
 $suff = $1;
}

$memory =
  sprintf( "%.2f", ( $memory / $samtools_sort_CPUs ) )
  . $suff;    # samtools sort uses -memory per CPU

my $pattern2 = $pattern;
$pattern2 =~ s/1/2/;

my @files = glob("*$pattern*.fastq");
push( @files, @ARGV );
my @verified_files;
for ( my $i = 0 ; $i < @files ; $i++ ) {
 if ( -s $files[$i] ) {
  push( @verified_files, $files[$i] );
 }
 else {
  warn "Skipping: did not find file " . $files[$i] . "\n";
 }
}
@files = @verified_files;
die "No files found!\n" unless @files;

my $build_cmd =
"$gmap_build_exec -D $gmap_dir -d $genome_dbname -T /tmp/pap056 -k 13 -b 10 -q 1 -e 0 $genome >/dev/null";
my $align_cmd =
"$gsnap_exec -B 5 -D $gmap_dir -d $genome_dbname --nthreads=$cpus -Q --npaths=50 --format=sam --sam-use-0M --no-sam-headers ";
$align_cmd .= " --nofails "        if $nofails;
$align_cmd .= " --fails-as-input " if !$nofails;

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
  my @del = glob("gsnap.$base.*");
  foreach (@del) { unlink($_); }
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
"$samtools_exec view -u -T $genome gsnap.$base.concordant_uniq | $samtools_exec sort -@ $samtools_sort_CPUs -l 9 -m $memory - gsnap.$base.concordant_uniq"
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
"$samtools_exec view -u -T $genome gsnap.$base.concordant_mult | $samtools_exec sort -@ $samtools_sort_CPUs -l 9 -m $memory - gsnap.$base.concordant_mult"
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
 print LOG &mytime . "CMD: $cmd\n" if $debug;
 chdir($dir)                       if $dir;
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  chdir($cwd) if $dir;
  &process_cmd_delete_fails( $dir, $delete_pattern ) if $delete_pattern;
  die "Error, cmd died with ret $ret\n";
 }
 chdir($cwd) if $dir;
 print LOG "Done\n";
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
