#!/usr/bin/env perl

=pod
=head1 NAME

align_rnaseq_gsnap.pl

=head1 USAGE

Run GSNAP on lots of RNASeq data.

Will automatically run against all files that match a pattern, eg if -pattern is '_1_' then "*_[12]_*" and also any files given in the command line.
Pairs are matched using the _[12]_ Files must only have one of _1_ or _2_ in their filename. The pattern can be changed with -pattern1 and -pattern2.

Mandatory:

 -fasta          :s  FASTA of genome
 -dbname         :s  Name of database for GMAP. Will create if it doesn't exist.
 
 
Optional:

 -input_dir      :s  Directory with read files (defaults to current working directory)
 -gmap_dir       :s  Where the GMAP databases are meant to live (def. ~/databases/gmap)
 -cpus           :i  Number of CPUs/threads (def. 6). I don't recommend more than 6 in a system that has 12 CPUs
 -help
 -pattern1           Pattern for automatching left pair files with *'pattern1'*.fastq (defaults to '_1_')
 -pattern2           Pattern for automatching right pair (defaults to '_2_')
 -nofail             Don't print out failures (I/O friendlyness if sparse hits expected). Otherwise captured as FASTQ
 -suffix             Build/use suffix array (fast, downweights SNPs, use for non-polymorphic genomes)
 -path_number        Maximum number of hits for the read pair. If more that these many hits, then nothing is returned (defaults to 50)
 -commands_only  :s  Don't run commands, instead write them out into a file as specified by the option. Useful for preparing jobs for ParaFly
 -split_input    :i  Split the input FASTQ files to these many subfiles. Good for running large RNASeq datasets. Needs -commands_only above
 -notpaired          Data are single end. Don't look for pairs (use -pattern1 to glob files)
 -intron_db      :s  GMAP intron splice database
 -intron_size    :i  Maximum intron length (def. 70,000)
 -memory             Memory for samtools sorting, use suffix G M b (def '35G')
 -verbose

=head1 AUTHORS

 Alexie Papanicolaou

        CSIRO Ecosystem Sciences
        alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

Copyright 2012-2014 the Commonwealth Scientific and Industrial Research Organization. 
See LICENSE file for license info
It is provided "as is" without warranty of any kind.

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
my (
     $input_dir,               $pattern2,      $debug,
     $genome,                  $genome_dbname, $nofails,
     $suffix,                  $help,          $intron_splice_db,
     $just_write_out_commands, $split_input,   $notpaired
);
my $cwd = `pwd`;
chomp($cwd);
my $gmap_dir           = $ENV{'HOME'} . '/databases/gmap/';
my $repeat_path_number = 50;
my $intron_length      = 70000;
my $cpus               = 6;
my $memory             = '35G';
my $pattern1            = '_1_';

&GetOptions(
             'debug'           => \$debug,'verbose'=>\$verbose,
             'fasta:s'         => \$genome,
             'dbname:s'        => \$genome_dbname,
             'intron_db:s'     => \$intron_splice_db,
             'gmap_dir:s'      => \$gmap_dir,
             'intron_size:i'   => \$intron_length,
             'cpus|threads:i'  => \$cpus,
             'memory:s'        => \$memory,
             'help'            => \$help,
             'pattern1:s'      => \$pattern1,
             'pattern2:s'      => \$pattern2,
             'nofail'          => \$nofails,
             'suffix'          => \$suffix,
             'input_dir:s'     => \$input_dir,
             'path_number:i'   => \$repeat_path_number,
             'commands_only:s' => \$just_write_out_commands,
             'split_input:i'   => \$split_input,
             'notpaired'       => \$notpaired,
);

pod2usage if $help;
pod2usage "No genome FASTA\n" unless $genome && -s $genome;
pod2usage "No GMAP genome database name\n" unless $genome_dbname;
pod2usage "GMAP database does not exist: $gmap_dir\n" unless -d $gmap_dir;
pod2usage "Split input requires the -commands_only option\n"
  if $split_input && !$just_write_out_commands;
pod2usage "Split input cannot be more than 100\n"
  if $split_input && $split_input > 100;
pod2usage "You only provided pattern1 and didn't specify -notpaired\n" if $pattern1 && !$notpaired && !$pattern2;
$input_dir = $cwd unless $input_dir;
my $samtools_sort_CPUs = int( $cpus / 2 ) > 2 ? int( $cpus / 2 ) : 2;
my $suff = "";

if ( $memory =~ s/([A-Z])$// ) {
 $suff = $1;
}

$memory =
  int( ( $memory / $samtools_sort_CPUs ) ) < 1
  ? '1G'
  : int( ( $memory / $samtools_sort_CPUs ) )
  . $suff;    # samtools sort uses -memory per CPU

unless ( $pattern2 || $notpaired ) {
 $pattern2 = $pattern1;
 $pattern2 =~ s/1/2/;
}

my @files = glob("$input_dir/*$pattern1*");
push( @files, "$input_dir/$pattern1" ) if -s "$input_dir/$pattern1";
push( @files, @ARGV );
my %verified_files;
for ( my $i = 0 ; $i < @files ; $i++ ) {
 next if basename($files[$i])=~/^gsnap/;
 if ( -s $files[$i] ) {
  $verified_files{$files[$i]} = 1;
 }
 else {
  warn "Skipping: did not find file " . $files[$i] . "\n";
 }
}
@files = sort keys %verified_files;
die "No files found!\n" unless @files;
print "Found these files:\n".join("\n",@files)."\n";

my ( $build_cmd, $align_cmd );
if ($suffix) {
 $build_cmd =
"$gmap_build_exec -D $gmap_dir -d $genome_dbname -k 13 -q 1 -e 0 $genome >/dev/null";
 $align_cmd =
"$gsnap_exec -B 4 -D $gmap_dir -d $genome_dbname --nthreads=$cpus --localsplicedist=$intron_length -N 1 -Q --npaths=$repeat_path_number --format=sam --sam-use-0M --no-sam-headers ";
}
else {
 $build_cmd =
"$gmap_build_exec -D $gmap_dir -d $genome_dbname -k 13 -q 1 -e 0 --no-sarray $genome >/dev/null";
 $align_cmd =
"$gsnap_exec --use-sarray=0 -B 4 -D $gmap_dir -d $genome_dbname --nthreads=$cpus --localsplicedist=$intron_length -N 1 -Q --npaths=$repeat_path_number --format=sam --sam-use-0M --no-sam-headers ";
}

$align_cmd .= " --nofails "                    if $nofails;
$align_cmd .= " --fails-as-input "             if !$nofails;
$align_cmd .= " -s $intron_splice_db "         if $intron_splice_db;
$align_cmd .= " --pairmax-rna=$intron_length " if !$notpaired;

open( CMD, ">$just_write_out_commands" ) if $just_write_out_commands;

if ($notpaired) {
 my @files_to_do = &checked_unpaired_files(@files);
 &align_unpaired_files(@files_to_do);
}
else {
 my @files_to_do = &checked_paired_files(@files);
 &align_paired_files(@files_to_do);
}

close(CMD) if $just_write_out_commands;
########################################
sub process_cmd {
 my ( $cmd, $dir, $delete_pattern ) = @_;
 print &mytime . "CMD: $cmd\n" if $debug || $verbose;
 undef($dir) if $dir && $dir eq '.';
 if ($just_write_out_commands) {
  print CMD "cd $dir ; $cmd ; cd $cwd;  " if $dir;
  print CMD "$cmd; "                      if !$dir;
 }
 else {
  chdir($dir) if $dir;
  my $ret = system($cmd);
  if ( $ret && $ret != 256 ) {
   chdir($cwd) if $dir;
   &process_cmd_delete_fails( $dir, $delete_pattern ) if $delete_pattern;
   die "Error, cmd died with ret $ret\n";
   chdir($cwd) if $dir;
  }
 }
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

sub checked_unpaired_files() {
 my @files = @_;
 my @files_to_do;
 foreach my $file ( sort @files ) {
  if ($file =~ /\.bz2$/ || $file =~ /\.gz$/){
	push( @files_to_do, $file );
  }
  elsif ($split_input) {
   print "Splitting data for unpaired $file\n";
   my $lines = `wc -l < $file`;
   chomp($lines);
   my $number_of_lines = int( ( $lines / 4 ) / $split_input );
   $number_of_lines *= 4;
   die
"Number of lines is not as expected for FASTQ ($number_of_lines / $lines)\n"
     unless $number_of_lines % 4 == 0;
   system("split -a 3 -d -l $number_of_lines $file $file. ");
   my @new_files = glob("$file.0??");
   print "\t Adding ";

   foreach my $f (@new_files) {
    if ( $f =~ /$file\.\d+$/ ) {
     print " $f";
     push( @files_to_do, $f );
    }
   }
   print "\n";
  }
  else {
   push( @files_to_do, $file );
  }
 }
 return @files_to_do;
}

sub checked_paired_files() {
 my @files = @_;
 my @files_to_do;
 foreach my $file ( sort @files ) {
  my $pair = $file;
  $pair =~ s/$pattern1/$pattern2/;
  next if $pair eq $file;
  unless ( -s $pair ) {
   warn "Didn't find pair of $file. Skipping\n";
   next;
  }
  if ($file =~ /\.bz2$/ || $file =~ /\.gz$/){
	push( @files_to_do, $file );
  }
  elsif ($split_input) {
   print "Splitting data for pairs $file & $pair\n";
   my $lines = `wc -l < $file`;
   chomp($lines);
   my $number_of_lines = int( ( $lines / 4 ) / $split_input );
   $number_of_lines *= 4;
   die
"Number of lines is not as expected for FASTQ ($number_of_lines / $lines)\n"
     unless $number_of_lines % 4 == 0;
   system("split -a 3 -d -l $number_of_lines $file $file. ");
   system("split -a 3 -d -l $number_of_lines $file $pair. ");
   my @new_files = ( glob("$file.0??"), glob("$pair.0??") );
   print "\t Adding ";

   foreach my $f (@new_files) {
    if ( $f =~ /$file\.\d+$/ ) {
     print " $f";
     push( @files_to_do, $f );
    }
   }
   print "\n";
  }
  else {
   push( @files_to_do, $file );
  }
 }
 return @files_to_do;
}

sub align_unpaired_files() {
 my @files_to_do = @_;
 foreach my $file ( sort @files_to_do ) {
  my $base = basename($file);
  my $group_id;
  if ( $split_input && $base =~ /\.\d+$/ ) {
   $base =~ s/$pattern1.+(\.0\d\d)$/$1/;
   $group_id = $base;
   $group_id =~ s/\.0\d\d$//;
  }
  else {
   $base =~ s/$pattern1.+//;
   $group_id = $base;
  }
  print "Processing $group_id ($file)\n";
  $base .= "_vs_$genome_dbname";
  if ( -s "gsnap.$base.log" ) {
   open( LOG, "gsnap.$base.log" );
   my @log = <LOG>;
   close LOG;
   next if $log[-1] && $log[-1] =~ /^GSNAP Completed/;
  }
  open( LOG, ">gsnap.$base.log" );
  &process_cmd($build_cmd) unless -d $gmap_dir . '/' . $genome_dbname;
  my $base_out_filename = $notpaired ? "gsnap.$base.unpaired"  : "gsnap.$base.concordant";
  my $file_align_cmd = $align_cmd;

  $file_align_cmd .= ' --bunzip2 ' if $file =~ /\.bz2$/; 
  $file_align_cmd .= ' --gunzip ' if $file =~ /\.gz$/; 

  $file_align_cmd .=
    " --split-output=gsnap.$base --read-group-id=$group_id $file  ";
  &process_cmd( $file_align_cmd, '.', "gsnap.$base*" )
    unless (    -s "$base_out_filename"."_uniq"
             || -s "$base_out_filename"."_uniq.bam" );
  unless ( -s "$base_out_filename"."_uniq.bam" ) {
   &process_cmd(
"$samtools_exec view -h -u -T $genome $base_out_filename"."_uniq | $samtools_exec sort -@ $samtools_sort_CPUs -l 9 -m $memory - $base_out_filename"."_uniq"
   );
   &process_cmd("$samtools_exec index $base_out_filename"."_uniq.bam");
   print LOG "\n$base_out_filename"."_uniq.bam:\n";
   &process_cmd(
    "$samtools_exec flagstat $base_out_filename"."_uniq.bam >> gsnap.$base.log"
   );
   unlink("$base_out_filename"."_uniq");
  }
  unless ( -s "$base_out_filename"."_mult.bam" ) {
   &process_cmd(
"$samtools_exec view -h -u -T $genome $base_out_filename"."_mult | $samtools_exec sort -@ $samtools_sort_CPUs -l 9 -m $memory - $base_out_filename"."_mult"
   );
   &process_cmd("$samtools_exec index $base_out_filename"."_mult.bam");
   print LOG "\n$base_out_filename"."_mult.bam:\n";
   &process_cmd(
    "$samtools_exec flagstat $base_out_filename"."_mult.bam >> gsnap.$base.log"
   );
   unlink("$base_out_filename"."_mult");
  }
  unless ( -s "$base_out_filename"."_uniq_mult.bam" ) {
   &process_cmd(
"$samtools_exec merge -@ $cpus -l 9 $base_out_filename"."_uniq_mult.bam $base_out_filename"."_uniq.bam $base_out_filename"."_mult.bam"
   );
   &process_cmd("$samtools_exec index $base_out_filename"."_uniq_mult.bam");
   print LOG "\n$base_out_filename"."_uniq_mult.bam:\n";
   &process_cmd(
"$samtools_exec flagstat $base_out_filename"."_uniq_mult.bam >> gsnap.$base.log"
   );
  }

  print LOG "\nGSNAP Completed!\n";
  close LOG;
  print CMD " rm -f $file " if $just_write_out_commands && $split_input;
  print CMD "\n" if $just_write_out_commands;
 }

}

sub align_paired_files() {
 my @files_to_do = @_;
 foreach my $file ( sort @files_to_do ) {
  my $pair = $file;
  $pair =~ s/$pattern1/$pattern2/;
  my $base = basename($file);
  my $group_id;
  if ( $split_input && $base =~ /\.\d+$/ ) {
   $base =~ s/$pattern1.+(\.0\d\d)$/$1/;
   $group_id = $base;
   $group_id =~ s/\.0\d\d$//;
  }
  else {
   $base =~ s/$pattern1.+//;
   $group_id = $base;
  }
  print "Processing $group_id ($file)\n";
  $base .= "_vs_$genome_dbname";
  if ( -s "gsnap.$base.log" ) {
   open( LOG, "gsnap.$base.log" );
   my @log = <LOG>;
   close LOG;
   next if $log[-1] && $log[-1] =~ /^GSNAP Completed/;
  }
  open( LOG, ">gsnap.$base.log" );
  &process_cmd($build_cmd) unless -d $gmap_dir . '/' . $genome_dbname;
  my $base_out_filename = $notpaired ? "gsnap.$base.unpaired"  : "gsnap.$base.concordant";
  my $file_align_cmd = $align_cmd;

  $file_align_cmd .= ' --bunzip2 ' if $file =~ /\.bz2$/; 
  $file_align_cmd .= ' --gunzip ' if $file =~ /\.gz$/; 

  $file_align_cmd .=
    " --split-output=gsnap.$base --read-group-id=$group_id $file $pair ";
  &process_cmd( $file_align_cmd, '.', "gsnap.$base*" )
    unless (    -s "$base_out_filename"."_uniq"
             || -s "$base_out_filename"."_uniq.bam" );
  unless ( -s "$base_out_filename"."_uniq.bam" ) {
   &process_cmd(
"$samtools_exec view -h -u -T $genome $base_out_filename"."_uniq | $samtools_exec sort -@ $samtools_sort_CPUs -l 9 -m $memory - $base_out_filename"."_uniq"
   );
   &process_cmd("$samtools_exec index $base_out_filename"."_uniq.bam");
   print LOG "\n$base_out_filename"."_uniq.bam:\n";
   &process_cmd(
    "$samtools_exec flagstat $base_out_filename"."_uniq.bam >> gsnap.$base.log"
   );
   unlink("$base_out_filename"."_uniq");
  }
  unless ( -s "$base_out_filename"."_mult.bam" ) {
   &process_cmd(
"$samtools_exec view -h -u -T $genome $base_out_filename"."_mult | $samtools_exec sort -@ $samtools_sort_CPUs -l 9 -m $memory - $base_out_filename"."_mult"
   );
   &process_cmd("$samtools_exec index $base_out_filename"."_mult.bam");
   print LOG "\n$base_out_filename"."_mult.bam:\n";
   &process_cmd(
    "$samtools_exec flagstat $base_out_filename"."_mult.bam >> gsnap.$base.log"
   );
   unlink("$base_out_filename"."_mult");
  }
  unless ( -s "$base_out_filename"."_uniq_mult.bam" ) {
   &process_cmd(
"$samtools_exec merge -@ $cpus -l 9 $base_out_filename"."_uniq_mult.bam $base_out_filename"."_uniq.bam $base_out_filename"."_mult.bam"
   );
   &process_cmd("$samtools_exec index $base_out_filename"."_uniq_mult.bam");
   print LOG "\n$base_out_filename"."_uniq_mult.bam:\n";
   &process_cmd(
"$samtools_exec flagstat $base_out_filename"."_uniq_mult.bam >> gsnap.$base.log"
   );
  }

  print LOG "\nGSNAP Completed!\n";
  close LOG;
  print CMD " rm -f $file $pair " if $just_write_out_commands && $split_input;
  print CMD "\n" if $just_write_out_commands;
 }
}
