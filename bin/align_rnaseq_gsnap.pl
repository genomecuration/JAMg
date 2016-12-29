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

 Input:
 -input_dir      :s  Directory with read files (defaults to current working directory)
 -gmap_dir       :s  Where the GMAP databases are meant to live (def. JAMG_PATH/databases/gmap)
 -pattern1           Pattern for automatching left pair files with *'pattern1'*.fastq (defaults to '_1_')
 -pattern2           Pattern for automatching right pair (defaults to '_2_')
 -piccard_0m         Ask gsnap to add 0M between insertions (only for piccard compatibility, issues with most other software)
 -filetype       :s  Only process files ending with this text. Do NOT use a wildcard (e.g no *fastq, just fastq)
 -split_input    :i  Split the input FASTQ files to these many subfiles (good for a single large readset). Needs -commands_only
 -commands_only  :s  Don't run commands, instead write them out into a file as specified by the option. Useful for preparing jobs for ParaFly
 -notpaired          Data are single end. Don't look for pairs (use -pattern1 to glob files)

 Speed:
 -cpus           :i  Number of CPUs/threads (def. 6). I don't recommend more than 6 in a system that has 12 CPUs
 -memory             Memory for samtools sorting, use suffix G M b (def '35G')
 -do_parallel    :i  Run these many alignments (if multiple input files) in parallel. Number of CPUs per alignment is -cpus divided by -do_parallel. Note, memory is not divided.
 -suffix             Build/use suffix array (fast, downweights SNPs, use for non-polymorphic genomes). Not suggested for RNAseq
 -build_only         Build genome (with suffix array) but don't do any alignments. Useful for building genome to be used many times
 -path_number        Maximum number of hits for the read pair. If more that these many hits, then nothing is returned (defaults to 10)
 -do_proportion  :i  Only process one sequence every this many reads (e.g. 1000). Good for doing a subset to build an intron DB.

 Output:
 -nofail             Don't print out failures (I/O friendlyness if sparse hits expected). Otherwise captured as FASTQ

 Introns:
 -intron_db      :s  GMAP intron or splice database (iit_store). Automatically found and used if it exists as a single *.iit in maps/ directory.
 -intron_hard        If -intron_db, conservative handling of junction reads and eliminate soft-clipping (--ambig-splice-noclip --trim-mismatch-score=0)
 -intron_size    :i  Maximum intron length (def. 70,000)

 Other:
 -verbose
 -help
 -large_genome       You have a very large genome (~ 4+ Gb). Use this option.

Note: For new genomes, it is possible to run this program once and use the output to build an -intron_db. Use the script augustus_RNAseq_hints.pl
 and pay attention to -min_jr_reads and maybe -intron_db_only. Make sure you use the same -genome FASTA. The output will be *junctions.all.hints.splice.gmap

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
use Carp;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;
use Time::localtime;
use File::Basename;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
use threads;
use Thread_helper;
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";

my ( $gmap_build_exec, $gsnap_exec, $samtools_exec,$bunzip2_exec,$bedtools_exec ) =
  &check_program( "gmap_build", "gsnap", "samtools",'bunzip2','bedtools' );
&samtools_version_check($samtools_exec);
my (
     $input_dir,               $pattern2,      $debug, $piccard_0m, $do_parallel,
     $genome,                  $genome_dbname, $nofails, $intron_hard,
     $suffix,                  $help,          $intron_splice_db,
     $just_write_out_commands, $split_input,   $notpaired, $verbose, $build_only
);
my $cwd = `pwd`;
chomp($cwd);
my $gmap_dir           = "$RealBin/../databases/gmap/";
my $repeat_path_number = 10;
my $intron_length      = 70000;
my $cpus               = 6;
my $memory             = '35G';
my $pattern1            = '_1_';
my $filetype = '';
my $do_proportion;
my $do_large_genome;
&GetOptions(
	     'do_parallel:i'   => \$do_parallel,
             'debug'           => \$debug,
	     'verbose'=>\$verbose,
             'fasta:s'         => \$genome,
             'dbname:s'        => \$genome_dbname,
             'intron_db:s'     => \$intron_splice_db,
             'intron_hard'     => \$intron_hard,
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
             'piccard_0m'      => \$piccard_0m,
	     'filetype:s'      => \$filetype,
             'build_only'      => \$build_only,
	     'do_proportion:i' => \$do_proportion,
	     'large_genome'    => \$do_large_genome,
);
pod2usage if $help;
pod2usage "No genome FASTA\n" unless $genome && -s $genome;
pod2usage "No GMAP genome database name\n" unless $genome_dbname;
pod2usage "GMAP database does not exist: $gmap_dir\n" unless -d $gmap_dir;
pod2usage "Split input requires the -commands_only option\n"
  if $split_input && !$just_write_out_commands;
pod2usage "Split input cannot be more than 100\n"
  if $split_input && $split_input > 100;
pod2usage "You only provided pattern1 or more than 1 input file but didn't specify -notpaired\n" if $pattern1 && !$notpaired && !$pattern2 && $ARGV[1];
$input_dir = $cwd unless $input_dir;
$do_parallel = 1 if $do_parallel && $do_parallel < 1;
if ($do_parallel && $do_parallel > 1){
	$cpus = int($cpus / $do_parallel);
}

( $gsnap_exec ) = &check_program( "gsnapl" ) if $do_large_genome;


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

my @files = glob("$input_dir/*$pattern1*$filetype");
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
die "No files found!\n" unless @files || $build_only;
print "Found these files:\n".join("\n",@files)."\n";

my ( $build_cmd, $align_cmd );
if ($suffix || $build_only) {
 warn "\n\nWARNING: Using -suffix with RNA-Seq gapped alignments is not recommended... it seems to be much (order of magnitude) slower in my tests\n\n\n" unless $build_only;
 $build_cmd = "$gmap_build_exec -D $gmap_dir -d $genome_dbname -e 0 $genome >/dev/null";
 $align_cmd =
"$gsnap_exec -B 5 -D $gmap_dir -d $genome_dbname --nthreads=$cpus --localsplicedist=$intron_length -N 1 -Q --npaths=$repeat_path_number --format=sam ";
}
else {
 $build_cmd = "$gmap_build_exec -D $gmap_dir -d $genome_dbname -e 0 --build-sarray=0 $genome >/dev/null";
 $align_cmd =
"$gsnap_exec --use-sarray=0 -B 5 -D $gmap_dir -d $genome_dbname --nthreads=$cpus  --localsplicedist=$intron_length -N 1 -Q --npaths=$repeat_path_number --format=sam ";
}

system($build_cmd) unless -d $gmap_dir . '/' . $genome_dbname;
system("$samtools_exec faidx $genome") unless -s "$genome.fai";
die "Failed to build genome ($genome.fai and $gmap_dir/$genome_dbname) " unless -s "$genome.fai" && -d "$gmap_dir/$genome_dbname";

if ($build_only){
	print "Build complete. User stop requested\n";
	exit(0);
}

if (!$intron_splice_db){
 my @checks = glob("$gmap_dir/$genome_dbname/$genome_dbname.maps/*iit");
 if (@checks && length(scalar(@checks)) == 1){
  $intron_splice_db = $checks[0];
  print "Will use the splice file $intron_splice_db\n";
  sleep(1);
 }
}


$align_cmd .= " --nofails "                    if $nofails;
$align_cmd .= " -s $intron_splice_db "         if $intron_splice_db;
$align_cmd .= " --ambig-splice-noclip --trim-mismatch-score=0 " if $intron_splice_db && $intron_hard;
$align_cmd .= " --pairmax-rna=$intron_length " if !$notpaired;
$align_cmd .= " --sam-use-0M " if $piccard_0m;
$align_cmd .= " --part=1/$do_proportion " if $do_proportion;

open( CMD, ">$just_write_out_commands" ) if $just_write_out_commands;

if ($notpaired) {
 my @files_to_do = &checked_unpaired_files(@files);
 if ($do_parallel && $do_parallel > 1){
	my $thread_helper = new Thread_helper($do_parallel);
	foreach my $f (@files_to_do){
	        my $thread        = threads->create('align_unpaired_files',($f));
	        $thread_helper->add_thread($thread);
		sleep(10);
	}
	$thread_helper->wait_for_all_threads_to_complete();
 }else{
	&align_unpaired_files(@files_to_do);
 }
}else {
 my @files_to_do = &checked_paired_files(@files);
 if ($do_parallel && $do_parallel > 1){
	my $thread_helper = new Thread_helper($do_parallel);
	foreach my $f (@files_to_do){
	        my $thread        = threads->create('align_paired_files',($f));
	        $thread_helper->add_thread($thread);
		sleep(10);
	}
	$thread_helper->wait_for_all_threads_to_complete();
 }else{
	&align_paired_files(@files_to_do);
 }
}

close(CMD) if $just_write_out_commands;
########################################
sub process_cmd {
 my ( $cmd, $dir, $delete_pattern ) = @_;
 print &mytime . "CMD: $cmd\n" if $debug || $verbose;
 undef($dir) if $dir && $dir eq '.';
 if ($just_write_out_commands) {
  print CMD "cd $dir && $cmd && cd $cwd;  " if $dir;
  print CMD "$cmd; "                      if !$dir;
 }
 else {
  chdir($dir) if $dir;
  my $ret = system($cmd);
  if ( $ret && $ret != 256 ) {
   chdir($cwd) if $dir;
   &process_cmd_delete_fails($cwd, $dir, $delete_pattern ) if $delete_pattern;
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
 my @progs = @_;
 my @paths;
 foreach my $prog (@progs) {
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
 my $cwd     = shift;
 my $dir     = shift;
 my $pattern = shift;
 $dir = $cwd if !$dir;
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
   die "Number of lines is not as expected for FASTQ ($number_of_lines / $lines)\n"
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
 my (@files_to_do) = @_;
 foreach my $file ( sort @files_to_do ) {
  my $qual_prot = &check_fastq_format($file);
  $qual_prot = $qual_prot eq 'fasta' ? '' : ' --quality-protocol='.$qual_prot;
  my $base = basename($file);
  my $group_id;
  if ( $split_input && $base =~ /\.\d+$/ ) {
   $base =~ s/$pattern1.+(\.0\d\d)$/$1/;
   if (!$base){$base = $pattern1; chop($base);}
   $group_id = $base;
   $group_id =~ s/\.0\d\d$//;
  }
  else {
   $base =~ s/$pattern1.+//;
   if (!$base){$base = $pattern1; chop($base);}
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
  my $base_out_filename = $notpaired ? "gsnap.$base.unpaired"  : "gsnap.$base.concordant";
  my $file_align_cmd = $align_cmd;

  $file_align_cmd .= ' --bunzip2 ' if $file =~ /\.bz2$/; 
  $file_align_cmd .= ' --gunzip ' if $file =~ /\.gz$/; 
  $file_align_cmd .= $qual_prot if $qual_prot;

  $file_align_cmd .= " --split-output=gsnap.$base --read-group-id=$base $file ";
  &process_cmd( $file_align_cmd, '.', "gsnap.$base*" )
    unless (    -s "$base_out_filename"."_uniq"
             || -s "$base_out_filename"."_uniq.bam" );
  unless ( -s "$base_out_filename"."_uniq.bam" || $just_write_out_commands) {
   &process_cmd("$samtools_exec view -h -u -T $genome $base_out_filename"."_uniq | $samtools_exec sort -@ $samtools_sort_CPUs -l 9 -m $memory -o $base_out_filename"."_uniq.bam -");
   &process_cmd("$samtools_exec index $base_out_filename"."_uniq.bam");

   ## For JBrowse
   &process_cmd("$bedtools_exec genomecov -split -bg -g $genome.fai -ibam $base_out_filename"
     ."_uniq.bam| sort -S 4G -k1,1 -k2,2n > $base_out_filename"."_uniq.coverage.bg");
   &process_cmd("bedGraphToBigWig $base_out_filename"."_uniq.coverage.bg $genome.fai $base_out_filename"."_uniq.coverage.bw") if `which bedGraphToBigWig`;

   print LOG "\n$base_out_filename"."_uniq.bam:\n";
   &process_cmd(
    "$samtools_exec flagstat $base_out_filename"."_uniq.bam >> gsnap.$base.log"
   );
   unlink("$base_out_filename"."_uniq");
  }
  unless ( -s "$base_out_filename"."_mult.bam" || $just_write_out_commands) {
   &process_cmd("$samtools_exec view -h -u -T $genome $base_out_filename"."_mult | $samtools_exec sort -@ $samtools_sort_CPUs -l 9 -m $memory -o $base_out_filename"."_mult.bam -");
   &process_cmd("$samtools_exec index $base_out_filename"."_mult.bam");
   print LOG "\n$base_out_filename"."_mult.bam:\n";
   &process_cmd(
    "$samtools_exec flagstat $base_out_filename"."_mult.bam >> gsnap.$base.log"
   );
   unlink("$base_out_filename"."_mult");
  }
# decided to remove as it was a resource hog
#  unless ( -s "$base_out_filename"."_uniq_mult.bam" ) {
#   &process_cmd("$samtools_exec merge -@ $cpus -l 9 $base_out_filename"."_uniq_mult.bam $base_out_filename"."_uniq.bam $base_out_filename"."_mult.bam"   );
#   &process_cmd("$samtools_exec index $base_out_filename"."_uniq_mult.bam");
#   print LOG "\n$base_out_filename"."_uniq_mult.bam:\n";
#   &process_cmd("$samtools_exec flagstat $base_out_filename"."_uniq_mult.bam >> gsnap.$base.log"   );
#  }

  print LOG "\nGSNAP Completed!\n" unless $just_write_out_commands;
  close LOG;
  print CMD " rm -f $file " if $just_write_out_commands && $split_input;
  print CMD "\n" if $just_write_out_commands;
 }

}

sub align_paired_files() {
 my @files_to_do = @_;
 foreach my $file ( sort @files_to_do ) {
  my $qual_prot = &check_fastq_format($file);
  $qual_prot = $qual_prot eq 'fasta' ? '' : ' --quality-protocol='.$qual_prot;
  my $pair = $file;
  $pair =~ s/$pattern1/$pattern2/;
  my $base = basename($file);
  my $group_id;
  if ( $split_input && $base =~ /\.\d+$/ ) {
   $base =~ s/$pattern1.+(\.0\d\d)$/$1/;
   if (!$base){$base = $pattern1; chop($base);}
   $group_id = $base;
   $group_id =~ s/\.0\d\d$//;
  }
  else {
   $base =~ s/$pattern1.+//;
   if (!$base){$base = $pattern1; chop($base);}
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
  my $base_out_filename = $notpaired ? "gsnap.$base.unpaired"  : "gsnap.$base.concordant";
  my $out_halfmapped = "gsnap.$base.halfmapping_uniq";

  my $file_align_cmd = $align_cmd;
  $file_align_cmd .= ' --bunzip2 ' if $file =~ /\.bz2$/; 
  $file_align_cmd .= ' --gunzip ' if $file =~ /\.gz$/; 
  $file_align_cmd .= $qual_prot if $qual_prot;

  $file_align_cmd .=    " --split-output=gsnap.$base --read-group-id=$base --allow-pe-name-mismatch $file $pair ";
  &process_cmd( $file_align_cmd, '.', "gsnap.$base*" )    unless (    -s "$base_out_filename"."_uniq" || -s "$base_out_filename"."_uniq.bam" );

  unless ( -s "$base_out_filename"."_uniq.bam" || $just_write_out_commands) {
   &process_cmd("$samtools_exec view -h -u -T $genome $base_out_filename"."_uniq | $samtools_exec sort -@ $samtools_sort_CPUs -l 9 -m $memory -o $base_out_filename"."_uniq.bam -");
   &process_cmd("$samtools_exec index $base_out_filename"."_uniq.bam");

   ## For JBrowse
   &process_cmd("$bedtools_exec genomecov -split -bg -g $genome.fai -ibam $base_out_filename"
     ."_uniq.bam| sort -S 4G -k1,1 -k2,2n > $base_out_filename"."_uniq.coverage.bg");
   &process_cmd("bedGraphToBigWig $base_out_filename"."_uniq.coverage.bg $genome.fai $base_out_filename"."_uniq.coverage.bw") if `which bedGraphToBigWig`;

   print LOG "\n$base_out_filename"."_uniq.bam:\n";
   &process_cmd( "$samtools_exec flagstat $base_out_filename"."_uniq.bam >> gsnap.$base.log" );
   unlink("$base_out_filename"."_uniq");
  }
  unless ( -s "$base_out_filename"."_mult.bam" || $just_write_out_commands) {
   &process_cmd("$samtools_exec view -h -u -T $genome $base_out_filename"."_mult | $samtools_exec sort -@ $samtools_sort_CPUs -l 9 -m $memory -o $base_out_filename"."_mult.bam -");
   &process_cmd("$samtools_exec index $base_out_filename"."_mult.bam");
   print LOG "\n$base_out_filename"."_mult.bam:\n";
   &process_cmd(
    "$samtools_exec flagstat $base_out_filename"."_mult.bam >> gsnap.$base.log"
   );
   unlink("$base_out_filename"."_mult");
  }
  if (!$just_write_out_commands && (-s $out_halfmapped && !-s "$out_halfmapped.bam")){
    &process_cmd("$samtools_exec view -h -u -T $genome $out_halfmapped | $samtools_exec sort -@ $samtools_sort_CPUs -l 9 -m $memory -o $out_halfmapped.bam -");
    &process_cmd("$samtools_exec index $out_halfmapped.bam");
    unlink($out_halfmapped) if -s "$out_halfmapped.bam";
  }
#  unless ( -s "$base_out_filename"."_uniq_mult.bam" ) {
#   &process_cmd("$samtools_exec merge -@ $cpus -l 9 $base_out_filename"."_uniq_mult.bam $base_out_filename"."_uniq.bam $base_out_filename"."_mult.bam"   );
#   &process_cmd("$samtools_exec index $base_out_filename"."_uniq_mult.bam");
#   print LOG "\n$base_out_filename"."_uniq_mult.bam:\n";
#   &process_cmd("$samtools_exec flagstat $base_out_filename"."_uniq_mult.bam >> gsnap.$base.log"   );
#  }

  print LOG "\nGSNAP Completed!\n" unless $just_write_out_commands;
  close LOG;
  print CMD " rm -f $file $pair " if $just_write_out_commands && $split_input;
  print CMD "\n" if $just_write_out_commands;
 }
}
sub check_fastq_format() {
 my $fq        = shift;
 my $max_seqs  = 100;
 my $max_lines = $max_seqs * 4;
 my ( @lines, $number, $counter );
 if ( $fq =~ /.bz2$/ ) {
  @lines = `$bunzip2_exec -dkc $fq|head -n $max_lines`;
 }
 else {
  @lines = `head -n $max_lines $fq`;
 }
 chomp(@lines);
 for ( my $k = 0 ; $k < @lines ; $k += 4 ) {
  my $ids = $lines[$k];
  if ($ids =~ /^>/){
        return 'fasta';
  }
  confess "$fq: Not in illumina format!\n" unless $ids =~ /^@/;
  $counter++;
  my $seq   = $lines[ $k + 1 ];
  my $idq   = $lines[ $k + 2 ];
  my $qual  = $lines[ $k + 3 ];
  my @quals = split( //, $qual );
  for ( my $i = 0 ; $i <= $#quals ; $i++ ) {
   $number = ord( $quals[$i] );
   if ( $number > 75 ) {
    warn "File $fq is solexa/illumina phred64 format!\n";
    return 'illumina';
   }
  }
  last if $counter >= $max_seqs;
 }
 return 'sanger';
}

