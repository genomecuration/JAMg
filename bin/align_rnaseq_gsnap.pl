#!/usr/bin/env perl


use strict;
use warnings;
use Time::localtime;
use File::Basename;

my $debug = 1;
my $cwd = `pwd`;chomp($cwd);
my $genome = shift;
my $genome_db = shift;
my $intron_splice_db = shift;
$genome = '/archive/pap056/work/assembly_qc/harm_csiro4b/assembly_csiro4b.scaffolds.public.fsa' unless $genome;
$genome_db = 'harm_csiro4b_unmasked' unless $genome_db;

my $gmap_dir = '/databases/gmap/';
my $intron_length = 70000;

my @files = glob("*_1_sequence*fastq");
die unless @files;

my $build_cmd = "gmap_build -D $gmap_dir -d $genome_db -T /tmp/pap056 -k 13 -b 10 -q 1 -e 0 --no-sarray $genome >/dev/null";
my $align_cmd = "gsnap --use-sarray=0 -B 5 -D $gmap_dir -d $genome_db --nthreads=10  --pairmax-rna=$intron_length  --localsplicedist=$intron_length -N 1 -Q --npaths=50 --format=sam --sam-use-0M --no-sam-headers --fails-as-input ";
$align_cmd .= " -s $intron_splice_db " if $intron_splice_db;

foreach my $file (sort @files){
	my $pair = $file;
	$pair =~s/_1_sequence/_2_sequence/;
	unless ($pair ne $file && -s $pair){
		warn "Didn't find pair of $file. Skipping\n";
		next;
	}
	my $base = basename($file);
	$base =~s/_1_sequence.+//;
        my $group_id = $base;
	$base.="_vs_$genome_db";
        if (-s "gsnap.$base.log"){
		open (LOG,"gsnap.$base.log");
		my @log = <LOG>;
		close LOG;
		next if $log[-1] && $log[-1]=~/^Done/;
		my @del = glob("gsnap.$base.*");
		foreach (@del){unlink($_);}
	}
	open (LOG,">gsnap.$base.log");
	&process_cmd($build_cmd) unless -d $gmap_dir.$genome_db;
	my $file_align_cmd = $align_cmd;
	$file_align_cmd.=" --split-output=gsnap.$base --read-group-id=$group_id $file $pair ";
	&process_cmd($file_align_cmd);
	close LOG;
}

########################################
sub process_cmd {
  my ( $cmd, $dir ) = @_;
  print &mytime."CMD: $cmd\n" if $debug;
  print LOG &mytime."CMD: $cmd\n" if $debug;
  chdir($dir) if $dir;
  my $ret = system($cmd);
  if ( $ret && $ret != 256 ) {
    chdir($cwd) if $dir;
    die "Error, cmd died with ret $ret\n";
  }
  chdir($cwd) if $dir;
  print LOG "Done\n";
  return;
}

sub mytime() {
  my @mabbr = qw(January February March April May June July August September October November December);
  my @wabbr = qw(Sunday Monday Tuesday Wednesday Thursday Friday Saturday);
  my $sec = localtime->sec() < 10 ? '0' . localtime->sec() : localtime->sec();
  my $min = localtime->min() < 10 ? '0' . localtime->min() : localtime->min();
  my $hour = localtime->hour() < 10 ? '0' . localtime->hour() : localtime->hour();
  my $wday = $wabbr[localtime->wday];
  my $mday = localtime->mday;
  my $mon = $mabbr[localtime->mon];
  my $year = localtime->year() + 1900;
  return "$wday, $mon $mday, $year: $hour:$min:$sec\t";
}

