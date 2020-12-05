#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");

my $cpus = 4;
my $sort_buffer = '1G';
my $tmpdir = $ENV{'TMP'};
$tmpdir = $ENV{'TMPDIR'} if !$tmpdir;
$tmpdir = '/tmp' if !$tmpdir;

my $repeat_gff_file = shift;
die "Give Repeat GFF file\n" unless $repeat_gff_file && -s $repeat_gff_file;

my $sort_exec = &check_sort_version;


my $repeat_hints_file=$repeat_gff_file.'.hints';

 open (IN,$repeat_gff_file);
 open (OUT,">$repeat_hints_file");
 while (my $ln=<IN>){
	next if $ln=~/^#/;
	my @data = split("\t",$ln);
	$data[2] = 'nonexonpart';
	$data[8] = 'src=RM;pri=6';
	print OUT join("\t",@data)."\n";
 }
 close OUT;
 close IN;

&process_cmd("$sort_exec -n -k 4,4 $repeat_hints_file| $sort_exec -s -n -k 5,5 | $sort_exec -s -n -k 3,3 | $sort_exec -s -k 1,1 -o $repeat_hints_file.");
rename("$repeat_hints_file.",$repeat_hints_file);

################3
sub check_sort_version(){
	my ($sort_exec) = &check_program('sort');
	my @v=`$sort_exec --version`;

	if ($v[0] && $v[0]=~/(\d+)\.(\d+)\s*$/){
		my $major = $1;
		my $minor = $2;
		if ($major >= 8 && $minor >= 6){
			return "$sort_exec -T $tmpdir --parallel $cpus -S $sort_buffer";
		}else{
			return "$sort_exec -S $sort_buffer -T $tmpdir";
		}
	}else{
		die "Sort of coreutils not found!";
	}
}

sub check_program() {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  pod2usage "Error, path to a required program ($prog) cannot be found\n\n"
    unless $path =~ /^\//;
  chomp($path);
  #$path = readlink($path) if -l $path;
  push( @paths, $path );
 }
 return @paths;
}
###
sub process_cmd {
 my ($cmd) = @_;
 print "CMD: $cmd\n";
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  die "Error, cmd died with ret $ret\n";
 }
 return $ret;
}
