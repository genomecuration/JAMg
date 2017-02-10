#!/usr/bin/env perl

# example perl script to run it automatically.

=pod

=head1 NAME

 preprocess_illumina_automatic.pl

=head1 USAGE

	-naming	    :s	=> Tell me (case-insensitive) if you data come from sra seqcenter1 or seqcenter2. Determines format of files:
			SRA 		"*_1_fastq*"	DEFAULT
			seqcenter1 	"*_1_sequence.fastq*"
			seqcenter2	"*R1_*.fastq*"
	-parallel   :i	=> How many files to process in parallel. Warning, each file will use 4 CPUs. Defaults to 1
  
  Any other options you want to pass to preprocess_illumina.pl, end them also at the end (e.g. -noadapt)

=cut

use strict;
use warnings;
use FindBin qw/$RealBin/;
use lib ("$RealBin/PerlLib");
use Time::localtime;
use threads;
use Thread_helper;
use Pod::Usage;
use Getopt::Long qw/:config pass_through/;

$ENV{PATH} .= ":$RealBin";

my ($is_sra,$is_seqcenter1,$is_seqcenter2) = (1,0,0);
my $parallel               = 1;
my ($naming,$do_help);

&GetOptions(
	'name|naming:s' => \$naming,
	'parallel|do_parallel:i' => \$parallel,
	'help'	=> \$do_help,
);
my $extras = join(' ',@ARGV);

die pod2usage if $do_help;


$is_sra = 1 if $naming && $naming=~/sra/i;
$is_seqcenter1 = 1 if $naming && $naming=~/seqcenter1/i;
$is_seqcenter2 = 1 if $naming && ($naming=~/seqcenter2/i || $naming=~/rama/i  ) ;

# change this to find all the files from left pairs
my @files;
@files = glob("*_1_fastq*") if $is_sra; #SRA
@files = glob("*_1_sequence.fastq*") if $is_seqcenter1;
@files = glob("*R1_*.fastq*") if $is_seqcenter2;

my %check;
foreach my $f (@files){
	$check{$f}++ unless $f=~/trimmomatic/;
}
foreach my $f (@files){
	delete ($check{$f.".bz2"}) if $check{$f.".bz2"};
}

@files = sort (keys  %check);

print "\nWill process these files:\n\t".join("\n\t",@files)."\n\n";

sleep(5);

my $thread_helper = new Thread_helper($parallel);
my @failed_cmds;
foreach my $f (sort @files){
	my $cmd = "$RealBin/preprocess_illumina.pl ";
	my $pair = $f;
	# change this to grab the pair's filename by substituting something 
	$pair=~s/_1_fastq/_2_fastq/ if $is_sra; #SRA
	$pair=~s/_1_sequence/_2_sequence/ if $is_seqcenter1;
	$pair=~s/R1_/R2_/ if $is_seqcenter2;
	if (!-s $pair || $pair eq $f){
		$cmd .= " '$f'";
	}else{
		$cmd .= " -paired $extras '$f' '$pair'";
	}
	# can change the cmd to something you like
        my $thread        = threads->create('process_cmd',$cmd);
        $thread_helper->add_thread($thread);
        sleep(10);
}

$thread_helper->wait_for_all_threads_to_complete();
print "Done\n";
print "\nThese commands failed:\n\t".join("\n\t",@failed_cmds) if $failed_cmds[0];



####################################
sub process_cmd {
 my $cmd = shift;
 print &mytime . "CMD: $cmd\n";
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
   warn "Error, cmd died with ret $ret\n";
   push(@failed_cmds,$cmd);
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

