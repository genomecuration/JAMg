#!/usr/bin/env perl

=pod

=head1 USAGE

	-threads|CPU :i Number of CPUs to use
	-a3m         :s A3M Directory to process
        -cs219          Option to also produce cs219 files
        -base_hh     :s Base directory with hhblits/hhsuite software

=cut

use strict;
use warnings;
use File::Basename;
use threads;
use Pod::Usage;
use Thread_helper;
use Getopt::Long;

my $threads = 10;
my ($a3m_dir,$also_do_cs219,$base_hhblits);

GetOptions(
	'threads|CPU:i'=>\$threads,
	'a3m:s'  => \$a3m_dir,
	'cs219' => \$also_do_cs219,
	'base_hh:s' => \$base_hhblits
);

pod2usage "No A3M directory given or found\n" unless $a3m_dir && -d $a3m_dir;
pod2usage "CS219 processing requires the base directory of HHblits\n" if $also_do_cs219 && !$base_hhblits;
pod2usage "Cannot find the base directory of HHblits\n" if $base_hhblits && !-d $base_hhblits;


my $hhmake_exec = `which hhmake` || die("$!"); chomp($hhmake_exec);
my $hhm_dir = 'hhm_dir';
my ($cs219_dir,$cs_command);
mkdir ($hhm_dir) unless -d $hhm_dir;

if ($also_do_cs219){
	($cs219_dir,$cs_command) = ('cs219_dir',"cstranslate -D $base_hhblits/lib/hh/data/context_data.lib -A $base_hhblits/lib/hh/data/cs219.lib -x 0.3 -c 4");
	mkdir $cs219_dir unless -d $cs219_dir;
}

my $thread_helper = new Thread_helper($threads);


my @files = glob("$a3m_dir/*a3m");
my $a3m_files_num = scalar(@files);
my $hmm_files_num = int(0);
my $cs_files_num =int(0);
my $counter;
foreach my $file (@files){
	$counter++;
	if ($counter % 1000 == 0){
		print "\r Processed $counter / ".scalar(@files)."                   ";
	}
	my $base = $file;
	$base = basename($base);
	$base =~s/\.a3m$//;
	#cs219
	if ($also_do_cs219){
		my $thread = threads->create('do_cs219',$file,$base);
		$thread_helper->add_thread($thread);
	}
	#HHM
	my $outfile = "$hhm_dir/".$base.'.hhm';
	unlink($outfile) if -s $outfile;
	
	if (-s $file >= 10000){
		my $thread = threads->create('do_job',$file,$outfile);
		$thread_helper->add_thread($thread);
	}
}

 $thread_helper->wait_for_all_threads_to_complete();
 my @failed_threads = $thread_helper->get_failed_threads();
 if (@failed_threads) {
  die "Error, " . scalar(@failed_threads) . " threads failed.\n";
  exit(1);
 }

print "\nDone!\n";


###
sub do_job(){
 my ($in,$out) = @_;
 my $err = system($hhmake_exec." -i $in -o $out -v 0 ");
 unless ($err && $err > 0){
	$hmm_files_num++;
 }
}


sub do_cs219(){
	my ($in,$base) = @_;
	my $err = system($cs_command." -i $in -o $cs219_dir/$base.seq219  >/dev/null ");
	unless ($err && $err > 0){
		$cs_files_num++;
	}
}
