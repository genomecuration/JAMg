#!/usr/bin/env perl

=pod

=head1 USAGE

Prepare data for Trinity Genome Guided. Run it on the base directory where align_rnaseq_gsnap.pl has been run.

Mandatory Options:

	   -bam   :s  => BAM file, co-ordinate sorted
	OR -sam   :s  => SAM file, co-ordinate sorted
        OR -files :s  => Many BAM files, co-ordinate sorted

Optional:

	-intron_max    :i => Maximum intron size
        -minimum_reads :i => Minimum number of reads required to process (defaults to 50)
        -small_cutoff  :i => Maximum file size of *.reads file to assign it as a 'small' and quick run (defaults to 1024^3, i.e. 1 megabyte)
        -medium_cutoff :i => As above but for assigning it as 'medium', defaults to 10 x small_cutoff (higher than this is going to be assigned as 'long')

 Note: Except for -intron, the defaults should be ok for most projects. 

=cut


use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Time::localtime;
use List::Util 'shuffle';
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin";
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/trinityrnaseq/util";

my $cwd   = `pwd`;
chomp($cwd);
my $intron_max_size = 70000;
my $minimum_reads = 50;
my $small_cut     = 1024 * 1024 * 1024;
my ($medium_cut,$bam_file,$sam_file,@bam_files,$delete_sam,@read_files);
my $debug = 1;
&GetOptions(
        'debug'         => \$debug,
	'minimum_reads:i' => \$minimum_reads,
	'small_cutoff:i' => \$small_cut,
	'medium_cutoff:i' => \$medium_cut,
	'bam:s'   => \$bam_file,
	'sam:s'   => \$sam_file,
	'files_bam:s{,}' => \@bam_files,
	'intron_max:i'   => $intron_max_size
);

my ($samtools_exec,$TGG_prep_exec,$TGG_exec) = &check_program('samtools','prep_rnaseq_alignments_for_genome_assisted_assembly.pl','GG_write_trinity_cmds.pl');

@read_files = `find . -maxdepth 6 -name "*.reads"`;
chomp(@read_files);
if (!$read_files[0]){
	if ($sam_file && -s $sam_file){
		print "Will use $sam_file as input\n";
	}elsif ($bam_file && -s $bam_file){
		$delete_sam = 1;
		$sam_file = $bam_file.'.sam';
		&process_cmd("$samtools_exec view $bam_file > $sam_file");
	}elsif (@bam_files){
		$delete_sam = 1;
		pod2usage ("No input files...\n") unless scalar(@bam_files) > 1;
		$bam_file = 'RNASeq_TGG_input.bam';
		$sam_file = 'RNASeq_TGG_input.sam';
		&process_cmd("$samtools_exec merge $bam_file ".join(" ",@bam_files));
	}else{
		pod2usage ("No input files...\n");
	}
	pod2usage ("Can't produce SAM file from input... Are they sorted by co-ordinate?\n") unless $sam_file && -s $sam_file;
	print "Will use $sam_file as input to TGG\n";
	&process_cmd("$TGG_prep_exec --coord_sorted_SAM  $sam_file -I $intron_max_size");
	unlink($sam_file) if $delete_sam;
	@read_files = `find . -maxdepth 6 -name "*.reads"`;
	chomp(@read_files);
}
print "Now searching for .reads files within a depth of 6 subdirectories and producing new *_trinity_GG commands...\n";
die "No read files found.\n" unless @read_files && scalar(@read_files) > 1;

open( OUT,   ">stilltodo.list" );
open( SMALL, ">stilltodo.list.k" );
open( MED,   ">stilltodo.list.m" );
open( LARGE, ">stilltodo.list.g" );
$medium_cut    = $small_cut * 10 if !$medium_cut;
my $counter;

foreach my $file (@read_files) {
 my $base = $file;
 $base =~ s/.trinity.reads$//;
 my $sam = $base . '.sam';
 next unless -s $sam;
 unless ( -s $file . ".out.Trinity.fasta" ) {
  my $size = -s $file;
  next unless $size;
  if ( $size < $small_cut ) {
   my $reads = `wc -l < $file`;
   chomp($reads);
   $reads /= 2;
   if ( $reads < $minimum_reads ) {
    warn
"$file: Fewer than $minimum_reads reads. Skipping small readset ($reads reads; $size bytes)\n";
    next;
   }
   else {
    print SMALL $file . "\n";
   }
  }
  elsif ( $size >= $small_cut && $size < $medium_cut ) {
   print MED $file . "\n";
  }
  elsif ( $size >= $medium_cut ) {
   print LARGE $file . "\n";
  }
  print OUT $file . "\n";
  $counter++;
 }
}
close OUT;
close SMALL;
close MED;
close LARGE;

system( "rm -f small_trinity_GG.cmds* medium_trinity_GG.cmds* large_trinity_GG.cmds*");

system("$TGG_exec --reads_list_file stilltodo.list.k --paired  > small_trinity_GG.cmds") if -s "stilltodo.list.k";
system("$TGG_exec --jaccard_clip --reads_list_file stilltodo.list.m --paired  > medium_trinity_GG.cmds") if -s "stilltodo.list.m";
system("$TGG_exec --jaccard_clip --reads_list_file stilltodo.list.g --paired  > large_trinity_GG.cmds") if -s "stilltodo.list.g";

unlink("stilltodo.list");
unlink("stilltodo.list.k");
unlink("stilltodo.list.m");
unlink("stilltodo.list.g");

if ( -s "small_trinity_GG.cmds" ) {
 open( IN, "small_trinity_GG.cmds" );
 my @array;
 while ( my $ln = <IN> ) {
  $ln =~ s/JM 2G --CPU 4/JM 1G --CPU 1/;
  chomp($ln);
  $ln .= " >/dev/null\n";
  push( @array, $ln );
 }
 close IN;
 open( OUT, ">small_trinity_GG.cmds." );
 print OUT shuffle(@array);
 close OUT;
 rename( "small_trinity_GG.cmds.", "small_trinity_GG.cmds" );
 system("split -d -a 3 -l 1000 small_trinity_GG.cmds small_trinity_GG.cmds.");
 print "Produced small_trinity_GG.cmds for running in a cluster (e.g. with ParaFly)\n";
}

if ( -s "medium_trinity_GG.cmds" ) {
 open( IN, "medium_trinity_GG.cmds" );
 my @array;
 while ( my $ln = <IN> ) {
  $ln =~ s/JM 2G --CPU 4/JM 3G --CPU 2/;
  chomp($ln);
  $ln .= " >/dev/null\n";
  push( @array, $ln );
 }
 close IN;
 open( OUT, ">medium_trinity_GG.cmds." );
 print OUT shuffle(@array);
 close OUT;
 rename( "medium_trinity_GG.cmds.", "medium_trinity_GG.cmds" );
 system("split -d -a 3 -l 100 medium_trinity_GG.cmds medium_trinity_GG.cmds.");
 print "Produced medium_trinity_GG.cmds for running in a cluster (e.g. with ParaFly)\n";
}

if ( -s "large_trinity_GG.cmds" ) {
 open( IN,  "large_trinity_GG.cmds" );
 open( OUT, ">large_trinity_GG.cmds." );
 while ( my $ln = <IN> ) {
  $ln =~ s/JM 2G --CPU 4/JM 3G --CPU 6/;
  chomp($ln);
  $ln .= " >/dev/null\n";
  print OUT $ln;
 }
 close IN;
 close OUT;
 rename( "large_trinity_GG.cmds.", "large_trinity_GG.cmds" );
 system("split -d -a 3 -l 1 large_trinity_GG.cmds large_trinity_GG.cmds.");
 print "Produced large_trinity_GG.cmds for running in a cluster (e.g. with ParaFly)\n";
}



#######################################################################

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

