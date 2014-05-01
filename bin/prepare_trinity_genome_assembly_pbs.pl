#!/usr/bin/env perl

=pod

=head1 USAGE

Prepare data for Trinity Genome Guided. Run it on the base directory where align_rnaseq_gsnap.pl has been run.

Mandatory Options:

	   -bam   :s  => BAM file, co-ordinate sorted
	OR -sam   :s  => SAM file, co-ordinate sorted
        OR -files :s  => Many BAM files, co-ordinate sorted

Optional:

	-help
	-intron_max      :i => Maximum intron size
        -minimum_reads   :i => Minimum number of reads required to process (defaults to 50)
        -small_cutoff    :i => Maximum file size of *.reads file to assign it as a 'small' and quick run (defaults to 1024^3, i.e. 1 megabyte)
        -medium_cutoff   :i => As above but for assigning it as 'medium', defaults to 10 x small_cutoff (higher than this is going to be assigned as 'long')
	-single_stranded :s => Give if single stranded library (if single end: F or R,  if paired:  FR or RF)
	-single_end         => Give this option if it is single end data.
	-cpu             :i => Number of CPUs (def 4)
	-memory          :s => Sorting memory to use, give as e.g. 20G (def 20G).
	-split_scaffolds :s => Split alignments to one per reference sequence

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
$ENV{PATH} .= ":$RealBin:$RealBin/";
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin";
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/trinityrnaseq/util";
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/trinityrnaseq/util/support_scripts";

my $cwd   = `pwd`;
chomp($cwd);
my $intron_max_size = 70000;
my $minimum_reads = 50;
my $small_cut     = 1024 * 1024 * 1024;
my ($medium_cut,@sam_files,@bam_files,$delete_sam,@read_files,$is_single_stranded,$single_end,$help,$do_split_scaffolds);
my $debug = 1;
my $cpus = 4;
my $memory = '20G';
&GetOptions(
	'memory:s'     => \$memory,
	'cpu|thread:i' => \$cpus,
	'help'         => \$help,
        'debug'         => \$debug,
	'minimum_reads:i' => \$minimum_reads,
	'small_cutoff:i' => \$small_cut,
	'medium_cutoff:i' => \$medium_cut,
	'sam:s{,}'   => \@sam_files,
	'bam|files_bam:s{,}' => \@bam_files,
	'intron_max:i'   => $intron_max_size,
	'single_stranded:s' => \$is_single_stranded,
	'single_end' => \$single_end,
	'split_scaffolds' => \$do_split_scaffolds
);
pod2usage if $help;
my ($samtools_exec,$TGG_prep_exec,$TGG_exec) = &check_program('samtools','prep_rnaseq_alignments_for_genome_assisted_assembly.pl','JAMG_TGG_cmds.pl');

@read_files = `find . -maxdepth 6 -name "*.reads"`;
chomp(@read_files);

############### PART ONE ###########################

if (!$read_files[0]){
	if (@bam_files){
		foreach my $bam_file (@bam_files){
			die "Cannot find BAM file $bam_file\n" unless -s $bam_file;
		}
		print "Converting BAMs:\n".join(" ",@bam_files)."\n";
		$delete_sam = 1;
		my $sam_file = 'RNASeq_TGG_input.sam';
		if ($do_split_scaffolds){
			foreach my $bam_file (@bam_files){
				push(@sam_files,&split_scaffold_sam($bam_file));
			}
		}else{
			if (scalar(@bam_files > 1 )){
				&process_cmd("$samtools_exec merge - ".join(" ",@bam_files)." |$samtools_exec view -@ $cpus -F4 - > $sam_file " ) unless -s $sam_file;
			}else{
				&process_cmd("$samtools_exec view -@ $cpus -F4 ".$bam_files[0]." > $sam_file") unless -s $sam_file;
			}
			pod2usage ("Can't produce SAM file from input... Are they sorted by co-ordinate?\n") unless @sam_files && -s $sam_files[0];
			push(@sam_files,$sam_file);
		}
	}

	pod2usage ("Can't find SAM files\n") unless @sam_files && -s $sam_files[0];
        print "Will use SAM as input:\n".join(" ",@sam_files)."\n";

	foreach my $sam_file (@sam_files){
		my $TGG_prep_cmd = "$TGG_prep_exec --coord_sorted_SAM  $sam_file -I $intron_max_size --sort_buffer $memory --CPU $cpus ";
		$TGG_prep_cmd .= " --SS_lib_type $is_single_stranded " if $is_single_stranded;
		$TGG_prep_cmd .= " --min_reads_per_partition $minimum_reads ";
		&process_cmd($TGG_prep_cmd);
		unlink($sam_file) if $delete_sam;
		@read_files = `find . -maxdepth 6 -name "*.reads"`;
		chomp(@read_files);
	}
}

############### PART TWO ###########################

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
my $TGG_exec_cmd = "$TGG_exec";
$TGG_exec_cmd .= " -paired " if !$single_end;
system("$TGG_exec_cmd -reads_list_file stilltodo.list.k > small_trinity_GG.cmds") if -s "stilltodo.list.k";
# large readsets with Paired end data need to make use of jaccard
$TGG_exec_cmd .= " -jaccard_clip " if !$single_end;
system("$TGG_exec_cmd -reads_list_file stilltodo.list.m > medium_trinity_GG.cmds") if -s "stilltodo.list.m";
# huge datasets need normalizing
$TGG_exec_cmd .= " -normalize ";
system("$TGG_exec_cmd -reads_list_file stilltodo.list.g  > large_trinity_GG.cmds") if -s "stilltodo.list.g";

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


sub split_scaffold_sam(){
	my $bam_file = shift;
	return unless $bam_file && -s $bam_file;
	print "Splitting alignments per genome sequence\n";
	my @scaff_sams;
	my @scaff_id_lines;
	&process_cmd($samtools_exec." index $bam_file" ) unless -s "$bam_file.bai";	
	my @bam_scaff_ids = `$samtools_exec view -H $bam_file |grep '^@SQ'`);
	foreach my $ln (@bam_scaff_ids){
		if ($ln=~/\bSN:(\S+)\b/){	
			push(@scaff_id_lines,$1);
		}
	}
	foreach my $scaff_id (@scaff_id_lines){
		my $sam_file = "$bam_file.$scaff_id.sam";
		&process_cmd($samtools_exec."  view -@ $cpus -F4 $bam_file $scaff_id > $sam_file" ) unless -s $sam_file;
		unlink($sam_file) if -e $sam_file && !-s $sam_file;
		push(@scaff_sams,$sam_file) if -s $sam_file;
	}
	die "No SAM files produced from genome sequence $genome_sequence\n" unless @scaff_sams && -s $scaff_sams[0];
	return @scaff_sams;
}
