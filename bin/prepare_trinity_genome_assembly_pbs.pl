#!/usr/bin/env perl

=pod

=head1 USAGE

Prepare data for Trinity Genome Guided. Run it on the base directory where align_rnaseq_gsnap.pl has been run.

Mandatory Options:

	   -bam |files  :s  => BAM files, co-ordinate sorted. They will be processed as one big file (unless -split)
	OR -sam         :s  => SAM files, co-ordinate sorted. They will be processed separately (may or may not be a good idea)

Optional:

	-help
	-intron_max      :i => Maximum intron size
	-boundary        :i => Number of reads to define a boundary (def. 2). For deep RNA-Seq in very large scaffolds or genomic DNA contamination, 
you need to increase this (e.g 10 or 25; see distribution of bg/bigwig file) but you may lose lowly expressed transcripts.
        -minimum_reads   :i => Minimum number of reads required to process (defaults to 50)
        -small_cutoff    :i => Maximum file size of *.reads file to assign it as a 'small' and quick run (defaults to 1024^3, i.e. 1 megabyte)
        -medium_cutoff   :i => As above but for assigning it as 'medium', defaults to 10 x small_cutoff (higher than this is going to be assigned as 'long')
	-single_stranded :s => Give if single stranded library (if single end: F or R,  if paired:  FR or RF)
	-single_end         => Give this option if it is single end data.
	-cpu             :i => Number of CPUs (def 4)
	-memory          :s => Sorting memory to use, give as e.g. 20G (def 20G).
	-split_scaffolds    => Split alignments to one per reference sequence
        -scaffold_thread :i => Number of parallel processes to use with -split_scaffolds (def 3)
	-skip_sam	    => Jump into the processing regardless existing .reads regardless of the -bam / -sam option

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

#threaded
use Thread_helper;    # threads -1

my $cwd   = `pwd`;
chomp($cwd);
my $boundary = 2;
my $scaffold_size_cutoff = 3000;
my $intron_max_size = 70000;
my $minimum_reads = 50;
my $small_cut     = 1024 * 1024 * 1024;
my ($medium_cut,@user_provided_sam_files,@user_provided_bam_files,$delete_sam,@read_files,$is_single_stranded,$single_end,$help,$do_split_scaffolds);
my $debug;
my $cpus = 4;
my $memory = '20G';
my $threads                = 3;
my $skip_sam;
pod2usage $! unless &GetOptions(
	'memory:s'     => \$memory,
	'boundary:i'   => \$boundary,
	'cpu|thread:i' => \$cpus,
	'help'         => \$help,
        'debug'         => \$debug,
	'minimum_reads:i' => \$minimum_reads,
	'small_cutoff:i' => \$small_cut,
	'medium_cutoff:i' => \$medium_cut,
	'sam:s{,}'   => \@user_provided_sam_files,
	'bam|files_bam:s{,}' => \@user_provided_bam_files,
	'intron_max:i'   => $intron_max_size,
	'single_stranded:s' => \$is_single_stranded,
	'single_end' => \$single_end,
	'split_scaffolds' => \$do_split_scaffolds,
	'scaffold_threads:i'    => \$threads,
	'skip_sam'           => \$skip_sam
);
pod2usage if $help;
my ($samtools_exec,$TGG_prep_exec,$TGG_exec) = &check_program('samtools','prep_rnaseq_alignments_for_genome_assisted_assembly.pl','JAMG_TGG_cmds.pl');

if ($skip_sam){
	@read_files = `find . -maxdepth 6 -name "*.reads"`;
	chomp(@read_files);
}
############### PART ONE ###########################

if (!$read_files[0]){
	if (@user_provided_bam_files){
		foreach my $user_provided_bam_file (@user_provided_bam_files){
			die "Cannot find BAM file $user_provided_bam_file\n" unless -s $user_provided_bam_file;
		}
		print "Converting BAMs:\n".join(" ",@user_provided_bam_files)."\n";
		$delete_sam = 1;
		my $user_provided_sam_file = 'RNASeq_TGG_input.sam';
		system("rm -f $user_provided_sam_file.*.sam");
		if (scalar(@user_provided_bam_files > 1 )){
			&process_cmd("$samtools_exec merge - ".join(" ",@user_provided_bam_files)." |$samtools_exec view -@ $cpus -F4 - > $user_provided_sam_file " ) unless -s $user_provided_sam_file;
		}else{
			&process_cmd("$samtools_exec view -@ $cpus -F4 ".$user_provided_bam_files[0]." > $user_provided_sam_file") unless -s $user_provided_sam_file;
		}
		pod2usage ("Can't produce SAM file from input... Are they sorted by co-ordinate?\n") unless -s $user_provided_sam_file;
		push(@user_provided_sam_files,$user_provided_sam_file);
	}

	my @sam_files_to_process;

	if ($do_split_scaffolds){
		print "Splitting scaffolds...\n";
		foreach my $orig_sam (@user_provided_sam_files){
			push(@sam_files_to_process,&split_scaffold_sam($orig_sam));
		}
	}else{
		@sam_files_to_process = @user_provided_sam_files;
	}


	pod2usage ("Can't find SAM files\n") unless @sam_files_to_process && -s $sam_files_to_process[0];
        print "Will use SAM as input:\n".join(" ",@sam_files_to_process)."\n";

 	my $thread_helper = new Thread_helper($threads);
	foreach my $file (@sam_files_to_process){
		my $thread = threads->create( 'process_sam_file', $file );
		$thread_helper->add_thread($thread);
	}
	$thread_helper->wait_for_all_threads_to_complete();
}


############### PART TWO ###########################

@read_files = `find . -maxdepth 6 -name "*.reads"`;
chomp(@read_files);
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
system("$TGG_exec_cmd -reads_list_file stilltodo.list.k -threads 2 > small_trinity_GG.cmds") if -s "stilltodo.list.k";


# large readsets with Paired end data need to make use of jaccard
$TGG_exec_cmd .= " -jaccard_clip " if !$single_end;
system("$TGG_exec_cmd -reads_list_file stilltodo.list.m -memory 4G -threads 2 > medium_trinity_GG.cmds") if -s "stilltodo.list.m";

# huge datasets need normalizing
$TGG_exec_cmd .= " -normalize ";
system("$TGG_exec_cmd -reads_list_file stilltodo.list.g  -memory 4G -threads 6 > large_trinity_GG.cmds") if -s "stilltodo.list.g";

unlink("stilltodo.list");
unlink("stilltodo.list.k");
unlink("stilltodo.list.m");
unlink("stilltodo.list.g");

if ( -s "small_trinity_GG.cmds" ) {
 open( IN, "small_trinity_GG.cmds" );
 my @array = <IN>;
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
 my @array = <IN>;
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
 my @array = <IN>;
 close IN;
 open( OUT, ">large_trinity_GG.cmds." );
 print OUT shuffle(@array);
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
 return $ret;
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
	my $sam_file = shift;
	my $delete = shift;
	die "No file $sam_file\n" unless $sam_file && -s $sam_file;
	my %fh_hash;
	my @scaff_id_lines;
	my $previous_scaff_id;
	open (SAM,$sam_file);
	while (my $ln=<SAM>){
		my @data = split("\t",$ln);
		next unless $data[8];
		my $scaff_id = $data[2];
		my $sam_file = "$sam_file.$scaff_id.sam";
		my $fh;
		if (!$fh_hash{$scaff_id}){
			open (my $fh1,">>$sam_file") || die $!;
			$fh_hash{$scaff_id} = $fh1;
			$fh = $fh1;
			push(@scaff_id_lines,$scaff_id);
		}else{
			$fh = $fh_hash{$scaff_id};
		}
		print $fh $ln;
		if ($previous_scaff_id && $scaff_id ne $previous_scaff_id){
			my $old_fh =  $fh_hash{$previous_scaff_id};
			close $old_fh if $old_fh;
		}
		$previous_scaff_id = $scaff_id;
	}

	my @scaff_sams;
	foreach my $scaff_id (keys %fh_hash){
		my $fh = $fh_hash{$scaff_id};
		close $fh;
		my $sam_file = "$sam_file.$scaff_id.sam";
		unlink($sam_file) if -e $sam_file && !-s $sam_file;
		push(@scaff_sams,$sam_file) if -s $sam_file;
	}

	my $total_scaffolds = scalar(@scaff_id_lines);
	print "Splitted alignments per genome sequence above $scaffold_size_cutoff b.p. ($total_scaffolds scaffolds)\n";
	die "No SAM files produced from genome\n" unless @scaff_sams && -s $scaff_sams[0];
	unlink($sam_file) if $delete;
	return @scaff_sams;
}

sub process_sam_file(){
	my $local_sam_file = shift;
	return unless -s $local_sam_file;
	my $local_boundary = $boundary;
	$local_boundary *= 10 if ($local_boundary == 2 && -s $local_sam_file > 1e9);
	my $TGG_prep_cmd = "$TGG_prep_exec --coord_sorted_SAM  $local_sam_file -I $intron_max_size --sort_buffer $memory -C $local_boundary ";
	$TGG_prep_cmd .= " --SS_lib_type $is_single_stranded " if $is_single_stranded;
	$TGG_prep_cmd .= " --min_reads_per_partition $minimum_reads ";
	my $res = &process_cmd($TGG_prep_cmd);
	unlink($local_sam_file) if $delete_sam && (!$res || $res == 256);
}
