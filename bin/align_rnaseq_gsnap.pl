#!/usr/bin/env perl

=pod
=head1 NAME

align_rnaseq_gsnap.pl

=head1 USAGE

Run GSNAP on lots of RNASeq data. See source code

=cut

use strict;
use warnings;
use Pod::Usage;
use Time::localtime;
use File::Basename;

my $debug = 1;
my $cwd   = `pwd`;
chomp($cwd);
my $genome           = shift;
my $genome_dbname    = shift;
my $intron_splice_db = shift;

my ( $gmap_build_exec, $gsnap_exec ) = &check_program( "gmap_build", "gsnap" );

$genome =
'/archive/pap056/work/assembly_qc/harm_csiro4b/assembly_csiro4b.scaffolds.public.fsa'
  unless $genome;
$genome_dbname = 'harm_csiro4b_unmasked' unless $genome_dbname;

my $gmap_dir      = '/databases/gmap/';
my $intron_length = 70000;

my @files = glob("*_1_sequence*fastq");
die unless @files;

my $build_cmd =
"$gmap_build_exec -D $gmap_dir -d $genome_dbname -T /tmp/pap056 -k 13 -b 10 -q 1 -e 0 --no-sarray $genome >/dev/null";
my $align_cmd =
"$gsnap_exec --use-sarray=0 -B 5 -D $gmap_dir -d $genome_dbname --nthreads=10  --pairmax-rna=$intron_length  --localsplicedist=$intron_length -N 1 -Q --npaths=50 --format=sam --sam-use-0M --no-sam-headers --fails-as-input ";
$align_cmd .= " -s $intron_splice_db " if $intron_splice_db;

foreach my $file ( sort @files ) {
	my $pair = $file;
	$pair =~ s/_1_sequence/_2_sequence/;
	unless ( $pair ne $file && -s $pair ) {
		warn "Didn't find pair of $file. Skipping\n";
		next;
	}
	my $base = basename($file);
	$base =~ s/_1_sequence.+//;
	my $group_id = $base;
	$base .= "_vs_$genome_dbname";
	if ( -s "gsnap.$base.log" ) {
		open( LOG, "gsnap.$base.log" );
		my @log = <LOG>;
		close LOG;
		next if $log[-1] && $log[-1] =~ /^Done/;
		my @del = glob("gsnap.$base.*");
		foreach (@del) { unlink($_); }
	}
	open( LOG, ">gsnap.$base.log" );
	&process_cmd($build_cmd) unless -d $gmap_dir . $genome_dbname;
	my $file_align_cmd = $align_cmd;
	$file_align_cmd .=
	  " --split-output=gsnap.$base --read-group-id=$group_id $file $pair ";
	&process_cmd($file_align_cmd);
	close LOG;
	unless ( -s "gsnap.$base.concordant_uniq.bam" ) {
		&process_cmd(
"samtools view -u -T $genome gsnap.$base.concordant_uniq|samtools sort -m 30000000000 - gsnap.$base.concordant_uniq"
		);
		&process_cmd("samtools index gsnap.$base.concordant_uniq.bam");
		&process_cmd(
"samtools flagstat gsnap.$base.concordant_uniq.bam >> gsnap.$base.log"
		);
	}
}

########################################
sub process_cmd {
	my ( $cmd, $dir ) = @_;
	print &mytime . "CMD: $cmd\n"     if $debug;
	print LOG &mytime . "CMD: $cmd\n" if $debug;
	chdir($dir)                       if $dir;
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
	my @mabbr =
	  qw(January February March April May June July August September October November December);
	my @wabbr = qw(Sunday Monday Tuesday Wednesday Thursday Friday Saturday);
	my $sec = localtime->sec() < 10 ? '0' . localtime->sec() : localtime->sec();
	my $min = localtime->min() < 10 ? '0' . localtime->min() : localtime->min();
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
		pod2usage
		  "Error, path to a required program ($prog) cannot be found\n\n"
		  unless $path =~ /^\//;
		chomp($path);
		$path = readlink($path) if -l $path;
		push( @paths, $path );
	}
	return @paths;
}

