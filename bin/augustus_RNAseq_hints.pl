#!/usr/bin/env perl

=pod

=head1 NAME

 augustus_RNAseq_hints.pl

=head1 USAGE

Create hint files for Augustus using RNASeq/EST. One is junction reads (excellent for introns), the other is RNASeq/EST coverage

 -dir|aug        The directory where Augustus is (if augustus is not already in your PATH).
 -bam|in         The input BAM file (co-ordinate sorted).
 -genome|fasta   The genome assembly FASTA file.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin";

my ($samtools_exec,$bedtools_exec) = &check_program('samtools','bedtools');
my ($augustus_exec);
#Options
my ($augustus_dir,$bamfile,$genome,$help);

GetOptions(
    'help'   => \$help,
    'dir|aug:s'  => \$augustus_dir,
    'bam|in:s'  => \$bamfile,
    'genome|fasta:s' => \$genome,
);

pod2usage if $help;
if ($augustus_dir){
	$augustus_dir = readlink($augustus_dir) if -l $augustus_dir;
	$augustus_exec = $augustus_dir.'/bin/augustus';
}else{
  $augustus_exec = &check_program('augustus');
  $augustus_dir = dirname(dirname($augustus_exec)) if $augustus_exec && !$augustus_dir;
}

pod2usage "Cannot find Augustus\n" unless $augustus_exec && -s $augustus_exec;

$ENV{PATH} .= ":$augustus_dir/bin:$augustus_dir/scripts";

pod2usage "Cannot find the BAM or genome FASTA file\n" unless $bamfile && -s $bamfile && $genome && (-s $genome || -s $genome.'.fai');

my $genome_idx_cmd = "$samtools_exec $genome";
&process_cmd($genome_idx_cmd) unless -s $genome.'.fai';
die "Cannot index genome $genome\n" unless -s $genome.'.fai';

my $junction_cmd = "$samtools_exec rmdup -S $bamfile - | $bedtools_exec bamtobed -bed12 | bed12_to_augustus_junction_hints.pl | $augustus_dir/scripts/join_mult_hints.pl > $bamfile.junctions.hints";
&process_cmd($junction_cmd) unless -s "$bamfile.junctions.hints";

my $coverage_cmd = "$bedtools_exec genomecov -split -bg -g $genome.fai -ibam $bamfile >  $bamfile.coverage.bg";
my $hint_cmd = "$augustus_dir/scripts/bedgraph2wig.pl --bedgraphfile=$bamfile.coverage.bg | $augustus_dir/scripts/wig2hints.pl | $augustus_dir/scripts/join_mult_hints.pl >  $bamfile.coverage.hints";
&process_cmd($coverage_cmd) unless -s "$bamfile.coverage.bg";
&process_cmd($hint_cmd) unless -s "$bamfile.coverage.hints";

if (-s "$bamfile.junctions.hints" && -s "$bamfile.junctions.hints"){
 print "Done!\n";
}else{
 die "Something went wrong....\n";
}
###
sub check_program(){
	my @paths;
	foreach my $prog (@_){
	        my $path = `which $prog`;
        	pod2usage "Error, path to a required program ($prog) cannot be found\n\n" unless $path =~ /^\//;
		chomp($path);
		$path = readlink($path) if -l $path;
		push(@paths,$path);
	}
	return @paths;
}
###
sub process_cmd {
  my ( $cmd) = @_;
  print "CMD: $cmd\n";
  my $ret = system($cmd);
  if ( $ret && $ret != 256 ) {
    die "Error, cmd died with ret $ret\n";
  }
  return $ret;
}
