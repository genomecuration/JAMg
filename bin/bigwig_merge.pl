#!/usr/bin/env perl


use strict;
use warnings;
use Data::Dumper;


my $genome_path = shift || die("Please provide the path to the genome FASTA\n");
my $verbose = 1;

my ( $bigwigmerge_exec,$bedGraphToBigWig_exec, $samtools_exec ) = &check_program( 'bigWigMerge','bedGraphToBigWig', 'samtools' );

&process_cmd($samtools_exec . " faidx $genome_path") unless -s "$genome_path.fai";

my @files = sort glob("*coverage.bw");
mkdir ('merged') if !-d 'merged';

my %hash;
foreach my $file (@files){
	if ($file=~/^(\S+)_L\d+(_\S+)$/){
		my $new_file = $1.$2;
#		print "Merging $file into merged/$new_file\n";
		push(@{$hash{$new_file}},$file);
	}
}



foreach my $outfile (sort keys %hash){
	my $infiles = join(' ',sort @{$hash{$outfile}});
	&process_cmd($bigwigmerge_exec." $infiles merged/$outfile.bg");
	&process_cmd($bedGraphToBigWig_exec. " merged/$outfile.bg merged/$outfile.bw $genome_path.fai merged/$outfile.bw");
	unlink("merged/$outfile.bg");
} 





####################################################
sub check_program() {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  die "Error, path to required $prog cannot be found\n"
    unless $path =~ /^\//;
  chomp($path);
  #$path = readlink($path) if -l $path;
  push( @paths, $path );
 }
 return @paths;
}

sub process_cmd {
 my ($cmd) = @_;
 print "CMD: $cmd\n" if $verbose;
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  die "Error, cmd died with ret $ret\n";
 }
 return $ret;
}

