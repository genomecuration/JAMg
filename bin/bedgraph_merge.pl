#!/usr/bin/env perl


use strict;
use warnings;
use Data::Dumper;

my $verbose = 1;
my ( $bigwigmerge_exec ) = &check_program( 'bigWigMerge' );

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
	&process_cmd($bigwigmerge_exec." ".$infiles." merged/$outfile");
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

