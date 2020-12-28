#!/usr/bin/env perl


use strict;
use warnings;
use Data::Dumper;

my $merge_mult_uniq = shift;
my $verbose = 1;
my ( $samtools_exec ) = &check_program( 'samtools' );

my @files = sort glob("*bam");
mkdir ('merged') if !-d 'merged';

my %hash;
foreach my $file (@files){
	my $pattern = $merge_mult_uniq ? '^(\S+)_L\d+(_\S+)\.\S+\.bam' : '^(\S+)_L\d+(_\S+)\.bam$';

	if ($file=~/$pattern/){
		my $new_file = $1.$2.'.bam';
		push(@{$hash{$new_file}},$file);
	}
}



foreach my $outfile (sort keys %hash){
	my $infiles = join(' ',sort @{$hash{$outfile}});
	&process_cmd($samtools_exec." merge --threads 10 -l 9 merged/$outfile ".$infiles);
	&process_cmd($samtools_exec." index merged/$outfile") if -s "merged/$outfile";
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

