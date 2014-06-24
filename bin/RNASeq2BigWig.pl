#!/usr/bin/env perl

=pod

=head1 USAGE

Create BigWig files for JBrowse

Mandatory options:

 -bam|in           s  The input BAM file (co-ordinate sorted).
 -genome|fasta     s  The genome assembly FASTA file.

NB: Needs bedGraphToBigWig from Kent's src code

 git clone git://genome-source.cse.ucsc.edu/kent.git


=head1 AUTHORS

 Alexie Papanicolaou

        CSIRO Ecosystem Sciences
        alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

Copyright 2012-2014 the Commonwealth Scientific and Industrial Research Organization.
See LICENSE file for license info
It is provided "as is" without warranty of any kind.

=cut


use strict;
use warnings;
use FindBin qw($RealBin);
use Pod::Usage;
use Getopt::Long;
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";

my ( @bamfiles, $genome, $help );

pod2usage $! unless &GetOptions(
            'help'              => \$help,
            'bam|in:s{,}'          => \@bamfiles,
            'genome|fasta:s'    => \$genome,
);

pod2usage if $help;
pod2usage "Cannot find the BAM or genome FASTA file\n"
  unless $bamfiles[0]
   && -s $bamfiles[0]
   && $genome
   && ( -s $genome || -s $genome . '.fai' );

my ( $samtools_exec, $bedtools_exec, $bigwig_exec ) = &check_program( 'samtools', 'bedtools','bedGraphToBigWig' );

&process_cmd("$samtools_exec faidx $genome") unless -s $genome . '.fai';
die "Cannot index genome $genome\n" unless -s $genome . '.fai';


my $master_bamfile;
if (scalar(@bamfiles == 1)){
  $master_bamfile = $bamfiles[0];
}else{
        foreach my $bamfile (@bamfiles){
                die "Cannot find $bamfile\n" unless -s $bamfile;
        }
        $master_bamfile = 'master_bamfile.bam';
        &process_cmd("$samtools_exec merge -\@2 -r $master_bamfile ".join(" ",@bamfiles)) unless -s $master_bamfile;
}


&process_cmd("$bedtools_exec genomecov -split -bg -g $genome.fai -ibam $master_bamfile| sort -S 1G -k1,1 -k2,2n > $master_bamfile.coverage.bg");
&process_cmd("$bigwig_exec $master_bamfile.coverage.bg $genome.fai $master_bamfile.coverage.bw");

print "Done. See $master_bamfile.coverage.bw\n";



###
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
###
sub process_cmd {
 my ($cmd) = @_;
 print "CMD: $cmd\n";
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  die "Error, cmd died with ret $ret\n";
 }
 return $ret;
}

