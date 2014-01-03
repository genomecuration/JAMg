#!/usr/bin/env perl

# dss gccgagaactccgctcgttctgtgcgttctcctgtcccaggtagggaagaggggctgccgggcgcgctctgcgccccgtttc
# ass gttagtatgcttctttaattttttttctccctgaaattataggaaccagatgttaaaaaattagaagaccaacttcaaggcg
# ass --------------------------ggctttgtctttgcagaatttatagagcggcagcacgcaaagaacaggtattacta

=pod

=head1 NAME

 get_splice_sites.pl

=head1 USAGE

        'hint:s' =>\$hint_file,
        'fasta:'s =>\$fasta_file

=head1 AUTHORS

 Alexie Papanicolaou

        CSIRO Ecosystem Sciences
        alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

Copyright 2012-2014 the Commonwealth Scientific and Industrial Research Organization. 
This software is released under the Mozilla Public License v.2.

It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.mozilla.org/MPL/2.0.

=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";
use Fasta_reader;

my ($hint_file,$fasta_file);

&GetOptions(
	'hint:s' =>\$hint_file,
	'fasta:s' =>\$fasta_file
);

pod2usage unless $hint_file && -s $hint_file && $fasta_file && -s $fasta_file;

my (%site_dss_count,%site_ass_count);
print "Parsing genome $fasta_file\n";
my $fasta_obj = new Fasta_reader($fasta_file);
my %fasta_data = $fasta_obj->retrieve_all_seqs_hash();
print "Parsing data $hint_file\n";
open (IN,$hint_file) || die;
open (OUT,">$hint_file.splice");

while (my $ln=<IN>){
	my @data = split("\t",$ln);
	next unless $data[6] && ($data[2] eq 'ass' || $data[2] eq 'dss');
	my ($ref,$type,$start,$stop,$strand) = ($data[0],$data[2],$data[3],$data[4],$data[6]);
	$start--; # convert to 0-base
	$stop--;
	die "Cannot find sequence for $ref in FASTA file $fasta_file\n" unless $fasta_data{$ref};
	my $site = uc(substr($fasta_data{$ref},$start,2));
	$site = &revcomp($site) if $strand eq '-';

	if ($type eq 'dss'){
		unless ($site eq 'GT' || $site eq 'GC'){
			warn "Unexpected Donor site $site found. Skipping\n";
			next;
		}
		$site_dss_count{$site} ++;
	}elsif($type eq 'ass'){
		if ($site ne 'AG'){
			warn "Unexpected Acceptor site $site found. Skipping\n";
			next;
		}
		$site_ass_count{$site} ++;
	}
	my $site_seq = lc(substr($fasta_data{$ref},($start-40),82)) ||next;
	$site_seq = &revcomp($site_seq) if $strand eq '-';
#debug
#	my $site_seq1 = lc(substr($fasta_data{$ref},($start-40),40));
#	my $site_seq2 = uc(substr($fasta_data{$ref},$start,2));
#	my $site_seq3 = lc(substr($fasta_data{$ref},($start+2),40));
#	my $site_seq = $site_seq1.$site_seq2.$site_seq3;
	print OUT "$type $site_seq\n";;
}
close IN;
close OUT;

print "Completed: Donors:\n";
print Dumper \%site_dss_count;
print "\tAcceptors:\n";
print Dumper \%site_ass_count;


sub revcomp {
  my $dna     = shift;
  my $revcomp = reverse($dna);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}


