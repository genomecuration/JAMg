#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;

=pod

USAGE

 Multiple files:

 samtools merge -u -@ 4 - ../*bam | samtools view - | rename_bam.pl rename.tsv | samtools view -b -T $GENOME_PATH -o bamout.bam -

 Single file:

 samtools view bamfile | rename_bam.pl rename.tsv | samtools view -b -T $GENOME_PATH -o bamout.bam -

=cut

my $list = shift || pod2usage;

my @files = @ARGV;
if (!$files[0]){
	push(@files,"/dev/stdin");
}

my %hash;
open (IN,$list)||die;
while (my $ln=<IN>){
	chomp($ln);
	my @data = split("\t",$ln);
	$hash{$data[0]} = $data[1] if $data[1];
}
close IN;

foreach my $bam (@files){

	open (BAM,$bam)||die;
	my $output = $bam ne '/dev/stdin' ? $bam.'.renamed' : "/dev/stdout";
	open (OUT,">$output")||die "$bam to $output: ".$!;

	while (my $ln=<BAM>){
		chomp($ln);
		my @data = split("\t",$ln);
		my $old_id = $data[2];
		my $new_id = $hash{$data[2]};
		#die "The new id for $old_id not found\n" unless $new_id || $old_id eq '*';
		next unless $new_id || $old_id eq '*';
		$data[2] = $new_id if $new_id;
		my $pair_old = $data[6];
		my $pair_new = $hash{$data[6]};
		$data[6] = $hash{$data[6]} if $hash{$data[6]}; 
		while ($data[-1]=~/(\w+):\d+/g){
			my $o = $1 || next;
			my $n = $hash{$o} || next;
			$data[-1]=~s/$o/$n/;
		}

		print OUT join("\t",@data)."\n";
	}
	close BAM;
	close OUT;
}

