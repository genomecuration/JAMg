#!/usr/bin/env perl

=pod

=head1 NAME hmmscandom2hints

=head1 USAGE

	-exonfile 	:s	FASTA with input to hmmscan
	-domains	:s{1,}	One or more DOMAIN output files from hmmscan
	-output		:s	Basename for output files
	-score_cut	:f	Cutoff BIT score for processing (defaults 10.0)

=cut

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;

my $exonfile;
my $base_output;
my @domainfiles;
my $score_cutoff = 10.0;


&GetOptions(
	'exonfile:s' => \$exonfile,
	'output:s' => \$base_output,
	'domains:s{1,}' => \@domainfiles,
	'score_cut:f' => \$score_cutoff,
);



my %retrotransposon_IDS = (
'PF14529'=>1,
'PF11474'=>1,
'PF00078'=>1,
'PF07727'=>1,
'PF07727'=>1,
'PF13456'=>1,
'PF06815'=>1,
'PF13655'=>1,
'PF06817'=>1,
'PF13966'=>1,
'PF02495'=>1,
'PF08475'=>1
);

if (@domainfiles && scalar(@domainfiles) == 1 && -d $domainfiles[0]){
	@domainfiles = glob($domainfiles[0]."/*dom");
}

pod2usage "Please provide the exon FASTA input to HHMSEARCH\n" unless $exonfile && -s $exonfile;
pod2usage "Please provide the files that have domain output of HHMSEARCH or a directory with .dom files\n" unless $domainfiles[0] && -s $domainfiles[0];

$base_output = "$exonfile.domains" if !$base_output;
pod2usage "Base output files already exist $base_output*\n" if -s "$base_output.gff3" || -s "$base_output.hints";


my %hash;

open (EXON,$exonfile);
print "Building lookup from $exonfile\n";

while (my $ln = <EXON>){
	next unless $ln=~/^>/;
	my $strand = $ln=~/REVERSE SENSE/ ? '-' : '+';
	if ($ln=~/^>(\S+)\s\[(\d+)\s\-\s(\d+)\]/){
		my $id = $1;
		my $start = ($strand eq '+') ? $2 : $3;
		my $end = ($strand eq '+') ? $3 : $2;
		$hash{$id} = "$strand ".$start."\t".$end;
	}
}
close EXON;


open (GFF3,">$base_output.gff3") ||die;
open (HINTS,">$base_output.hints") ||die;

foreach my $domainfile (@domainfiles){
	next unless -s $domainfile;
	print "Processing domain HMMSEARCH file $domainfile\n";
	open (IN, $domainfile)||die $!;

  while (my $ln=<IN>){
	next if $ln=~/^#/;
	die "Wrong file format for hhmsearch output ($domainfile)\n" if $ln=~/^>/;
	chomp($ln);
	my @data = split(/\s+/,$ln);
	next unless $data[13];
	my $score = $data[7];
	next unless $score >= $score_cutoff;
	my $name = $data[0];
	my $desc = $data[-1];
	next if $name=~/transcriptase/i;
	next if $desc=~/transcriptase/i;
	my $hit = $data[1];
	$hit = $name if $hit eq '-';
	my $check = $hit;
	if ($check=~/^PF/){
		$check=~s/\.\d+$//;
		next if $retrotransposon_IDS{$check};
	}

	my $id = $data[3] || next;
	my $start_end = $hash{$id} || next;
	my $strand;
	if ($start_end=~s/^(\S)\s//){
		$strand = $1;
	}else{next;}
	$id =~s/_\d+$//;

	my ($start,$seqend) = split("\t",$start_end);
	$start+=$data[19];
	my $end = $start + ($data[20] - $data[19]);
	next if $end > $seqend;

	print GFF3 "$id\thmmsearch\tprotein_match\t".$start."\t$end"
	."\t$score\t$strand\t.\tID=$hit;Name=$name\n";

	print HINTS "$id\thmmsearch\tCDSpart\t"     .$start."\t$end"
	."\t$score\t$strand\t.\tsrc=HU;grp=$hit;prio=2\n";
  }
  close IN;
}

close GFF3;
close HINTS;
print "Done, see $base_output.* \n";
