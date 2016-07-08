#!/usr/bin/perl

=pod

=head1 NAME

genome_to_taxonomy.pl

=head1 USAGE

For each sequence, identify genes and what taxonomy they represent

 Mandatory:
  -fasta          :s  FASTA of genome (contigs or scaffolds)
 
 Optional:
  -cpus           :i  Number of CPUs/threads (def. 6).

=head1 AUTHORS

 Alexie Papanicolaou

        Hawkesbury Institute for the Environment
        alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

Copyright 2015 the Western Sydney University
See LICENSE file for license info
It is provided "as is" without warranty of any kind.


=cut

use strict;
use warnings;
use threads;
use List::Util 'shuffle';
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
use Thread_helper;
use Bio::LITE::Taxonomy::NCBI::Gi2taxid;

$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/:$RealBin/../3rd_party/augustus/bin/";

my ($blastp_exec,$blast_taxonomy_summary_exec,$augustus_exec,$ParaFly_exec) = &check_program('blastp','blast_taxonomy_summary.pl','augustus','ParaFly');
my $cpus = 6;
my $fasta;

GetOptions(
	'cpus|threads=i' => \$cpus,
	'fasta|in=s'=>\$fasta

);

$fasta = shift if !$fasta;

pod2usage unless $fasta;

if (!-d "$fasta.checks"){
	mkdir("$fasta.checks");
	my $orig_sep = $/;
	$/ = ">";
	open (FASTA,$fasta);
	while (my $record=<FASTA>){
		chomp($record);
		next unless $record;
		my @lines = split("\n",$record);
		my $id_line = shift(@lines);
		my $id = $id_line;
		$id =~s/\s+.+$//;
		open (SEQ,">$fasta.checks/$id.seq");
		print SEQ ">$id\n".join("\n",@lines);
		close SEQ;
	}
	close (FASTA);
	$/ = $orig_sep;
}
chdir("$fasta.checks") || die $!;

if (!-s "cmds"){
	my @scaff_files = glob("*seq");
	open (CMDS,">cmds");
	foreach my $file (@scaff_files){
		print CMDS "if [ ! -s $file.data ]; then $augustus_exec --genemodel=complete "
                ." --AUGUSTUS_CONFIG_PATH=$RealBin/../3rd_party/augustus/config"
	        ." --alternatives-from-sampling=false --progress=false --gff3=on --UTR=off --noInFrameStop=true --species=fly $file > $file.data ; fi \n";
	}
	close CMDS;
}

system ("$ParaFly_exec -c cmds -CPU $cpus -shuffle -failed_cmds cmds.failed");

my @aug_files = glob("*data");
die unless @aug_files;
my $switch = -s "../augustus_simple_preds.gff" ? 1 : undef;
open (MASTERGFF,">../augustus_simple_preds.gff") unless $switch;
my $i=1;
my $total = scalar(@aug_files);
foreach my $file (@aug_files){
	print "Processing file $i of $total        \r";$i++;
	next unless $file && -s $file || -s "$file.pep";
	my $counter = 1;
	my $id = $file;
	$id=~s/.seq.data//;
	open (IN,$file);
	open (OUT,">$file.pep");
	while (my $ln=<IN>){
		print MASTERGFF $ln unless $switch;
		next unless $ln=~/^\#/;
		if ($ln=~/^\# protein sequence/){
			$ln=~s/\#\s+protein sequence = \[//;
			chomp($ln);
			my $seq = $ln;
			if ($seq=~/\]$/){
				$seq=~s/\]$//;
			}else{
				while (my $ln2=<IN>){
					print MASTERGFF $ln2 unless $switch;
					$ln2=~s/^\#\s+//;
					chomp($ln2);
					$seq .=$ln2;
					last if $seq=~s/\]$//;
				}
			}
			print OUT ">$id.gene.$counter length=".length($seq)."\n$seq\n";
			$counter++;
		}
	}
	close IN;
	close OUT;
	if (!-s "$file.pep"){
		#warn "File $file had no protein data\n";
		unlink("$file.pep");
		next;
	}
}
close MASTERGFF unless $switch;

my @pep_files = glob("*.pep");
my $thread_helper = new Thread_helper(int($cpus/2));
$i=0;
$total = scalar(@pep_files);
foreach my $file (@pep_files){
	next unless $file && -s $file;
	if (!-s "$file.trim"){
		open (IN,$file);
		my $orig_sep = $/;
		$/=">";
		my @records; # shuffle sequences
		while (my $record=<IN>){
			chomp($record);
			next if !$record || length($record) < 135;
			push(@records,'>'.$record);
		}
		@records = shuffle(@records);
		my $counter=int(0);
		open (OUT,">$file.trim");
		foreach my $record (@records){
			print OUT $record;
			$counter++ if length($record) >= 135; # assume up to 35 chars of headers
			last if $counter == 10;
		}
		close IN;
		close OUT;
		$/ = $orig_sep;
	}
	$file .= ".trim";
	next unless -s $file;
        my $thread        = threads->create('run_search',$file);
        $thread_helper->add_thread($thread);
	$i++;
	print "Processed file $i of $total        \r";
}
$thread_helper->wait_for_all_threads_to_complete() if $i && $i>0;
chdir("../");

print "\nDone.\n";

#################################################################
sub run_search(){
	my $file = shift;
	my $outfile = $file."_vs_swissprot";
	my $check = &check_blast_out($outfile);
	system("$blastp_exec -db $RealBin/../databases/blastdb/uniprot_sprot -evalue 1e-30 -query $file -out $outfile -parse_deflines -num_threads 2 -outfmt 5 -max_target_seqs 50") unless $check == 1;
	system("$blast_taxonomy_summary_exec -xml -top 10 -in $outfile >/dev/null 2>/dev/null");
}

sub check_blast_out(){
	my $file = shift;
	my $success=int(0);
	return int(0) unless -s $file;
	open (IN,$file) || die $!;
	while (my $ln=<IN>){
		if ($ln =~/\<\/BlastOutput\>/){
			$success++;
		}
	}
	close IN;
	return $success;
}

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
