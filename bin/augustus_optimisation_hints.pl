#!/usr/bin/env perl

=pod

=head1 NAME

 augustus_optimisation_hints.pl

=head1 USAGE

Create hint files for Augustus optimisation files using those created for the genome sequence. 

Mandatory options:

 -hints|in         s{1,}  	The hints file you want to process
 -fasta		   s  		The FASTA file as created by prepare_golden (e.g. final_golden_genes.gff3.nr.golden.optimization.good.gb.fasta)

=head1 AUTHORS

 Alexie Papanicolaou

        Hawkesbury Institute for the Environment
        alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

Copyright 2016 the Western Sydney University
See LICENSE file for license info
It is provided "as is" without warranty of any kind.

=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util qw(sum);
use Pod::Usage;
use File::Basename;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";
use Fasta_reader;

#Options
my ( @hintfiles, $fasta_file, $help );

my $cpus = 4;
my $sort_buffer = '5G';
my $tmpdir = $ENV{'TMP'};
$tmpdir = $ENV{'TMPDIR'} if !$tmpdir;
$tmpdir = '/tmp' if !$tmpdir;
pod2usage $! unless &GetOptions(
            'help'              => \$help,
            'hints|in:s{,}'     => \@hintfiles,
            'fasta:s'           => \$fasta_file,
	    'cpus:i'            => \$cpus,
	    'memory:s'          => \$sort_buffer,
	    'tmp:s'             => \$tmpdir
);

my $sort_exec = &check_sort_version;

pod2usage if $help;

pod2usage "Cannot find the hints or genome FASTA file\n"
  unless $hintfiles[0]
   && -s $hintfiles[0]
   && $fasta_file
   && ( -s $fasta_file);

print STDERR "Acquiring reference data\n";
my (%ref_data);
my $orig_sep = $/;
$/ =">"; 
#==> complete_reference_training_set.fasta <==
#>Herato0101_11477519-11486756
my $total_count=int(0);
open (FASTA,$fasta_file)||die;
while (my $record=<FASTA>){
	chomp($record);next unless $record;
	my @lines = split("\n",$record);
	my $id = shift (@lines);
	if ($id=~/^(\S+)_(\d+)-(\d+)$/){
		my ($ref,$start,$stop) = ($1,$2,$3);# always + strand
		$ref_data{$ref}{$start}{$stop}++;
		$total_count++;
	}

}
close FASTA;
$/ = $orig_sep;
print STDERR "Found $total_count values\n\n";
# ==> master_bamfile2.bam.rnaseq.hints <==
#Herato0101      RNASeq  exonpart        436     469     21      .       .       src=RCOV;pri=4
$|=1;


foreach my $hint_file (@hintfiles){
	print STDERR "Parsing hints file $hint_file\n";
	my $line_max = `wc -l < $hint_file`;chomp($line_max);
	my $line_count = int(0);
	open (IN,$hint_file);
	while (my $ln =<IN>){
		$line_count++;
		my @data = split("\t",$ln);
		next unless $data[8];
		if ($data[3] > $data[4]){
			my $t = $data[4];
			$data[4] = $data[3];
			$data[3] = $t;
			$data[6] = '-';
		}
		die "Unexpected start > stop in hints file:\n$ln\n" if $data[3] > $data[4] ;
		my ($ref,$hint_start,$hint_stop) = ($data[0],$data[3],$data[4]);
		# if we have a reference sequence match
		if ($ref_data{$ref}){
			# and the hint is part of our dataset
			foreach my $gene_start (keys %{$ref_data{$ref}}){
				next if $hint_start < $gene_start;
				foreach my $gene_stop (keys %{$ref_data{$ref}{$gene_start}}){
					last if $hint_stop > $gene_stop;
					if ($hint_start >= $gene_start && $hint_stop <= $gene_stop){
						my $opt_hint_start = $hint_start - $gene_start +1;
						my $opt_hint_stop  = $hint_stop - $gene_start +1;
						my @opt_data = @data;
						$opt_data[0] = $ref."_$gene_start"."-$gene_stop";
						$opt_data[3] = $opt_hint_start;
						$opt_data[4] = $opt_hint_stop;
						print join("\t",@opt_data);
					}
				}
			}
		}
		my $prog = int(($line_count / $line_max) * 100);
		print STDERR "$prog %  \r" if $line_count % 10000 == 0;
	}
	close IN;
	print STDERR "100 % \n";
}
$|=0;


######################################
sub check_sort_version(){
	my ($sort_exec) = &check_program('sort');
	my @v=`$sort_exec --version`;

	if ($v[0] && $v[0]=~/(\d+)\.(\d+)\s*$/){
		my $major = $1;
		my $minor = $2;
		if ($major >= 8 && $minor >= 6){
			return "$sort_exec -T $tmpdir --parallel $cpus -S $sort_buffer";
		}else{
			return "$sort_exec -S $sort_buffer -T $tmpdir";
		}
	}else{
		die "Sort of coreutils not found!";
	}
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

###
sub process_cmd {
 my ($cmd) = @_;
 print STDERR "CMD: $cmd\n";
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  die "Error, cmd died with ret $ret\n";
 }
 return $ret;
}

sub merge_hints(){
 my $file = shift;
 open (IN,$file);
 open (OUT,">$file.merged");
 my (@current_line,@previous_line);
 while (my $ln=<IN>) {
    @current_line = split ("\t",$ln);
    if (!@previous_line){
        @previous_line = @current_line;
    }elsif(($current_line[0] eq $previous_line[0]) && ($current_line[2] eq $previous_line[2]) && 
    (($current_line[3] >= $previous_line[3]) && ($current_line[4] <= $previous_line[4]))
      && ($current_line[6] eq $previous_line[6])){
     # update previous_line by adding current to it
        chomp($previous_line[8]);
        $previous_line[8] =~ s/(grp=[^;]*);*//;
        my $grp = $1;
        $grp .= ';' if $grp;
        $grp = '' if !$grp;
        my ($lm,$m)=(1,1);
        if ($previous_line[8] =~ /mult=(\d+);/){
            $lm = $1;
            $previous_line[8] =~ s/mult=\d+;//;
        }
        if ($current_line[8] =~ /mult=(\d+);/){
            $m = $1;
        }
        $previous_line[8] = "mult=" . ($lm+$m) . ";$grp" . $previous_line[8]."\n";
     
    }elsif (
	    !(($current_line[0] eq $previous_line[0]) && ($current_line[2] eq $previous_line[2]) 
	&& ($current_line[3] == $previous_line[3]) && ($current_line[4] == $previous_line[4])  
	&& ($current_line[6] eq $previous_line[6]))){
        print OUT join("\t",@previous_line);
        @previous_line = @current_line;
    }
    
     else {
        # update previous_line by adding current to it
        chomp($previous_line[8]);
        $previous_line[8] =~ s/(grp=[^;]*);*//;
        my $grp = $1;
        $grp .= ';' if $grp;
        $grp = '' if !$grp;
        my ($lm,$m)=(1,1);
        if ($previous_line[8] =~ /mult=(\d+);/){
            $lm = $1;
            $previous_line[8] =~ s/mult=\d+;//;
        }
        if ($current_line[8] =~ /mult=(\d+);/){
            $m = $1;
        }
        $previous_line[8] = "mult=" . ($lm+$m) . ";$grp" . $previous_line[8]."\n";
    }
 }
  print OUT join("\t",@previous_line) if (@previous_line);
  close IN;
  close OUT;
  unlink($file);
  rename($file.'.merged',$file);
}

sub touch() {
 my $file = shift;
 system("touch $file");
}
