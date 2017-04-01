#!/usr/bin/perl

=pod

=head1 NAME

 split_gff_to_genbank.pl

=head1 USAGE

 Please provide a GFF (without any ##gff-version metadata etc) and a FASTA file
 Third option could be the GFF record separator (defaults to two newlines)


=cut


use strict;
use Pod::Usage;
use warnings;

my $gff_file = shift;
my $fasta_file = shift ;
my $gff_record = shift;
$gff_record = "\n\n" if !$gff_record ;
die "Please provide a GFF (without any ##gff-version metadata etc) and a FASTA file\n" unless $gff_file && -s $gff_file && $fasta_file && -s $fasta_file;
my ($seqret_exec) = ('seqret');

my $orig_sep = $/;

print "Processing FASTA $fasta_file\n";
my %seq_hash;
$/ = ">";
open (FASTA,$fasta_file)||die;
while (my $record = <FASTA>){
	chomp($record);
	next unless $record;
	my @lines = split("\n",$record);
	my $ref_id = shift(@lines);
	if ($ref_id=~/(\S+)/){
		$ref_id = $1;
	}else{
		die "Invalid format for $fasta_file at $record\n";
	}
	$ref_id=~s/\s+//g;
	my $seq = join("\n",@lines);
	$seq_hash{$ref_id} = $seq;
}
close FASTA;
$/=$orig_sep;

print "Processing GFF $gff_file\n";
my %gff_hash;
open (GFF,$gff_file)||die;
$/ = $gff_record;
while (my $record = <GFF>){
	chomp($record);
	next unless $record;
	my @lines = split("\n",$record);
	my @data = split("\t",$lines[0]);
	my $ref_id = $data[0];
	$ref_id=~s/\s+//g;
	$gff_hash{$ref_id} .= $record.$gff_record;
}
close GFF;
$/=$orig_sep;

print "Producing output file in $gff_file.dir\n";
mkdir("$gff_file.dir") unless -d "$gff_file.dir";
my $counter = int(0);
my $total = scalar(keys %gff_hash);
foreach my $ref_id (keys %gff_hash){
	$gff_hash{$ref_id}=~s/\s+$//;
	die "Cannot find reference sequence $ref_id in FASTA\n" unless $seq_hash{$ref_id};
	open (FSA,">$gff_file.dir/$ref_id.fsa");
	print FSA ">$ref_id\n".$seq_hash{$ref_id};
	close FSA;
	open (OUT,">$gff_file.dir/$ref_id.gff3");
	print OUT "##gff-version 3\n".$gff_hash{$ref_id}."\n##FASTA\n>$ref_id\n".$seq_hash{$ref_id};
	close OUT;
	if (-s "$gff_file.dir/$ref_id.gff3"){
		system("$seqret_exec $gff_file.dir/$ref_id.gff3 $gff_file.dir/$ref_id.gb -feature -osformat2 genbank -sformat1 gff3 -auto");
	}
	$counter++;
	print "COMPLETED $counter / $total    \r" if $counter %10 == 0;
}
system("cat $gff_file.dir/*gb > $gff_file.gb");
print "\nDone!\n";


#############################
sub check_program() {
 my @progs = @_;
 my @paths;
 foreach my $prog (@progs) {
  my $path = `which $prog`;
  pod2usage "Error, path to a required program ($prog) cannot be found\n\n"
    unless $path =~ /^\//;
  chomp($path);
  $path = readlink($path) if -l $path;
  push( @paths, $path );
 }
 return @paths;
}

