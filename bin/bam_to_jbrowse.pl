#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use JSON;

# expect: filename\tall other columns

my $xyplot = 1;
my $tsv = shift || die ("Give me a TSV and an input JSON\n");
my $injson_file = shift || die ("Give me a TSV and an input JSON\n");
my $jsondir = dirname($injson_file) || die $!;
die unless -d $jsondir;
mkdir("$jsondir/bam/") if !-d "$jsondir/bam/";

my $bam_dir = shift;
$bam_dir = './' if !$bam_dir;

my @bam_files = sort glob("gsnap.*.merged.bam");
die "No gsnap.*.merged.bam files found\n" unless scalar(@bam_files)>0;


my ($json_track_hashref,@json_lines);
if (-s $injson_file){
	open (IN, $injson_file) || die $!;
	my @a = <IN>;
	my $str = join('',@a);
	$json_track_hashref = decode_json $str;
	close IN;
}else{
	die "The $injson_file is empty\n";
}




my %hash;
open (IN, $tsv) || die $!;
my @headers = split("\t",<IN>);
die  "TSV file should have a header line. First column should be filename or file"
     ." and second column should be name or friendly name (case insensitive)\n" unless $headers[0]=~/^file\s?n?a?m?e?$/i && $headers[1]=~/^f?r?i?e?n?d?l?y?\s*name$/i;
chomp($headers[-1]);
$headers[0] = 'filename';
$headers[1] = 'friendlyname';


while (my $ln=<IN>){
	chomp($ln);next unless $ln;
	next if $ln=~/^#/ || $ln=~/^\s*$/;
	my @data = split("\t",$ln);
	die "Weird format for $tsv\n" if !$data[1];

	my $filename = $data[0];
	my ($bam_mask,$bam_filename);
	if ($filename =~/^(\S+)_L\d+_R\d+/){
		$bam_mask = $1.'_vs_';
	}
	foreach my $bam (@bam_files){
		$bam_filename = $bam if (!$bam_filename && $bam=~/gsnap\.$bam_mask/);
	}

	if ($bam_filename && -s $bam_filename){
		for (my $i=0;$i<scalar(@data);$i++){
			$hash{$bam_filename}{$headers[$i]} = $data[$i] if $data[$i] && $data[$i]!~/^\s*$/;
		}
	}

}
close IN;

print "Copying bigwig and backing up in $jsondir\n";
rename($injson_file,$injson_file.".bak");

#die Dumper \%hash;
foreach my $bam_filename (sort keys %hash){

#         "urlTemplate" => $url_prefix."/$organism/data/bam/".$file,

  my %hash_item = (
         "type" => "JBrowse/View/Track/Alignments2",
         "key" => "RNA-Seq_gsnap-align_" . $hash{$bam_filename}{'friendlyname'},
         "label" => "RNA-Seq of " . $hash{$bam_filename}{'friendlyname'},
         "storeClass" => "JBrowse/Store/SeqFeature/BAM",
         "urlTemplate" => "bam/$bam_filename",
         "category" => "Read alignments"
  );

   for (my $i=2;$i<scalar(@headers);$i++){
	$hash_item{"metadata"}{$headers[$i]} = $hash{$bam_filename}{$headers[$i]} if $hash{$bam_filename}{$headers[$i]};
   }
   push(@json_lines,\%hash_item);
   link($bam_filename,"$jsondir/bam/$bam_filename");
   die $! unless -s "$jsondir/bam/$bam_filename";
}


push(@{$json_track_hashref->{'tracks'}},@json_lines);
open (OUT,">$injson_file") || die $!;
print OUT to_json($json_track_hashref, {utf8 => 1, pretty => 1});
close (OUT);

print "Done see $injson_file\n";




