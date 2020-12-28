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
mkdir("$jsondir/bw/") if !-d "$jsondir/bw/";

my $bw_dir = shift;
$bw_dir = './' if !$bw_dir;

my @bw_files = sort glob("gsnap.*.uniq.coverage.bw");
die "No gsnap.*.uniq.coverage.bw files found\n" unless scalar(@bw_files)>0;


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
	my ($bw_mask,$bw_filename);
	if ($filename =~/^(\S+)_L\d+_R\d+/){
		$bw_mask = $1.'_vs_';
	}
	foreach my $bw (@bw_files){
		$bw_filename = $bw if (!$bw_filename && $bw=~/gsnap\.$bw_mask/);
	}

	if ($bw_filename && -s $bw_filename){
		for (my $i=0;$i<scalar(@data);$i++){
			$hash{$bw_filename}{$headers[$i]} = $data[$i] if $data[$i] && $data[$i]!~/^\s*$/;
		}
	}

}
close IN;

print "Copying bigwig and backing up in $jsondir\n";
rename($injson_file,$injson_file.".bak");

#die Dumper \%hash;
foreach my $bw_filename (sort keys %hash){

#         "urlTemplate" => $url_prefix."/$organism/data/bigwig/".$align_file
   my %hash_item = (
         "bicolor_pivot" => "zero",
         "storeClass" => "JBrowse/Store/SeqFeature/BigWig",
         "label" => "RNA-Seq of " . $hash{$bw_filename}{'friendlyname'}  . " coverage",
         "key" => "RNA-Seq_gsnap-coverage_" . $hash{$bw_filename}{'friendlyname'},
         "autoscale" => "local",
         "category" => "RNA-Seq coverage graph",
 	 "urlTemplate" => "bw/$bw_filename",
	);

	if ($xyplot){
	        $hash_item{"type"} = "JBrowse/View/Track/Wiggle/XYPlot";
	}else{
	        $hash_item{"type"} = "JBrowse/View/Track/Wiggle/Density";
	}

	for (my $i=2;$i<scalar(@headers);$i++){
		$hash_item{"metadata"}{$headers[$i]} = $hash{$bw_filename}{$headers[$i]} if $hash{$bw_filename}{$headers[$i]};
	}
   push(@json_lines,\%hash_item);
   link($bw_filename,"$jsondir/bw/$bw_filename");
   die $! unless -s "$jsondir/bw/$bw_filename";
}


push(@{$json_track_hashref->{'tracks'}},@json_lines);
open (OUT,">$injson_file") || die $!;
print OUT to_json($json_track_hashref, {utf8 => 1, pretty => 1});
close (OUT);

print "Done see $injson_file\n";
