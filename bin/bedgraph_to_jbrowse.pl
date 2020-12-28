#!/usr/bin/env perl


use strict;
use warnings;
use Data::Dumper;
use JSON;

# expect: filename\tall other columns

my $xyplot = 1;
my $tsv = shift || die ("Give me a TSV\n");
my $dir = shift;
$dir = './' if !$dir;

my @bw_files = sort glob("gsnap.*.uniq.coverage.bw");
die "No gsnap.*.uniq.coverage.bw files found\n" unless scalar(@bw_files)>0;

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
	die "I don't know which file is this mask (gsnap.$bw_mask) for:\n   $ln\n" if !$bw_filename;
	
	for (my $i=0;$i<scalar(@data);$i++){
		$hash{$bw_filename}{$headers[$i]} = $data[$i] if $data[$i] && $data[$i]!~/^\s*$/;
	}


}
close IN;

#die Dumper \%hash;
my (@json_array);
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
   push(@json_array,\%hash_item);
}

print to_json(\@json_array, {utf8 => 1, pretty => 1});
