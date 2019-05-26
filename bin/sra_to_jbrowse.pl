#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use JSON;

my $aligns_dir = 'aligns';
my $url_prefix = 'https://s3b.hortgenomics.science/ntg/jbrowse';
my $injson_file = 'data/trackList.json';
my $plot = 1;

my $sra_description_file = shift || die("Provide SRA description table file\n");
die "No INPUT JSON FILE found ($injson_file)" unless $injson_file && -s $injson_file;

open (IN, $injson_file);
my @a = <IN>;
my $str = join('',@a);
my $json_track_hashref = decode_json $str;
close IN;
my @json_lines;

open (SRA,$sra_description_file) || die $!;
#open (OUT,"$sra_description_file.json");

my $header_str = <SRA>;
chomp($header_str);
my @headers = split("\t",$header_str);
my %header_lookup;
for (my $i=0;$i<@headers;$i++){
  $header_lookup{$headers[$i]} = $i;
}
# use $data[$header_lookup{''}]

while (my $ln=<SRA>){
	chomp($ln);
	next unless $ln;
	my @data = split("\t",$ln);
	my $organism = $data[$header_lookup{'Organism'}];
	$organism =~s/\s+/_/g;

	next if !-s $aligns_dir."/gsnap.".$data[$header_lookup{'Run'}]."_vs_".$organism.".concordant_uniq.coverage.bw";

	my %hash_item = (
         "bicolor_pivot" => "zero",
         "storeClass" => "JBrowse/Store/SeqFeature/BigWig",
         "label" => $data[$header_lookup{'Assay_Type'}]."_gsnap-coverage_".$data[$header_lookup{'Run'}],
         "autoscale" => "local",
         "key" => $data[$header_lookup{'Assay_Type'}]." of ".$data[$header_lookup{'Run'}]." coverage",
         "category" => $data[$header_lookup{'Assay_Type'}]." coverage graph",
         "urlTemplate" => $url_prefix."/$organism/data/bigwig/gsnap.ERR760722_vs_Citrus_maxima.concordant_uniq.coverage.bw",
	);

	if ($plot){
	        $hash_item{"type"} = "JBrowse/View/Track/Wiggle/XYPlot";
	}else{
	        $hash_item{"type"} = "JBrowse/View/Track/Wiggle/Density";
	}
	$hash_item{"metadata"}{"Description"} = $data[$header_lookup{'Assay_Type'}]." public data";
	$hash_item{"metadata"}{"Organism"} = $data[$header_lookup{'Organism'}] if $data[$header_lookup{'Organism'}];
	for (my $i=0;$i<@data;$i++){
		next if ($headers[$i] eq "DATASTORE_region" || $headers[$i] eq "Consent" || $headers[$i] eq "DATASTORE_filetype" || $headers[$i] eq "DATASTORE_provider");
		if ($data[$i] && $data[$i]=~/^[A-Za-z]{2,}/){
			$hash_item{"metadata"}{$headers[$i]} = $data[$i];
		}
		#not really used:
		#push(@{$hash{$headers[$i]}},$data[$i]);
	}

	push(@json_lines,\%hash_item);
}
close SRA;

push(@{$json_track_hashref->{'tracks'}},@json_lines);

print to_json($json_track_hashref, {utf8 => 1, pretty => 1});


####

