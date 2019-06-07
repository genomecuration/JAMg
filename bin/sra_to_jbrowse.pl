#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use JSON;

my $aligns_dir = 'aligns';
my $url_prefix = 'https://s3b.hortgenomics.science/ntg/jbrowse';
my $injson_file = 'data/trackList.json';
my $xyplot = 1;
my $organism;

my $sra_description_file = shift || die("Provide SRA description table file\n");
die "No INPUT JSON FILE found ($injson_file)" unless $injson_file && -s $injson_file;

print STDERR "Using $aligns_dir for finding alignments\nURL Prefix is $url_prefix\nINPUT JSON is $injson_file\n";
print STDERR "Coverage graphs will be XY plots\n" if $xyplot;
print STDERR "Coverage graphs will be Density plots\n" if !$xyplot;


open (IN, $injson_file);
my @a = <IN>;
my $str = join('',@a);
my $json_track_hashref = decode_json $str;
close IN;
my @json_lines;

open (SRA,$sra_description_file) || die $!;
open (OUT,">$sra_description_file.json");

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
	if (!$organism || $organism=~/from\s*data/i){
		$organism = $data[$header_lookup{'Organism'}];
		$organism =~s/\s+/_/g;
		print STDERR "Organism set as $organism\n";
	}
	
	my $align_file = "gsnap.".$data[$header_lookup{'Run'}]."_vs_".$organism.".concordant_uniq.coverage.bw";
	# single-end
	if (!-s $aligns_dir.'/'.$align_file){
		$align_file = "gsnap.".$data[$header_lookup{'Run'}]."_vs_".$organism.".unpaired_uniq.coverage.bw";
	}

	if (!-s $aligns_dir.'/'.$align_file){
		warn "File not found for: ".$data[$header_lookup{'Run'}]."\n";
		next;
	}

	my %hash_item = (
         "bicolor_pivot" => "zero",
         "storeClass" => "JBrowse/Store/SeqFeature/BigWig",
         "label" => $data[$header_lookup{'Assay_Type'}]."_gsnap-coverage_".$data[$header_lookup{'Run'}],
         "autoscale" => "local",
         "key" => $data[$header_lookup{'Assay_Type'}]." of ".$data[$header_lookup{'Run'}]." coverage",
         "category" => $data[$header_lookup{'Assay_Type'}]." coverage graph",
         "urlTemplate" => $url_prefix."/$organism/data/bigwig/".$align_file
	);

	if ($xyplot){
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
	}

	push(@json_lines,\%hash_item);
}
close SRA;

my $junction_file_base = 'master_bamfile.bam';
&process_bw($junction_file_base.".coverage.bw",\@json_lines,"RNA-Seq Global coverage","rnaseq_global_coverage");
&process_bam($junction_file_base.".junctions.bam.sorted",\@json_lines,"RNA-Seq junction reads","rnaseq_junctions");

if ($ENV{'GENOME_PATH'} && -s $ENV{'GENOME_PATH'}){
	my $gaps_file = $ENV{'GENOME_PATH'} . '.gap.wig.bw';
	system('genome_gaps_to_bed.pl '.$ENV{'GENOME_PATH'} ) if !-s $gaps_file;
	&process_bw($gaps_file,\@json_lines,"Genome gaps","genome_caps") if -s $gaps_file;
}


push(@{$json_track_hashref->{'tracks'}},@json_lines);

# print to_json($json_track_hashref, {utf8 => 1, pretty => 1});
print OUT to_json($json_track_hashref, {utf8 => 1, pretty => 1});
close OUT;
print "Done, see $sra_description_file.json\n";

####

sub process_bw(){
  my $file = shift;
  my $out_arrayref = shift;
  my $key = shift;
  my $label = shift;

  if (!-s $aligns_dir.'/'.$file){
    warn "BigWig not found ($aligns_dir/$file)\n";
    next;
  }
	my %hash_item = (
         "bicolor_pivot" => "zero",
         "storeClass" => "JBrowse/Store/SeqFeature/BigWig",
         "label" => $label,
         "autoscale" => "local",
         "key" => $key,
         "category" => "Global coverage graphs",
	 "type" => "JBrowse/View/Track/Wiggle/Density",
         "urlTemplate" => $url_prefix."/$organism/data/bigwig/".$file
	);
   push(@$out_arrayref,\%hash_item);

}


sub process_bam(){
  my $file = shift;
  my $out_arrayref = shift;
  my $key = shift;
  my $label = shift;

  if (!-s $aligns_dir.'/'.$file){
    warn "BAM not found ($aligns_dir/$file)\n";
    next;
  }
  my %hash_item = (
         "type" => "JBrowse/View/Track/Alignments2",
         "key" => $key,
         "label" => $label,
         "storeClass" => "JBrowse/Store/SeqFeature/BAM",
         "urlTemplate" => $url_prefix."/$organism/data/bam/".$file,
         "category" => "Read alignments"
  );

   push(@$out_arrayref,\%hash_item);
}

