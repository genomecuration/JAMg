#!/usr/bin/env perl

=pod

=head1 NAME

process_ncbi_assemblies_html.pl

=head1 USAGE

 Use one of:
    -file|in	:s	=> A assembly table file acquired from NCBI (see DESCRIPTION)
 OR (easiest):
    -taxid	:i	=> NCBI TaxID to download (could be a whole taxon, e.g. Panarthropoda)

=head1 DESCRIPTION

 Give a Taxon or file:

 wget 'https://www.ncbi.nlm.nih.gov/assembly/organism/browse-organism/4751/latest/?p$l=Ajax&page=1&pageSize=2500&sortCol=1&sortDir=1' -O 4751.html

input is expected to be like
<tr class="odd">
  <td>[Candida] aaseri (ascomycetes)</td>
  <td>
    <a href="https://www.ncbi.nlm.nih.gov/assembly/GCA_002068075.1/">cb02</a>
  </td>
  <td>KRIBB</td>
  <td>03/22/2017</td>
  <td>full</td>
  <td>Scaffold</td>
  <td>latest</td>
  <td>representative genome</td>
</tr>

=cut

use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;
use Pod::Usage;
use File::Basename;
use HTML::Strip;
$|=1;

our $VERSION = '0.1';
my ($file,$taxid);

GetOptions(
	'file|in:s' => \$file,
	'taxid:i' =>\$taxid
);

$file = shift unless $file;
if (!$file && $taxid){
	system("wget -q 'https://www.ncbi.nlm.nih.gov/assembly/organism/browse-organism/$taxid/latest/?p\$l=Ajax&page=1&pageSize=99999&sortCol=1&sortDir=1' -O $taxid.html >/dev/null");
	$file = "$taxid.html";
}elsif ($file && $file=~/^(\d+/){
	$taxid = $1;
	print "From file $file : TaxID assumed to be $taxid\n";
}

die "Please give a file to process\n" unless $file && -s $file;

my $data_hash_ref = &process_ncbi_genome_assembly($file);


############################################################################

sub process_ncbi_genome_assembly(){
	my $file = shift;
	my %out_hash;
	open (OUT,">$file.out.tsv");

	my $orig_sep = $/;
	$/ = '<tr';
	open (IN,$file);
	my $total_data_line = <IN>;
	my ($total,$counter)=(int(0),int(0));
	if ($total_data_line=~/\s(\d+)\s-->/){
		$total = $1;
	}
	my $hs = HTML::Strip->new();
	while (my $record=<IN>){
		chomp($record);
		next unless $record;
		my $url;
		if ($record=~/(www.ncbi.nlm.nih.gov[^"]+)/){
			$url = 'https://'.$1;
		}
		next unless $url;
		$record = $hs->parse( $record );
		my @lines = split("\n",$record);
		my $skip = shift(@lines);
		my @data = ($taxid);
		foreach my $ln (@lines){
			next if !$ln || $ln=~/^\s+$/;
			$ln=~s/^\s+//;
			$ln=~s/\s+$//;
			push(@data,$ln);
		}
		push(@data,$url);
		push(@data,&parse_download($url));
		$counter++;
		print OUT join("\t",@data)."\n";
#die Dumper @data;
		print "Processed: $counter / $total       \r";
	}
	$hs->eof;
	close (IN);
	close (OUT);
	$/ = $orig_sep;
	return \%out_hash;
}



sub parse_download(){
	my $url = shift;
	my $document = `wget -q $url -O /dev/stdout 2>/dev/null`;
	my ($assembly_ftp,$n50);
	my @lines = split("\n",$document);
	for (my $i=0;$i<scalar(@lines);$i++){
#		last if $i < scalar(@lines);
		my $ln = $lines[$i];
		if (!$assembly_ftp && $ln=~/a href="(ftp:\/\/ftp.ncbi.nlm.nih.gov\/genomes\/all\/[^"]+")/){
			$assembly_ftp = dirname($1);
                }elsif ($ln=~/<td>Scaffold N50<\/td>\D+([\d,]+)<\/td>/){
			$n50 = $1;
		}elsif (!$n50 && $ln=~/<td>Contig N50<\/td>\D+([\d,]+)<\/td>/){
                        $n50 = $1;
                }
	}
	$n50=~s/,//g;
	return ($assembly_ftp,$n50);
}

