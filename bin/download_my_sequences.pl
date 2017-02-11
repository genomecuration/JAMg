#!/usr/bin/env perl

=pod

=head1 NAME

download_my_sequences.pl - Download sequences from NCBI

=head1 VERSION

 Version 0.0001

=head1 USAGE 

Download sequences from NCBI

        -o|outfile  :s  => Outfile filename (def. ncbitax or given ''taxonomy'')
        -f|format   :s  => Format of outfile, ''fasta'' or ''gi'' for accessions only (def. fasta)
        -m|molecule :s  => ''protein'' (UNIPROT ); ''nucleotide'' (EMBL); ''cDNA'' (EMBL) or ''GSS'' (EMBL)
        -s|species  :s  => Latin name or NCBI TaxID of species
        -t|taxonomy :s  => A taxonomy if searching proteins
        -q|query    :s  => A specific query in NCBI format, e.g. 'chimpanzee[orgn]+AND+biomol+mrna[prop]'

=head1 DESCRIPTION

 Examples

        download_my_sequences.pl -s 'Drosophila melanogaster' -m cDNA -f fasta
        download_my_sequences.pl -t Insecta -t fasta -m protein   => downloads Insecta division of Uniprot
        download_my_sequences.pl -q '([uniprot-Taxonomy:Insecta]&[uniprot-Organelle:mitochondrion])' -t fasta     => As above but downloads mitochondrial sequences only
        Or more give exact complicated query: download invertebrate rDNA genes which are complete, avoid electronic identified sequences
         download_my_sequences.pl -q '(([embl-Division:INV]&(([embl-Description:rDNA]|[embl-Description:rRNA])&[embl-Description:gene]))&([embl-Description:complete]![embl-Description:similar]))'

=head1 AUTHORS

 Alexie Papanicolaou 1

        1 CSIRO Ecosystem Sciences, Australia
        alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

This software is released under the GNU General Public License version 3 (GPLv3).
It is provided "as is" without warranty of any kind. You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html. Please note that incorporating the whole software or parts of$

=head1 BUGS & LIMITATIONS

None known

=cut

use strict;
use warnings;
use LWP::Simple;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";

my ( $query, $query_text, $molecule,$species_latin,$ncbi_taxid,$taxonomy,$outfile );
my $format     = "fasta";
my $debug;
pod2usage $! unless &GetOptions(
	'outfile:s'  => \$outfile,
	'format:s'  => \$format,
	'query:s'    => \$query,
	'molecule:s' => \$molecule,
	'species:s'	=> \$species_latin,
	'taxonomy:s'	=> \$taxonomy,
	'verbose|debug' => \$debug,
);

# Sanity checks
if ( !$query && ((!$species_latin || !$taxonomy) && !$molecule)) { pod2usage; }
elsif ($taxonomy && $molecule!~/protein/i){die ("Taxonomy search is only possible with protein molecules...\n");}
elsif ($species_latin){
	$ncbi_taxid=&get_ncbitaxid($species_latin);
	if (!$ncbi_taxid){die ("Failed to find NCBI taxid for $species_latin\n");}
}
if (!$outfile && $taxonomy){$outfile="$taxonomy.$molecule.$format";}
elsif (!$outfile && $ncbi_taxid){$outfile="$ncbi_taxid.$molecule.$format";}

my $service_request= 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?';
my $database;

if ($molecule && $molecule =~ /protein/i){
	$database = 'protein';
}elsif ($molecule=~/cdna/i){
	$database = 'nucest';
}elsif (!$molecule || $molecule =~ /nucleotide/i){
	$database = 'nuccore';
}else{
	 die "I don't know what to do!\n";
}

$service_request .= "db=$database&usehistory=y&term=";


if ($query_text){
	print "Using user-provided query\n$query_text\n";
}elsif (!$query_text && $ncbi_taxid ) {
	$query_text = 'txid'.$ncbi_taxid.'[Organism:noexp]';
} elsif (!$query_text && $taxonomy) {
	$query_text = $taxonomy.'[porgn]';
} else {
	die "I don't know what to do!\n";
}
print "Making this query: $query_text from NCBI's $database\n";
$service_request.=$query_text;
print "Fetching $service_request\n" if $debug;
#http://www.ncbi.nlm.nih.gov/books/NBK25498/#chapter3.Application_3_Retrieving_large
#post the esearch URL
my $output = get($service_request);

#parse WebEnv, QueryKey and Count (# records retrieved)
my $web = $1 if ($output =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/);
my $count = $1 if ($output =~ /<Count>(\d+)<\/Count>/);

die "No results found\n" unless $count && $count >0;
print "Will download $count results into $outfile\n";

#open output file for writing
open(OUT, ">$outfile") || die "Can't write to file $outfile!\n";


#retrieve data in batches of 500
my $retmax = 500;
my $efetch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=$database&WebEnv=$web&query_key=$key&retmax=$retmax&rettype=$format&retmode=text&";
print "Fetching using efetch: $efetch_url\n" if $debug;

$|=1;
for (my $retstart = int(0); $retstart < $count; $retstart += $retmax) {
	print "Retrieved $retstart / $count     \r";
        my $fetch = $efetch_url . "retstart=$retstart";
        my $efetch_out = get($fetch);
        print OUT $efetch_out if $efetch_out;
	print "Retrieved $retstart / $count     \r";
}
close OUT;



####################
sub get_ncbitaxid ($){
        my $species=shift;
        my $ncbi_taxid;
        if ($species=~/^\d+$/){$ncbi_taxid=$species;}
        elsif ( $species =~ /^([A-Za-z]{2})[a-z]+\s/ ) {
		my $url = "https://www.ncbi.nlm.nih.gov/taxonomy/?term=\%22$species\%22\&report=taxid&format=text";
		print "Fetching $url\n" if $debug;
                system("wget -q \'$url' -O species.query"
                );
                open( SPECIES, "species.query" );
                while ( my $line = <SPECIES> ) {
                        if ( $line =~ /^<pre>(\d+)/ ) {
                                $ncbi_taxid = $1;
                        }
                }
                close(SPECIES);
                unlink("species.query") unless $debug;
        }
        return $ncbi_taxid;
} 

