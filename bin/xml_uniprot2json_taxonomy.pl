#!/usr/bin/env perl

=pod

=head1 NAME

 xml_uniprot2json_taxonomy.pl

=head1 USAGE

	-xml|in   	   :s	Input XML file from Uniprot
	-stop|max|limit	   :i 	Max number of entries to get (def to 0 which is all)
	-taxonomy          :s	Fetch only those having this matching Taxon string
	-pretty|human		JSON output is human readable


=head1 AUTHORS

 Alexie Papanicolaou 1
	alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

This software is released under the GNU General Public License version 3 (GPLv3).
It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html.
Please note that incorporating the whole software or parts of its code in proprietary software
is prohibited under the current license.

=head1 BUGS & LIMITATIONS

None known so far.

=cut

use strict;
use warnings;
use XML::LibXML::Reader;
use XML::Simple;
use JSON;
#use Data::Dumper;
use Getopt::Long;
use Pod::Usage;

my ($input_xml,$taxonomy_filter,$make_pretty,$json_obj);
my $max_limit = int(0);

GetOptions(
	'in|xml=s'	    	=> \$input_xml,
	'stop|max|limit:i'	=>\$max_limit,
	'taxonomy:s' 		=> \$taxonomy_filter,
	'pretty'		=> \$make_pretty
);
pod2usage "No input XML file!\n" if !$input_xml;

$input_xml = shift if !$input_xml;

my $outfile = $input_xml;
$outfile .= $taxonomy_filter if $taxonomy_filter;
$outfile .= $max_limit if $max_limit;
$outfile .= '.json';

my $reader = XML::LibXML::Reader->new(location => $input_xml) or die "Cannot read $input_xml: $!";
$json_obj = $make_pretty ? JSON::XS->new->canonical->utf8->pretty : JSON::XS->new->canonical->utf8;

open my $fh, '>', $outfile or die "Could not open file $outfile: $!";
print $fh "{\n";

my $entry_count = 0;

while ($reader->read) {
    if ($reader->nodeType == XML_READER_TYPE_ELEMENT) {
        if ($reader->name eq 'entry') {
            my $has_taxonomy = 0;
            my $xml_text = $reader->readOuterXml();
            my $xml_obj = XMLin($xml_text, ForceArray => 1,KeyAttr=>[]);
#  The default value for 'KeyAttr' is ['name', 'key', 'id'].
	    my $taxa = $xml_obj->{organism}->[0]->{lineage}->[0]->{taxon};
	    foreach my $t (@{$taxa}){
		$has_taxonomy = 1 if ($taxonomy_filter && $t eq $taxonomy_filter);
	    }
	    next unless $has_taxonomy;
	    my $acc= $xml_obj->{accession};
               $acc = @{$acc}[0];
	    print $fh ",\n" if $entry_count > 0;
	    print $fh '"' . $acc . '": ' . $json_obj->encode($xml_obj);
	    $entry_count++;
	}
    }
   last if $max_limit && $entry_count == $max_limit;
}

print $fh "\n}";
close $fh;


