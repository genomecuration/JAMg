#!/usr/bin/env perl

use strict;
use warnings;

my $debug = 1;

my $species_latin = shift;
die unless $species_latin;

my $ncbi_taxid=&get_ncbitaxid($species_latin);
if (!$ncbi_taxid){die ("Failed to find NCBI taxid for $species_latin\n");}
print "$species_latin: NCBI TaxID $ncbi_taxid\n";

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

