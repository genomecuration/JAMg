#!/usr/local/bin/perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES});
use Eutilities_ncbi;

use Getopt::Long qw(:config no_ignore_case bundling);



my $usage = qq {

    ################################################################################
    #   Entrez Efetch                                                              #
    #                                                                              #
    #  required:                                                                   #
    #    --database|D    :Entrez database     ie.  "Pubmed"                        #
    #    --query|q       :Entrez query        ie.  "haas bj[auth]"                 #
    #                                                                              #
    #  optional:                                                                   #
    #                                                                              #
    #  --chunkSize|C   :number of records to retrieve per call (default 500)       #
    #  --DEBUG|d       :debug mode, verbose                                        #
    #                                                                              #
    ################################################################################


};

my ($query, $db);

our $DEBUG = 0;

my $chunk_size = 500;
&GetOptions ('query|q=s' => \$query,
             'db|D=s' => \$db,
             'chunk_size|C=i' => \$chunk_size,
             'DEBUG|d' => \$DEBUG,
             );

unless ($query && $db) {
    die $usage;
}

$db = lc $db;



my %REPORT_STYLES = ( pubmed => { 
                                    retmode => "text",
                                    rettype => "abstract",
                                },
                      
                      nucleotide => {
                                      retmode => "text",
                                      rettype => "fasta",
                                  },

                      protein => {
                                        retmode => "text",
                                        rettype => "fasta",
                                    },

                      # ... add others as we need them.
                      #     could also make this command-line configurable.
                      
                      
                      );



my %params = ( db => $db,
               term => $query,
               email => "$ENV{USER}\@tigr.org",
               tool => "entrez_efetch.pl",
               usehistory => 'y',
               );

my %results = esearch (%params);


efetch_batch ( db=>$db,
               query_key => $results{query_key},
               WebEnv => $results{WebEnv},
               retmode => $REPORT_STYLES{$db}->{retmode},
               rettype => $REPORT_STYLES{$db}->{rettype},
               retmax => $chunk_size,
               tool => $params{tool},
               email => $params{email},
               );



exit(0);

