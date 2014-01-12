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
    #    --database|D          :Entrez database     ie.  "nucleotide"              #
    #    --accessions_file|F   :file containing all accessions to retrieve.        #
    #                                                                              #
    #  optional:                                                                   #
    #                                                                              #
    #  --chunkSize|C   :number of records to retrieve per call (default 500)       #
    #  --DEBUG|d       :debug mode, verbose                                        #
    #                                                                              #
    ################################################################################


};

my ($filename, $db);

our $DEBUG = 0;

my $chunk_size = 500;
&GetOptions ('accessions_file|F=s' => \$filename,
             'db|D=s' => \$db,
             'chunk_size|C=i' => \$chunk_size,
             'DEBUG|d' => \$DEBUG,
             );

unless ($filename && $db) {
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


open (my $fh, $filename) or die "Error, cannot open $filename";

my $num_per_query = 1000;
my $count = 0;
my $tmpfile = "$$.epost";
open (my $ofh, ">$tmpfile") or die $!;
while (<$fh>) {
    $count++;
    print $ofh $_;
    if ($count % $num_per_query == 0) {
        close $ofh;
        &run_epost($tmpfile);
        open ($ofh, ">$tmpfile") or die $!;
        $count = 0;
    }
}
if ($count) {
    &run_epost($tmpfile);
}
unlink $tmpfile;

exit(0);

####
sub run_epost {
    my $idfile = shift;
    
    print STDERR "processing $idfile\n";
    
    system "cat $idfile";
        
    my %params = ( db => $db,
                   id => $idfile,
                   email => "$ENV{USER}\@tigr.org",
                   tool => "entrez_efetch.pl",
                   # usehistory => 'y',
                   );
    
    my %results = epost_file (%params);
    my %summary = esummary( db=>$db,
                            query_key => $results{query_key},
                            WebEnv => $results{WebEnv},
                            #retmode => 'text',
                            tool => $params{tool},
                            email => $params{email},
                            );
    
    foreach my $uid (keys %summary) {
        print "$uid\n";
        my $data_href = $summary{$uid};
        foreach my $key (keys %$data_href) {
            my $val = $data_href->{$key};
            print "\t$key\t$val\n";
        }
    }
    
   
}

