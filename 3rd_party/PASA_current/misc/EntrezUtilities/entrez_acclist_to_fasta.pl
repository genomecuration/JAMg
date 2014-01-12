#!/usr/local/bin/perl

use strict;
use warnings;
use lib ($ENV{EUK_MODULES});
use Eutilities_ncbi;
use Data::Dumper;
use URI::Escape;
use Fasta_reader;
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
    #  -n              :number of accessions to search per UID request.  (default 100)
    ################################################################################


};

my ($filename, $db);

our $DEBUG = 0;

# open (STDERR, ">&STDOUT");

my $num_accs_per_uid_request = 100;

my $chunk_size = 500;
&GetOptions ('accessions_file|F=s' => \$filename,
             'db|D=s' => \$db,
             'chunk_size|C=i' => \$chunk_size,
             'DEBUG|d' => \$DEBUG,
             'n=i' => \$num_accs_per_uid_request,
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
my @accessions = `cat $filename`;
chomp @accessions;

## convert the accessions into UIDs


my $num_accs_total = scalar (@accessions);
my %accessions_fetching;
foreach my $accession (@accessions) {
    $accessions_fetching{$accession} = 1;
}


my @UIDs;
for (my $i = 0; $i <= $#accessions; $i += $num_accs_per_uid_request) {
    my @acc_list;
    for my $j ($i..($i+$num_accs_per_uid_request-1)) {
        my $acc = $accessions[$j];
        push (@acc_list, $acc) if ($j <= $#accessions && $acc =~ /\w/);
    }
    
    my $num_accs_to_query = scalar (@acc_list);

    if ($num_accs_to_query < 1) { next; }
    
    my $term = join ("[accession] OR ", @acc_list);
    $term .= "[accession]";
    
    if ($num_accs_per_uid_request <= 10) {
        print STDERR "searching: $term\n";
    }

    $term = uri_escape($term);
    
    my %params = ( db => $db,
                   term => $term,
                   email => "$ENV{USER}\@tigr.org",
                   tool => "entrez_efetch.pl",
                   usehistory => 'y',
                   retmax => 10 * $num_accs_per_uid_request, # beware, any hard limit is probably a bad idea...
                   );

    print STDERR "* retrieving $num_accs_per_uid_request from $i of $num_accs_total\n";
    
    my %results = esearch (%params);
        
    my $uids_aref = $results{uids};
    if (ref $uids_aref eq "ARRAY") {
        my @retrieved_uids = @$uids_aref;
        my $num_retrieved = scalar (@retrieved_uids);
        print STDERR "-retrieved $num_retrieved UIDs for $num_accs_to_query accessions\n"; 
        push (@UIDs, @$uids_aref);
    }
    else {
        print STDERR "-Sorry, no UIDs retrieved.\n";
    }
}

my $num_per_query = 1000;
my $count = 0;
my $tmpfile = "$$.epost";
open (my $ofh, ">$tmpfile") or die $!;

my $num_UIDs_total = scalar (@UIDs);
  
print STDERR "*** Received $num_UIDs_total UIDs for $num_accs_total\n";

## write them to a file for safe keeping:
{

    my $outputfile = "$filename.uids";
    open (my $fh, ">$outputfile") or die $!;
    foreach my $uid (@UIDs) {
        print $fh $uid . "\n";
    }
    close $fh;
    print STDERR "-UIDs archived to file: $outputfile\n";

}

my @uids_to_search;
foreach my $uid (@UIDs) {

    $count++;
    #print STDERR $count . "\n";
    push (@uids_to_search, $uid);
    if ($count % $num_per_query == 0) {
        &get_sequences(@uids_to_search);
        $count = 0;
        @uids_to_search = (); # clear it.
    }
}
if (@uids_to_search) {
    &get_sequences(@uids_to_search);
}

if (%accessions_fetching) {
    print STDERR "*WARNING* Not all accessions were retrieved from genbank.  \n";
    my $accs_file = "missing_accessions.$$.list";
    open (my $accs_fh, ">$accs_file") or die "Error, cannot write to $accs_file";
    print $accs_fh join ("\n", keys %accessions_fetching) . "\n";
    close $accs_fh;
    print STDERR "See missing accessions in file: $accs_file\n";
    exit(1);
}

exit(0);


####
sub get_sequences {
    my @uids = @_;

    print STDERR "get_sequences() . " . scalar (@uids) . " uids.\n";
    
    my %uids_needed;
    foreach my $uid (@uids) {
        $uids_needed{$uid} = 1;
    }
    
    while (%uids_needed) {
        
        my $uid_file = "tmp.$$.file";
        open (my $fh, ">$uid_file") or die "Error, cannot create file $uid_file";
        
        my @uids_here = keys %uids_needed;
        
        foreach my $uid (@uids_here) {
            print $fh "$uid\n";
        }
        close $fh;
        
        print STDERR "Attempting sequence retrieval for " . scalar (@uids_here) . " UIDs.\n";
        
        my $fasta = &run_epost($uid_file);
        $fasta =~ s/Error: The resource is temporarily unavailable//g;
        
        my $tmp_fasta_file = "tmp.$$.fasta";
        open (my $fasta_fh, ">$tmp_fasta_file") or die "Error, cannot open file $tmp_fasta_file for writing";
        print $fasta_fh $fasta;
        close $fasta_fh;

        my $fasta_reader = new Fasta_reader($tmp_fasta_file);
        while (my $seqobj = $fasta_reader->next() ) {
            my $fasta_format = $seqobj->get_FASTA_format();
            
            my $accession = $seqobj->get_accession();
            if ($accession =~ /^gi\|(\d+)/) {
                my $uid = $1;
                delete $uids_needed{$uid};
                
                ## retrieve the accession:
                my @acc_parts = split (/\|/, $accession);
                my $core_acc = $acc_parts[3];
                $core_acc =~ s/\.\d+$//;
                if ($accessions_fetching{$core_acc}) {
                    print $fasta_format;
                    delete $accessions_fetching{$core_acc};
                }
            }
        }
        if (%uids_needed) {
            my @remaining = keys %uids_needed;
            my $num_remaining = scalar (@remaining);
            print STDERR "Error, didn't receive response for $num_remaining uids: @remaining\ntrying again.\n";
        }
    }
    
    print STDERR "** found all UIDs in fasta response.\n";
    
    return;
    

}


####
sub run_epost {
    my $idfile = shift;
    
    # print STDERR "processing $idfile\n";
    
    my %params = ( db => $db,
                   id => $idfile,
                   email => "$ENV{USER}\@tigr.org",
                   tool => "entrez_efetch.pl",
                   usehistory => 'y',
                   );
    
    my %results = epost_file (%params);
        

    # efetch_batch
    
    my $results = efetch_batch  ( db=>$db,
                                  query_key => $results{query_key},
                                  WebEnv => $results{WebEnv},
                                  retmode => $REPORT_STYLES{$db}->{retmode},
                                  rettype => $REPORT_STYLES{$db}->{rettype},
                                  retstart=> 1,
                                  retmax => $chunk_size,
                                  tool => $params{tool},
                                  email => $params{email},
                                  returnOutput => 1
                                  );
    
    

        
    return $results;
}

