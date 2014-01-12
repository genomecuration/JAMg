#!/usr/local/bin/perl

use strict;
use warnings;

use LWP::UserAgent;
use HTTP::Cookies;
use HTTP::Request;
use File::Basename;
use URI::Escape;


my $usage = "usage: $0 UID_list_file\n\n";

my $uid_file = $ARGV[0] or die $usage;

##########################################################
## create a UID string, comma-delimited
my @uids = `cat $uid_file`;
chomp @uids;
my $uid_listing = join (',', @uids);
my $num_uids = scalar (@uids);

###########################################################
## post the UID list to the server:
my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/epost.fcgi";
my $url_params = "db=nucleotide&id=$uid_listing";

print "// Posting: " . scalar(@uids) . " uids.\n\n";

my $ua = new LWP::UserAgent;
$ua->agent("epost_file/1.0 " . $ua->agent);
my $req = new HTTP::Request POST => "$url";
$req->content_type('application/x-www-form-urlencoded');
$req->content("$url_params");
my $response = $ua->request($req);


my $raw = $response->content;
print $raw;

# get the web history key
$raw  =~ /<QueryKey>(\d+)<\/QueryKey>.*<WebEnv>(\S+)<\/WebEnv>/s;

my $query_key = $1 or die "$raw\n error, no query_key";
my $WebEnv = $2 or die "$raw\n error, no webenv";


#########################################################
# retrieve the list back from the server.
{
    
    print "\n\n// Getting the count of posted UIDs.\n\n";
    
    my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi";
    my $url_params = "db=nucleotide&"
        . "usehistory=y&"
        . "WebEnv=$WebEnv&"
        . "term=%23$query_key"; # note this use of the query_key looks like a hack...
    
    my $url_string = "$url?$url_params";
    

    print "URL: $url_string\n";
    
    my $response = $ua->get($url_string);
    print $response->content;
}


########################################################
# get the sequences
{
    
    print "\n\n// Retrieving the sequences\n\n";
    
    my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi";
    my $url_params = "db=nucleotide&"
        . "usehistory=y&"
        . "retmode=text&rettype=fasta&"
        . "retstart=0&"
        . "retmax=$num_uids&"
        . "query_key=$query_key&"
        . "WebEnv=$WebEnv";
    
    my $url_string = "$url?$url_params";
    
    print "$url_string\n";
    
    my $response = $ua->get($url_string);
    
    print $response->content;
}


exit(0);






