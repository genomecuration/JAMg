#!/usr/bin/env perl

use strict;
use warnings;

use LWP::UserAgent;
use HTTP::Cookies;
use HTTP::Request;
use File::Basename;
use URI::Escape;

my $usage = "usage: $0 accession_list_file [get|post]\n\n";

my $acc_list_file = $ARGV[0] or die $usage;
my $method = $ARGV[1] || "get";

unless ($method =~ /^(get|post)$/) {
    die "Error, must specify method of retrieval as 'get' or 'post' ";
}

# get accession list
my @accs = `cat $acc_list_file`;
chomp @accs;

unless (@accs) {
    die "Error, no accessions found in file $acc_list_file";
}

my $num_accs = scalar (@accs);
print STDERR "Requesting UIDs for $num_accs accessions\n";

# set up user agent.
my $ua =  LWP::UserAgent->new;
$ua->cookie_jar(HTTP::Cookies->new);
$ua->timeout(1000);


## create URL and parameters:
my $url = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi";
# build params:
my $url_params = "db=nucleotide"
    . "&retmax=$num_accs"
    . "&email=" . "$ENV{USERNAME}" . "@" . "$ENV{HOSTNAME}"
    . "&tool=" . basename($0);

my $search_term = join (" OR ", @accs);
$url_params .= "&term=" . uri_escape($search_term);



## Retrieve data via HTTP

my $response_text = "";

if ($method eq 'get') {
    ## GET method
    my $response = $ua->get($url . "?" . $url_params);
    if ($response->is_success) {
        $response_text = $response->content;
    }
    else {
        print STDERR $response->status_line;
    }
}
else {
    ## POST method
    my $request = new HTTP::Request (POST => $url);
    $request->content_type('application/x-www-form-urlencoded');
    $request->content($url_params);

    my $response = $ua->request($request);
    if ($response->is_success) {
        $response_text = $response->content;
    }
    else {
        print STDERR $response->status_line;
    }
}

while ($response_text =~ m|<Id>(\d+)</Id>|g) {
    print "$1\n";
}


exit(0);



