#!/usr/local/bin/perl

use strict;

use CGI;
use CGI::Carp qw(fatalsToBrowser);

my $cgi = new CGI();
my $image = $cgi->param('image');

print $cgi->header(-type=>'image/gif');

open (FILE, $image) or die $!;
while (<FILE>) {
    print;
}
close FILE;

exit(0);
