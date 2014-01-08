#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;

my %EV_TYPES;
my @entries; ## list of (ev_type_to_weights_href, accuracy_values_href), ...

while (<STDIN>) {
    if (/^Scoring weights:/) {
        my %ev_type_to_weights;
        
        ## parse weight values:
        my $line = <STDIN>;
        while ($line !~ /--------/) {
            $line =~ s/^\s+//;
            chomp $line;
            my ($ev_type, $weight) = split (/\s+/, $line);
            $ev_type_to_weights{$ev_type} = $weight;
            $EV_TYPES{$ev_type} = 1;
			$line = <STDIN>;
        }
        
        ## get accuracy values:
        until ($line =~ /^\#\#\# /) {
            $line = <STDIN>;
        }
        my %accuracy_values;
        while ($line =~ /^\#\#\# /) {
            chomp $line;
            my ($trash, $token, $counts, $eqsign, $value) = split (/\s+/, $line);
            if ($token =~ /genes|exons/) {
                $accuracy_values{$token} = $value;
            }
            $line = <STDIN>;
        }
        
        push (@entries, [\%ev_type_to_weights, \%accuracy_values] );
    }
}

my @ev_types = keys %EV_TYPES;
my @accuracy_types = qw (genes: exons:);
foreach my $entry (@entries) {
    my ($ev_types_to_weights_href, $accuracy_values_href) = @$entry;
    
    my $output_text = "";
    foreach my $ev_type (@ev_types) {
        my $weight = $ev_types_to_weights_href->{$ev_type} || 0;
        $output_text .= "$ev_type\t$weight\t";
    }
    foreach my $accuracy_type (@accuracy_types) {
        my $accuracy_value = $accuracy_values_href->{$accuracy_type};
		unless (defined $accuracy_value) { die "Error, no accuracy value for $accuracy_type of " . Dumper ($entry); }
        $output_text .= "$accuracy_type\t$accuracy_value\t";
    }
	
    chomp $output_text;
    print "$output_text\n";
}


exit(0);

