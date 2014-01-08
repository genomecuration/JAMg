#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;

my %ev_types;
my %acc_tokens;
my @rounds;


my $mode;
while (<STDIN>) {
    s/\'//g;
    if (/^MODE: (\S+)/) {
        $mode = $1;
    }
    elsif (/=>\d/) {
        s/\s+//g;
        my @entries = split (/,/);
        my %ev_and_weights;
        foreach my $entry (@entries) {
            my ($ev_type, $weight) = split (/=>/, $entry);
            $ev_types{$ev_type} = 1;
            $ev_and_weights{$ev_type} = $weight;
        }
        

        my $accuracies = <STDIN>;
        unless ($accuracies =~ /exon_acc/) { next; }
        my %acc_vals;
        $accuracies =~ s/\s//g;
        my @stats = split (/,/, $accuracies);
        foreach my $stat (@stats) {
            my ($accuracy_type, $accuracy_value) = split (/=>/, $stat);
            $acc_vals{$accuracy_type} = $accuracy_value;
            $acc_tokens{$accuracy_type} = 1;
        }

        push (@rounds, { mode => $mode,
                         weights => \%ev_and_weights,
                         accuracy => \%acc_vals, });
        
    }
}

my @ev_type_list = sort keys %ev_types;
my @accuracy_tokens = qw (exon_accuracy gene_accuracy BJHscore);

foreach my $round (@rounds) {
    
    my $mode = $round->{mode};
    my $weights = $round->{weights};
    my $accuracy = $round->{accuracy};
    # print Dumper ($round);

    my $weight_combo = join (",", sort keys %$weights);
    
    my $weights_text = "";
    
    foreach my $ev_type (@ev_type_list) {
        my $weight = $weights->{$ev_type};
        unless (defined $weight) { $weight = ""; }

        $weights_text .= "\t$ev_type\t$weight";
    }
    
    my $accuracy_text = "";
    foreach my $acc_token (@accuracy_tokens) {
        my $acc_value = $accuracy->{$acc_token};
        $accuracy_text .= "\t$acc_token\t$acc_value";
    }

    my $score = $accuracy->{BJHscore};
    
    print "$mode\t$score\tWEIGHTS:$weight_combo$weights_text\tACCURACY:$accuracy_text\n";
}


exit(0);



