#!/usr/bin/env perl

use strict;
use warnings;
use GD::Graph::bars;
use CGI;
use Cwd;

umask 0000;

my $usage = "usage: $0 weight_and_accuracy_summary.log\n\n";

my $weight_summary = $ARGV[0] or die $usage;

my $html_file = "$weight_summary.html";
my $images_dir = "images";
unless (-d $images_dir) {
    mkdir $images_dir;
}

my $image_counter = 0;
my @experiments;

open (my $fh, $weight_summary) or die $!;
while (<$fh>) {
    chomp;
    my ($mode, $evidence_and_weights, $accuracies) = split (/\t/);
    
    my %weights;
    my @combos = split (/,/, $evidence_and_weights);
    foreach my $combo (@combos) {
        my ($ev_type, $weight) = split (/=>/, $combo);
        $weights{$ev_type} = $weight;
    }

    my %accuracy_vals;
    foreach my $accuracy_type (split (/,/, $accuracies)) {
        my ($attribute_measured, $accuracy_value) = split (/=>/, $accuracy_type);
        $accuracy_vals{$attribute_measured} = $accuracy_value;
    }

    ## create graph:
    my $graph = GD::Graph::bars->new();
    my @ev_types = sort keys %weights;
    my @weights;
    foreach my $ev_type (@ev_types) {
        my $weight = $weights{$ev_type};
        push (@weights, $weight);
    }

    my $data = GD::Graph::Data->new([ \@ev_types, \@weights ]);
    
    $graph->set( 
                 x_label         => 'Weight',
                 y_label         => 'Evidence Types',
                 #title           => 'A Simple Bar Chart',
                 #y_max_value     => 8,
                 #y_tick_number   => 8,
                 #y_label_skip    => 2,
                 
                 x_labels_vertical => 1,
                 
                 # shadows
                 #bar_spacing     => 8,
                 #shadow_depth    => 4,
                 #shadowclr       => 'dred',
                 
                 #transparent     => 0,
                 );

    $graph->plot($data);
    
    $image_counter++;
    my $image_filename = "$images_dir/img.$$.$image_counter";
    
    open (my $imgfh, ">$image_filename") or die $!;
    binmode $imgfh;
    print $imgfh $graph->gd()->png();
    close $imgfh;
    
    push (@experiments, { image_filename => $image_filename,
                          weights => \%weights,
                          accuracies => \%accuracy_vals,
                      }
          );
    
}
close $fh;


@experiments = reverse sort { $a->{accuracies}->{BJHscore} <=> $b->{accuracies}->{BJHscore} } @experiments;

open ($fh, ">$html_file") or die $!;
my $cgi = new CGI();
print $fh $cgi->start_html("-title" => "EVM accuracy given weight combinations");

print $fh "<table>\n";

my $counter = 0;
foreach my $experiment (@experiments) {
    
    $counter++;
    
    my $accuracy_href = $experiment->{accuracies};
    my $image_filename = $experiment->{image_filename};
    my $weights_href = $experiment->{weights};

    print $fh "<tr><td>$counter</td><td><img src=\"$image_filename\" alt=image_$image_filename></td>\n";
    ## enumerate the weight values:
    my $weights_table_text = "";
    foreach my $ev_type (keys %$weights_href) {
        $weights_table_text .= "$ev_type = " . $weights_href->{$ev_type} . "<br>";
    }
    print $fh "<td><font size=-2>$weights_table_text</font></td>\n";

    my $accuracy_table = "";
    foreach my $accuracy_type (keys %$accuracy_href) {
        $accuracy_table .= "$accuracy_type = " . $accuracy_href->{$accuracy_type} . "<br>";
    }

    print $fh "<td>$accuracy_table</td></tr>\n";
    print $fh "<tr><td colspan=3>&nbsp;</td></tr>\n";
}

print $fh "</table>\n" . $cgi->end_html();


## launch the browser:
my $cmd = "firefox file://" . cwd() . "/$html_file";
system $cmd;



exit(0);

        

