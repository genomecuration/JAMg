#!/usr/bin/env perl

use strict;
use warnings;



my %combo_types;  # structure: 

## get the percent gene and percent exon accuracy stats:
foreach my $file (<type.*>) {
    my $text = `cat $file`;
    $text =~ m|\#\#\# genes: \(\d+/\d+\) = ([\d\.]+)| or die "Error, cannot extract genes from $text\n";
    my $percent_genes = $1;
    
    $text =~ m|\#\#\# exons: \(\d+/\d+\) = ([\d\.]+)| or die "Error, cannot extract exons from $text\n";
    my $percent_exons = $1;
    
    $file =~ s/type\.//;
    $combo_types{$file}->{genes} = sprintf "%.1f", $percent_genes;
    $combo_types{$file}->{exons} = sprintf "%.1f", $percent_exons;

    #print "$file\t$percent_exons\t$percent_genes\n";
}

###################################
## Get the exon type accuracies:
my @exon_types = qw (single  internal  initial  terminal);
foreach my $file (<*.exon_type_counts.summary>) {

    my %data;

    my @entries = `cat $file`;
    chomp @entries;
    
    foreach my $entry (@entries) {
        my ($exon_type, $gene_type, $count) = split (/\s+/, $entry);
        
        $exon_type = "$exon_type;$gene_type";
        
        $data{$exon_type} = $count;

    }

    $file =~ s/\.exon_type_counts.summary//;

    foreach my $exon_type (@exon_types) {
        
        my $template_count = $data{"$exon_type;template"};
        my $num_correct = $data{"$exon_type;EVM"};

        my $percent = sprintf ("%.1f", $num_correct / $template_count * 100);
        
        #print "$file\t$exon_type\t$percent\n";
        
        $combo_types{$file}->{exon_accuracies}->{$exon_type} = $percent;

    }
    #print "\n";
}

print "#type\texons\tgenes\tsingle\tinternal\tinitial\tterminal\n";

foreach my $type (sort {$combo_types{$a}->{genes}<=>$combo_types{$b}->{genes}} keys %combo_types) {
    
    my $struct = $combo_types{$type};
    print "$type\t" . $struct->{exons} . "\t" . $struct->{genes} .
        "\t" . $struct->{exon_accuracies}->{single} .
        "\t" . $struct->{exon_accuracies}->{internal} .
        "\t" . $struct->{exon_accuracies}->{initial} .
        "\t" . $struct->{exon_accuracies}->{terminal} . "\n";
}


exit(0);


    
                  
