#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use File::Basename;

my $usage = "usage: $0 evaluate.list prediction_file_A [prediction_file_B, ...]\n\n";

my $evaluate_file = $ARGV[0] or die $usage;
shift @ARGV;

my @pred_gffs = @ARGV;
unless (@pred_gffs) { die $usage; }

foreach my $pred_gff (@pred_gffs) {
    $pred_gff = basename($pred_gff);
}

my %exons_info;
my %coords_to_type;
my %genes_to_type;


open (my $fh, $evaluate_file) or die $!;
while (<$fh>) {
    s/\s//g;
    
    my $dir = $_;
    
    my $template = $dir . "/template.gff3";

    my @template_structs = &parse_gene_structs($template);

    if (scalar @template_structs != 1) {
        die "Error, $template does not have exactly one gene!";
    }
    
    my @pred_structs;
    foreach my $pred_gff (@pred_gffs) {
        my $prediction_file = $dir . "/$pred_gff";
        
        push (@pred_structs, &parse_gene_structs($prediction_file));
    }
    
    my $template_gene_struct = shift @template_structs;

    &compare_template_to_others($template_gene_struct, \@pred_structs);
        
}
close $fh;


exit(0);


####
sub parse_gene_structs {
    my $gff3_file = shift;

    my %data;

    open (my $fh, $gff3_file) or die $!;
    
    while (<$fh>) {

        chomp;
        unless (/\w/) { next; }
        my @x = split (/\t/);

        my $feat_type = $x[2];
        my $orient = $x[6];
        
        my ($lend, $rend) = ($x[3], $x[4]);
        my $gene_info = $x[8];
        my $pred_type = $x[1];
        my $contig = $x[0];

        if ($feat_type ne 'CDS') { next; }

        my $gene_id = $pred_type . "/$gene_info";

        my $gene_struct = $data{$gene_id};
        unless (ref $gene_struct) {
            $gene_struct = $data{$gene_id} = { coords => [],
                                               orient => $orient,
                                               pred_type => $pred_type,
                                               gene_id => $gene_id,
                                               contig => $contig,
                                           };
        }
        
        push (@{$gene_struct->{coords}}, [$lend, $rend]);
    }

    close $fh;

    return (values %data);
}

####
sub compare_template_to_others {
    my ($template_gene_struct, $pred_structs_aref) = @_;

    my $contig = $template_gene_struct->{contig};

    my ($min_lend, $max_rend) = &get_min_max_coords($template_gene_struct);

    ## separate each gene prediction by type:
        
    my %exon_coords_to_exon_type;

    foreach my $coord_pair (@{$template_gene_struct->{coords}}) {
        my ($lend, $rend) = @$coord_pair;
        $exon_coords_to_exon_type{"$contig:$lend-$rend"} = "template ";
    }
        
    foreach my $pred (@$pred_structs_aref) {
        
        my $pred_type = $pred->{pred_type};
        
        foreach my $coordset (@{$pred->{coords}}) {
            my ($lend, $rend) = @$coordset;
            if ($lend < $max_rend && $rend > $min_lend) { #overlap
                $exon_coords_to_exon_type{"$contig:$lend-$rend"} .= "$pred_type ";
            }
        }
    }
    
    foreach my $exon_coordname (keys %exon_coords_to_exon_type) {
        print "$exon_coordname\t" . $exon_coords_to_exon_type{$exon_coordname} . "\n";
    }
    

    
    return;


    
}


####
sub get_min_max_coords {
    my ($gene_struct) = @_;

    my $coords_list = $gene_struct->{coords};
    
    my @coords;
    foreach my $coord_pair (@$coords_list) {
        push (@coords, @$coord_pair);
    }

    @coords = sort {$a<=>$b} @coords;

    my $min_lend = shift @coords;
    my $max_rend = pop @coords;

    return ($min_lend, $max_rend);
}


    

        
        
    
