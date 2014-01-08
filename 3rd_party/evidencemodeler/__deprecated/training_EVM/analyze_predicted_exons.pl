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


print "-writing summary files: exons.summary and genes.summary\n\n";
{
    open (my $fh, ">exon_coords_to_type.summary") or die $!;
    foreach my $coords_dat (keys %coords_to_type) {
        print $fh $coords_dat . "\t" . $coords_to_type{$coords_dat} . "\n";
    }
    close $fh;
}

{
    open (my $fh, ">exon_type_counts.summary") or die $!;
    foreach my $exon_type (keys %exons_info) {
        my $preds_href = $exons_info{$exon_type};
        foreach my $pred_type (keys %$preds_href) {
            print $fh "$exon_type\t$pred_type\t" . $preds_href->{$pred_type} . "\n";
        }
    }
    close $fh;
}

{
    open (my $fh, ">genes_to_type.summary") or die $!;
    foreach my $gene_contig_id (keys %genes_to_type) {
        print $fh "$gene_contig_id\t" . $genes_to_type{$gene_contig_id} . "\n";
    }
    close $fh;
}


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

    $contig .= ",$min_lend-$max_rend";

    ## separate each gene prediction by type:
    
    my %type_to_pred_list;
    foreach my $pred (@$pred_structs_aref) {
        
        my $pred_type = $pred->{pred_type};
        
        my $pred_list = $type_to_pred_list{$pred_type};
        unless (ref $pred_list) {
            $pred_list = $type_to_pred_list{$pred_type} = [];
        }

        push (@$pred_list, $pred);
    }

    ## assign each template exon a type:
    my @cds_coords = @{$template_gene_struct->{coords}};
    my $orient = $template_gene_struct->{orient};

    my @cds_structs;

    if (scalar @cds_coords == 1) {
        my $coordset = shift @cds_coords;
        my ($lend, $rend) = @$coordset;
        push (@cds_structs, { lend => $lend, rend => $rend, type => 'single' });
    }
    else {
        ## separate to get initial, terminal, and internals:
        @cds_coords = sort {$a->[0]<=>$b->[0]} @cds_coords;
        if ($orient eq '-') {
            @cds_coords = reverse @cds_coords;
        }
        
        my $initial_cds_coords = shift @cds_coords;
        push (@cds_structs, { lend => $initial_cds_coords->[0], rend => $initial_cds_coords->[1], type => "initial" } );
        
        my $terminal_cds_coords = pop @cds_coords;
        push (@cds_structs, { lend => $terminal_cds_coords->[0], rend => $terminal_cds_coords->[1], type => "terminal" } );
        
        foreach my $remaining_coordset (@cds_coords) {
            push (@cds_structs, { lend => $remaining_coordset->[0], rend => $remaining_coordset->[1], type => "internal" } );
        }

    }

    
    ## check exon accuracy

    foreach my $cds_struct (@cds_structs) {
        my ($lend, $rend, $type) = ($cds_struct->{lend}, $cds_struct->{rend}, $cds_struct->{type});
        $coords_to_type{"$type;$contig;$lend;$rend"} = "template ";
        $exons_info{$type}->{"template"}++;
        
        foreach my $pred_type (sort keys %type_to_pred_list) {
            
            my $pred_list = $type_to_pred_list{$pred_type};
            
          PRED_STRUCTS:
            foreach my $pred_struct (@$pred_list) {
                
                my $coords_aref = $pred_struct->{coords};

                foreach my $coordset (@$coords_aref) {
                    my ($coords_lend, $coords_rend) = @$coordset;
                    
                    if ($coords_lend == $lend && $coords_rend == $rend) {

                        $coords_to_type{"$type;$contig;$lend;$rend"} .= "$pred_type ";
                        $exons_info{$type}->{$pred_type}++;
                        last PRED_STRUCTS;

                    }
                }
            }
        }
    }

    ## check gene accuracy:
    
    my $template_coords = $template_gene_struct->{coords};
    $genes_to_type{$contig} = "template ";
    foreach my $pred_type (sort keys %type_to_pred_list) {
            
        my $pred_list = $type_to_pred_list{$pred_type};
        
      PRED_STRUCTS:
            foreach my $pred_struct (@$pred_list) {
                my $coords_aref = $pred_struct->{coords};
                if (&same_gene_structure($template_coords, $coords_aref)) {
                    $genes_to_type{$contig} .= "$pred_type ";
                    last PRED_STRUCTS;
                }
            }
    }

    return;

}
    
####
sub same_gene_structure {
    my ($coords_A_aref, $coords_B_aref) = @_;

    my %data;
    foreach my $coordset (@$coords_A_aref) {
        my ($lend, $rend) = @$coordset;
        $data{"$lend-$rend"} = 1;
    }

    foreach my $coordset (@$coords_B_aref) {
        my ($lend, $rend) = @$coordset;
        if (exists $data{"$lend-$rend"} ) {
            delete $data{"$lend-$rend"};
        }
        else {
            return (0); # not same
        }
    }
    
    if (%data) {
        return (0); # not same
    }
    else {
        return (1); # same
    }
        
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


    

        
        
    
