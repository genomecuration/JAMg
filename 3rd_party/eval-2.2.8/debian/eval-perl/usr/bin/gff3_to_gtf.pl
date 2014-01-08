#! /usr/bin/perl 

use Bio::DB::GFF;
use GTF;
use Carp;
use strict;

my $usage = "$0 gff3-file\n";
@ARGV == 1 or die $usage;

my $gff = Bio::DB::GFF->new( -adaptor => "memory",
                             -gff    => $ARGV[0]);

my $gtf = GTF::new();

for my $gff_gene ($gff->features("gene")) {
    my $strand;
    if($gff_gene->strand == 1) { $strand = '+' }
    else                       { $strand = '-' }
    my $gtf_gene = GTF::Gene::new(find_id($gff_gene),
        $gff_gene->ref, $gff_gene->source, $strand);
    for my $gff_tx ($gff_gene->features("mRNA")) {
        my $gtf_tx = GTF::Transcript::new(find_id($gff_tx));
        for my $gff_feat 
                ( $gff_gene->features("exon"), $gff_gene->features("CDS")) {
            $gff_feat->absolute(1);
            my $phase = $gff_feat->phase || "0";
            my $score = $gff_feat->score || ".";
            my $gtf_feat = GTF::Feature::new($gff_feat->method,
               $gff_feat->start, $gff_feat->end, $score, $phase);
            $gtf_tx->add_feature($gtf_feat);
        }
        $gtf_gene->add_transcript($gtf_tx);
        $gtf_tx->correct_gff_stop_codon;
        $gtf_tx->fix_start_codon(1);
        $gtf_tx->fix_stop_codon(1);
        $gtf_tx->create_utr_objects_from_exons;
        $gtf_tx->fix_utr3_stop_codon_overlap;
        $gtf_tx->fix_frame;
        $gtf_tx->{Exons} = [];
    }
    $gtf->add_gene($gtf_gene);
}

$gtf->output_gtf_file(\*STDOUT);

sub find_id {
    my ($feat) = @_;

    #gff3 attempt first
    my %attrs3 = $feat->attributes;
    my $id = "";
    for my $k ("ID", "Name", "Parent") {
        if(defined ($attrs3{$k})) {
            $id = $attrs3{$k} ;
            last;
        }
    }

    
    # sometimes it does not come in the above form, but instead as Note attrs.
    # i.e. Note => ID=..., Note => Parent=...
    # this is generall for gff2-formatted files.
    # currently we just hope for an ID to always be there.
    my @attrs = $feat->attributes('Note');
    if($id !~ m/\S/) {
        for my $val (@attrs) {
            if($val =~ m/ID=(.+)$/) {
                $id = $1;
                last;
            }
        }
    }
    $id .= ".".$feat->id;
    return $id;
}
