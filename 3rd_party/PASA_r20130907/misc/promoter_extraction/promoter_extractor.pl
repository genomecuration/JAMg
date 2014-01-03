#!/usr/bin/perl

use strict;
use lib ($ENV{EUK_MODULES}, $ENV{EGC_SCRIPTS});
use CdbTools;
use Egc_library;

my $usage = "usage: $0 alignments_described_file  genome_fasta_file\n\n";

my $alignments_described_file = $ARGV[0] or die $usage;
my $genome_fasta_file = $ARGV[1] or die $usage;


my $PROMOTER_LENGTH = 2000;

my %asmbl_id_to_assembly_struct;


open (my $fh, $alignments_described_file) or die "Error, cannot open $alignments_described_file";
while (<$fh>) {
    chomp;
    my ($pasa_acc, $cdna_list, $genomic_acc, $alignment_txt) = split (/\t/);
    
    $pasa_acc =~ s/alignment_assembly://;
    $genomic_acc =~ s/genomic_acc://;
    
    my $alignment_struct = &parse_alignment_struct($alignment_txt);

    my $struct_list_aref = $asmbl_id_to_assembly_struct{$genomic_acc};
    unless (ref $struct_list_aref) {
        $struct_list_aref = $asmbl_id_to_assembly_struct{$genomic_acc} = [];
    }

    ## pack more info into the struct:
    $alignment_struct->{alignment_text} = $alignment_txt;
    $alignment_struct->{genomic_acc} = $genomic_acc;
    $alignment_struct->{pasa_acc} = $pasa_acc;
    $alignment_struct->{cdna_list} = $cdna_list;
    

    push (@$struct_list_aref, $alignment_struct);

}

## extract promoter regions:
foreach my $genomic_acc (keys %asmbl_id_to_assembly_struct) {
    my $struct_list_aref = $asmbl_id_to_assembly_struct{$genomic_acc};

    my $genomic_seq = &cdbyank_linear($genomic_acc, $genome_fasta_file);

    &extract_promoters(\$genomic_seq, $struct_list_aref);

}

exit(0);


####
sub extract_promoters {
    my ($genomic_seq_ref, $alignment_structs_aref) = @_;
    
    foreach my $alignment_struct (@$alignment_structs_aref) {
        
        ## unpack the struct:
        my ($lend, $rend, $orient) = ($alignment_struct->{lend}, $alignment_struct->{rend}, $alignment_struct->{orient});
        my $alignment_text = $alignment_struct->{alignment_text};
        my $genomic_acc = $alignment_struct->{genomic_acc};
        my $pasa_acc = $alignment_struct->{pasa_acc};
        my $cdna_list = $alignment_struct->{cdna_list};


        my $promoter_seq = "";
        if ($orient eq '+') {
            ## want upstream sequence:
            
            my $start_pos = $lend - $PROMOTER_LENGTH;
            if ($start_pos < 1) {
                $start_pos = 1;
            }

            my $extract_length = $lend - $start_pos;
            
            $promoter_seq = substr($$genomic_seq_ref, $start_pos -1, $extract_length);
        }
        else {
            ## minus strand
            # want downstream seq and revcomp it.
            my $start_pos = $rend + 1;
            my $extract_length = $PROMOTER_LENGTH;

            $promoter_seq = substr($$genomic_seq_ref, $start_pos - 1, $extract_length);
            $promoter_seq = &reverse_complement($promoter_seq);

        }
        
        $promoter_seq = &make_FASTA_format($promoter_seq);
        chomp $promoter_seq;
        
        print ">${pasa_acc}_prom$PROMOTER_LENGTH $genomic_acc $alignment_text $cdna_list\n$promoter_seq\n";
    }
    
    return;
    
}



####
sub parse_alignment_struct {
    my ($alignment_txt) = @_;

    my $orientation;

    $alignment_txt =~ m|orient\(a([\+\-])/s([\+\-\?])\)| or die "Error, couldn't match orientation info.";
    my ($aligned_orient, $spliced_orient) = ($1, $2);

    if ($spliced_orient ne '?') {
        $orientation = $spliced_orient;
    }
    else {
        $orientation = $aligned_orient;
    }

    my @coords;
    while ($alignment_txt =~ m|(\d+)\(\d+\)-(\d+)\(\d+\)|g ) {
        push (@coords, $1, $2);
    }
    @coords = sort {$a<=>$b} @coords;
    my $lend = shift @coords;
    my $rend = pop @coords;

    return ( { orient => $orientation,
               lend => $lend,
               rend => $rend,
           }
             );
}




    






