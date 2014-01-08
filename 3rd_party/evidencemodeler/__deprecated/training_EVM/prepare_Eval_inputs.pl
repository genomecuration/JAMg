#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../PerlLib");
use CdbTools;
use Gene_obj_indexer;
use GFF3_utils;
use Getopt::Long qw(:config no_ignore_case bundling);
use Carp;
use Fasta_reader;

umask(0000);

my $usage = <<_EOUSAGE_;

####################################################################################################
#  
#  --dir_listing          file containing list of EVM training directories (ie. all_train.entries.evaluate)
#
#  --template_ev_type     evidence type for the reference gene type.
#  --other_ev_type        evidence type (ie. EVM or Genscan)
#
#  --template_filename    name of reference gene structure file 
#  --other_filename        name of gff3 file in the EVM training directory.
#  --genome_filename      name of the genome fasta sequence file
#
#  --intergenic_flank     number of bases to include flanking the template gene structure.
#####################################################################################################


_EOUSAGE_

    ;

my ($dir_listing, $template_ev_type, $other_ev_type, 
    $template_filename, $other_filename, $intergenic_flank,
    $genome_filename);

&GetOptions ( 'dir_listing=s' => \$dir_listing,
              'template_ev_type=s' => \$template_ev_type,
              'other_ev_type=s' => \$other_ev_type,
              'template_filename=s' => \$template_filename,
              'other_filename=s' => \$other_filename,
              'intergenic_flank=i' => \$intergenic_flank,
              'genome_filename=s' => \$genome_filename,
              );


unless ($dir_listing 
        && $template_ev_type && $template_filename 
        && $other_ev_type && $other_filename
        && $genome_filename 
        && $intergenic_flank) { 
    die $usage;
}


main: {

    ## open output files for writing.
    my $reference_GTF_filename = "reference.GTF";
    open (my $reference_fh, ">$reference_GTF_filename") or die "Error, cannot write to $reference_GTF_filename";
    
    my $reference_fasta_filename = "reference.GTF.fasta";
    open (my $reference_genome_fh, ">$reference_fasta_filename") or die "Error, cannot write to $reference_fasta_filename";
    
    my $other_GTF_filename = "${other_ev_type}.GTF";
    open (my $other_fh, ">$other_GTF_filename") or die "Error, cannot write to $other_GTF_filename";
    
    

    open (my $fh, $dir_listing) or die "Error, cannot open $dir_listing";

    my $postfix_counter = 0;

    while (<$fh>) {
        my $dirname = $_;
        $dirname =~ s/\s+//g;
        
        ## get the template gene object.
        my $template_gene_obj_indexer = new Gene_obj_indexer( {create => "template.$$.inx"} );
        my $template_gff3_file = "$dirname/$template_filename";
        &GFF3_utils::index_GFF3_gene_objs($template_gff3_file, $template_gene_obj_indexer);
        
        ## should only be one gene.
        my @genes = $template_gene_obj_indexer->get_keys();
        if (scalar @genes != 1) {
            confess "Error, got more than one gene from $template_gff3_file";
        }
        
        my $template_gene_id = $genes[0];
        my $template_gene_obj = $template_gene_obj_indexer->get_gene($template_gene_id);
        my ($gene_lend, $gene_rend) = sort {$a<=>$b} $template_gene_obj->get_model_span();
        

        ## get the genome sequence
        my $genome_file = "$dirname/$genome_filename";
        my $fasta_reader = new Fasta_reader($genome_file);
        my $seqobj = $fasta_reader->next();
        my $genome_sequence = $seqobj->get_sequence();

        
        ## get the other gene object.
        my $other_gene_obj_indexer = new Gene_obj_indexer( { create => "other.$$.inx"} );
        my $other_filename = "$dirname/$other_filename";
        &GFF3_utils::index_GFF3_gene_objs($other_filename, $other_gene_obj_indexer);
        
        my @other_gene_ids = $other_gene_obj_indexer->get_keys();
        ## want just those of the proper source type
        my @other_gene_objs;
        foreach my $other_gene_id (@other_gene_ids) {
            my $gene_obj = $other_gene_obj_indexer->get_gene($other_gene_id);
            # print $gene_obj->toString();
            
            my $source = $gene_obj->{source};
            if ($source eq $other_ev_type) {
                push (@other_gene_objs, $gene_obj);
            }
        }
        
        ## write GTF to temp files
        my $tmp_reference_out = "tmp.$$.reference";
        open (my $tmp_ref_fh, ">$tmp_reference_out") or die "Error, cannot write to $tmp_reference_out";
        my $tmp_other_out = "tmp.$$.$other_ev_type";
        open (my $tmp_other_fh, ">$tmp_other_out") or die "Error, cannot write to $tmp_other_out";
        
        
        $postfix_counter++; ## add postfix to the contig name and the gene/model identifiers.
        #                     this is done primarily since we're later shifting everything down
        #                     to a restricted set of coordinates and models from different regions
        #                     of the same contig will collide in the new sequence space.
        
        $template_gene_obj->{asmbl_id} .= ".$postfix_counter";
        my $ref_contig_id = $template_gene_obj->{asmbl_id};
        $template_gene_obj->{TU_feat_name} .= ".$postfix_counter";
        $template_gene_obj->{Model_feat_name} .= ".$postfix_counter";
        
        &write_GTF($template_gene_obj, \$genome_sequence, "reference", $tmp_ref_fh);
                
        foreach my $other_gene_obj (@other_gene_objs) {
            $other_gene_obj->{asmbl_id} .= ".$postfix_counter";
            $other_gene_obj->{TU_feat_name} .= ".$postfix_counter";
            $other_gene_obj->{Model_feat_name} .= ".$postfix_counter";
            
            &write_GTF($other_gene_obj, \$genome_sequence, $other_ev_type, $tmp_other_fh);
        }
        
        close $tmp_ref_fh;
        close $tmp_other_fh;
        
        ## extract just the ranges within the genic region and flank
        my $region_lend = $gene_lend - $intergenic_flank;
        $region_lend = 1 if ($region_lend < 1);
        my $region_rend = $gene_rend + $intergenic_flank;

        &extract_region_from_GTF([$region_lend, $region_rend], $tmp_reference_out, $reference_fh);
        &extract_region_from_GTF([$region_lend, $region_rend], $tmp_other_out, $other_fh);
        
        my $subseq = substr($genome_sequence, $region_lend - 1, $region_rend - $region_lend + 1);
        $subseq =~ s/(\S{60})/$1\n/g;
        chomp $subseq;
        
        print $reference_genome_fh ">$ref_contig_id\n$subseq\n";
        
        unlink ($tmp_reference_out, $tmp_other_out, "other.$$.inx", "template.$$.inx");
        

    }
}

exit(0);

####
sub extract_region_from_GTF {
    my ($range_aref, $filename, $filehandle) = @_;

    my ($range_lend, $range_rend) = @$range_aref;

    open (my $fh, $filename) or die "Error, cannot open file $filename";
    while (<$fh>) {
        
        unless (/\w/) { next; }
        my @x = split (/\t/);
        ## make sure it overlaps the range
        my ($lend, $rend) = ($x[3], $x[4]);
        unless ($lend <= $range_rend && $rend >= $range_lend) {
            # no overlap
            next; 
        }

        ## adjust coordinates to range
        if ($lend < $range_lend) {
            $lend = $range_lend;
        }
        if ($rend > $range_rend) {
            $rend = $range_rend;
        }

        ## adjust to start at position 1
        $lend -= ($range_lend -1);
        $rend -= ($range_lend -1);

        $x[3] = $lend;
        $x[4] = $rend;

        print $filehandle join ("\t", @x);
    }
    close $fh;
    
    return;
}


####
sub write_GTF {
    my ($gene_obj, $genome_sequence_ref, $ev_type, $filehandle) = @_;
        
    eval {
        print $filehandle $gene_obj->to_GTF_format($genome_sequence_ref, source => $ev_type);
    };
    
    if ($@) {
        # something's not exactly perfect about this gene structure.  Report GTF in pseudogene format
        my $feat_name = $gene_obj->{Model_feat_name};
        print STDERR "Warning: prediction $feat_name, $ev_type is imperfect and encodes an interrupted ORF\n";
        $gene_obj->{is_pseudogene} = 1;
        print $filehandle $gene_obj->to_GTF_format($genome_sequence_ref, source => $ev_type);
    }

    return;
}
