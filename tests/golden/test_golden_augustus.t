#!/usr/bin/env perl
# Unit test for Golden::Augustus.
#
# Per plan §4.10: Golden::Augustus has no write_gb sub (verified by grep);
# the spec's fall-back is "test parse_gb + gb2gff3 round-trip". This test
# (3 cases) exercises that round-trip on a synthetic 2-locus GenBank input
# with multi-exon CDS on both strands.

use strict;
use warnings;
use Test::More tests => 7;
use File::Temp qw(tempdir);
use FindBin qw($Bin);

use lib "$Bin/../../PerlLib";
use Bio::SeqIO;
use Bio::Seq;
use Bio::SeqFeature::Generic;
use Bio::Location::Split;
use Bio::Location::Simple;
use Golden::Filter;   # provides Golden::Filter::process_cmd + sort_gff3 used internally by Augustus.pm
use Golden::Augustus qw(parse_gb gb2gff3);

# parse_gb shells out via $main::cdbfasta_exec — locate the binary so the
# round-trip can complete. cdbfasta lives in the SIF at /opt/jamg/bin.
chomp(my $cdbfasta = `which cdbfasta 2>/dev/null`);
SKIP: {
    skip "cdbfasta not on PATH (run under SIF)", 7 unless $cdbfasta && -x $cdbfasta;

    $main::cdbfasta_exec = $cdbfasta;
    my $tmp = tempdir(CLEANUP => 1);

    # ----------------------------------------------------------------------
    # Build a 2-locus GenBank: locus_plus has a multi-exon CDS on the +
    # strand; locus_minus has a multi-exon CDS on the - strand.
    # ----------------------------------------------------------------------
    my $gb_file = "$tmp/test.gb";
    {
        my $out = Bio::SeqIO->new(-file => ">$gb_file", -format => "genbank");

        # Locus 1: plus strand
        my $seq1 = Bio::Seq->new(-id => "locus_plus", -seq => "ATGC" x 200, -alphabet => "dna");
        my $loc1 = Bio::Location::Split->new();
        $loc1->add_sub_Location(Bio::Location::Simple->new(-start => 10, -end => 100, -strand => 1));
        $loc1->add_sub_Location(Bio::Location::Simple->new(-start => 200, -end => 300, -strand => 1));
        my $f1 = Bio::SeqFeature::Generic->new(-primary_tag => "CDS", -location => $loc1);
        $seq1->add_SeqFeature($f1);
        $out->write_seq($seq1);

        # Locus 2: minus strand
        my $seq2 = Bio::Seq->new(-id => "locus_minus", -seq => "GATC" x 200, -alphabet => "dna");
        my $loc2 = Bio::Location::Split->new();
        $loc2->add_sub_Location(Bio::Location::Simple->new(-start => 50, -end => 150, -strand => -1));
        $loc2->add_sub_Location(Bio::Location::Simple->new(-start => 250, -end => 400, -strand => -1));
        my $f2 = Bio::SeqFeature::Generic->new(-primary_tag => "CDS", -location => $loc2);
        $seq2->add_SeqFeature($f2);
        $out->write_seq($seq2);
    }
    ok(-s $gb_file, "GenBank fixture written");

    # ----------------------------------------------------------------------
    # Case 1: parse_gb returns a hashref keyed by locus id with exon arrays.
    # ----------------------------------------------------------------------
    my $href = parse_gb($gb_file);
    ok(exists $href->{locus_plus},  "case 1: parse_gb has locus_plus key");
    ok(exists $href->{locus_minus}, "case 1: parse_gb has locus_minus key");
    is($href->{locus_plus}{strand},  1,  "case 1: plus locus strand=1");
    is($href->{locus_minus}{strand}, -1, "case 1: minus locus strand=-1");

    # ----------------------------------------------------------------------
    # Case 2: gb2gff3 round-trip writes a valid sorted GFF3 with exon rows.
    # ----------------------------------------------------------------------
    my $gff_out = "$tmp/out.gff3";
    gb2gff3($href, $gff_out);
    ok(-s $gff_out, "case 2: gb2gff3 produces non-empty GFF3");

    # ----------------------------------------------------------------------
    # Case 3: GFF rows reference both loci and emit at least one exon row each.
    # The function emits per-exon lines plus a span row using the script's
    # canonical "\tGB\texon\t" format.
    # ----------------------------------------------------------------------
    my $content = do { open my $r, '<', $gff_out or die $!; local $/; <$r> };
    my $plus_rows  = () = ($content =~ /^locus_plus\tGB\texon\t/gm);
    my $minus_rows = () = ($content =~ /^locus_minus\tGB\texon\t/gm);
    ok(
        $plus_rows >= 2 && $minus_rows >= 2,
        "case 3: both loci emit >=2 exon rows each (plus=$plus_rows, minus=$minus_rows)"
    );
}
