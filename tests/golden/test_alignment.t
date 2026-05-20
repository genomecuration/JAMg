#!/usr/bin/env perl
# Unit test for Golden::Alignment::correct_exonerate_gff (the parser leg of
# Alignment.pm — the actual `run_exonerate` orchestrator launches exonerate
# which is the external boundary, see spec §4.10).
#
# Cases (3):
#   1. correct_exonerate_gff on a protein2genome fixture produces a non-empty
#      <input>.corrected.gff3.
#   2. correct_exonerate_gff dies when the fixture header says cdna2genome but
#      the driver flag $main::is_cdna is false (mismatched flag).
#   3. correct_exonerate_gff dies when fixture header says protein2genome but
#      $main::is_cdna is true.

use strict;
use warnings;
use Test::More tests => 5;
use File::Temp qw(tempdir);
use FindBin qw($Bin);

# Load the module from the repo tree.
use lib "$Bin/../../PerlLib";
use Golden::Alignment qw(correct_exonerate_gff);

my $fixture = "$Bin/../../test_suite/scripts/mini-exonerate-output.txt";
ok(-s $fixture, "fixture mini-exonerate-output.txt present");

my $tmp = tempdir(CLEANUP => 1);

# ---------------------------------------------------------------------------
# Driver globals needed by correct_exonerate_gff. Set the minimum required.
# ---------------------------------------------------------------------------
$main::identical_fraction_cutoff = 95;
$main::similar_fraction_cutoff   = 95;
$main::mismatch_cutoff           = 10;
$main::is_cdna                   = 0;     # protein2genome fixture
$main::liberal_cutoffs           = 0;
$main::only_complete             = 0;
$main::no_single_exon            = 0;
$main::nodataprint               = 1;
$main::failed_cutoff             = 0;
$main::scaffold_seq_hashref      = { ctgA => 'N' x 1000 };
$main::genome_sequence_file      = "$tmp/dummy.fa";
$main::stop_after_correction     = 0;
{
    open my $fh, '>', $main::genome_sequence_file or die $!;
    print $fh ">ctgA\n", ('N' x 1000), "\n";
    close $fh;
}

# ---------------------------------------------------------------------------
# Case 1: protein2genome fixture + is_cdna=0 -> <input>.corrected.gff3 written.
# Pre-create the .corrected.passed file so the trailing create_golden_gffs
# call short-circuits and does not try to read .gene/.pep helpers.
# ---------------------------------------------------------------------------
{
    my $in = "$tmp/case1-exon.txt";
    system("cp", $fixture, $in) == 0 or die "cp: $?";
    # short-circuit the golden-set creation: just give it any non-empty file.
    open my $passed, '>', "$in.corrected.passed" or die $!;
    print $passed "stub\n";
    close $passed;

    correct_exonerate_gff($in);
    my $out = "$in.corrected.gff3";
    ok(-s $out, "case 1: corrected.gff3 produced");
    my $content = do { open my $r, '<', $out or die $!; local $/; <$r> };
    like($content, qr/\tgene\t/, "case 1: corrected.gff3 contains a gene line");
}

# ---------------------------------------------------------------------------
# Case 2: protein2genome fixture but $main::is_cdna=1 -> dies with explicit error.
# ---------------------------------------------------------------------------
{
    local $main::is_cdna = 1;
    my $in = "$tmp/case2-exon.txt";
    system("cp", $fixture, $in) == 0 or die "cp: $?";
    open my $passed, '>', "$in.corrected.passed" or die $!;
    print $passed "stub\n";
    close $passed;

    eval { correct_exonerate_gff($in); };
    like($@, qr/protein2genome.+but -cdna was requested/i,
        "case 2: protein2genome + is_cdna=1 dies with explanatory error");
}

# ---------------------------------------------------------------------------
# Case 3: fixture with cdna2genome header + $main::is_cdna=0 -> dies.
# ---------------------------------------------------------------------------
{
    my $in = "$tmp/case3-exon.txt";
    open my $w, '>', $in or die $!;
    open my $r, '<', $fixture or die $!;
    while (my $ln = <$r>) {
        $ln =~ s/protein2genome/cdna2genome/g;
        print $w $ln;
    }
    close $r;
    close $w;
    open my $passed, '>', "$in.corrected.passed" or die $!;
    print $passed "stub\n";
    close $passed;

    local $main::is_cdna = 0;
    eval { correct_exonerate_gff($in); };
    like($@, qr/cdna2genome.+but -cdna was not requested/i,
        "case 3: cdna2genome + is_cdna=0 dies with explanatory error");
}
