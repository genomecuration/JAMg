#!/usr/bin/env perl
# Unit test for bin/augustus_RNAseq_hints.pl.
# Signature: augustus_RNAseq_hints.pl -bam <BAM> -genome <FASTA> -cpus N [other flags]
# Side effects: runs samtools, bedtools; produces <BAM>.coverage.hints
# in Augustus exonpart hints format.
#
# (Plan note: plan-text said "src=E features"; the script actually emits
# `src=RCOV` for coverage hints — we test the script's real behavior.)
#
# This test must run in the SIF (samtools + bedtools) and uses $TMP from
# the environment for the script's TMP requirement.

use strict;
use warnings;
use Test::More tests => 4;
use File::Temp qw(tempdir);
use FindBin qw($Bin);

my $script = "$Bin/../../bin/augustus_RNAseq_hints.pl";
my $repo   = "$Bin/../..";
my $bam    = "$repo/test_suite/mini-rnaseq.bam";
my $genome = "$repo/test_suite/mini-genome.fasta";

SKIP: {
    skip "mini fixtures missing", 4 unless -s $bam && -s $genome;

    my $tmp = tempdir(CLEANUP => 1);
    # script needs TMP/TMPDIR exported
    local $ENV{TMP} = $tmp;
    # copy inputs to a writable dir; the script writes side-files
    # alongside the bam and the fasta.
    my $bam_copy    = "$tmp/mini-rnaseq.bam";
    my $genome_copy = "$tmp/mini-genome.fasta";
    system("cp", $bam, $bam_copy)       == 0 or die "cp bam: $?";
    system("cp", $genome, $genome_copy) == 0 or die "cp genome: $?";
    # bai might be needed for samtools
    system("cp", "$bam.bai", "$bam_copy.bai") if -e "$bam.bai";

    # ------------------------------------------------------------------
    # Case 1: -no_junctions (skips junction subroutine) keeps the run lean.
    # The script produces a <BAM>.coverage.hints file in Augustus hints format.
    # ------------------------------------------------------------------
    my $rc = system("perl $script -bam $bam_copy -genome $genome_copy "
                  . "-cpus 1 -no_junctions >$tmp/run.log 2>&1");
    is($rc, 0, "case 1: exits 0 (see $tmp/run.log on failure)");
    my $hints = "$bam_copy.coverage.hints";
    ok(-s $hints, "case 1: coverage hints file produced");

    # ------------------------------------------------------------------
    # Case 2: hints file is Augustus-compatible: 9 tab-separated columns and
    # the script's expected `src=RCOV` token.
    # ------------------------------------------------------------------
    my $content = do { open my $r, '<', $hints or die $!; local $/; <$r> };
    my @lines   = grep { !/^\s*$/ } split /\n/, $content;
    cmp_ok(scalar(@lines), '>', 0, "case 2: hints file contains feature lines");
    like($content, qr/\tsrc=RCOV;pri=4/,
        "case 2: hints carry the expected src=RCOV;pri=4 attribute");
}
