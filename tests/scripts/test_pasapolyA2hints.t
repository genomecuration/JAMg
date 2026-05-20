#!/usr/bin/env perl
# Unit test for bin/pasapolyA2hints.pl.
# Signature: pasapolyA2hints.pl <polyAsites.fasta> [--radius=N]
# Behaviour: each polyA-FASTA header matching the regex
#   ^(.*)-([\d]*)_([+\-])\s([\d]*)\s+transcripts:\s([^,]*)
# is converted to a GFF line with feature 'tts' and the 'src=E;EST=...;pri=4'
# attributes. Writes to STDOUT.
#
# (Plan note: plan-text said "exonpart"; actual feature is 'tts'. We test the
# script's real output.)

use strict;
use warnings;
use Test::More tests => 4;
use File::Temp qw(tempdir);
use FindBin qw($Bin);

my $script = "$Bin/../../bin/pasapolyA2hints.pl";
my $tmp = tempdir(CLEANUP => 1);

# ---------------------------------------------------------------------------
# Build a tiny polyA-fasta with two records.
# ---------------------------------------------------------------------------
my $fa = "$tmp/polyA.fasta";
{
    open my $fh, '>', $fa or die $!;
    print $fh ">ctgA-1234_+ 8 transcripts: asmbl_1,asmbl_2,asmbl_3\nNNNNNNNNNN\n";
    print $fh ">ctgB-5678_- 3 transcripts: asmbl_5\nNNNNNNNN\n";
    close $fh;
}

# ---------------------------------------------------------------------------
# Case 1: default radius -> two 'tts' GFF lines on STDOUT with src=E.
# ---------------------------------------------------------------------------
{
    my $out_f = "$tmp/case1.out";
    my $rc = system("perl $script $fa > $out_f 2>/dev/null");
    is($rc, 0, "case 1: exits 0");
    my $content = do { open my $r, '<', $out_f or die $!; local $/; <$r> };
    my @tts_lines = grep { /\ttts\t/ } split /\n/, $content;
    ok(scalar(@tts_lines) >= 2, "case 1: produced >= 2 tts lines (got " . scalar(@tts_lines) . ")");
    like($content, qr/src=E/, "case 1: hints carry src=E attribute");
}

# ---------------------------------------------------------------------------
# Case 2: --radius=10 widens the start/end window. With ctgA-1234 plus radius
# 10 the GFF row should report start=1224 end=1244.
# ---------------------------------------------------------------------------
{
    my $out_f = "$tmp/case2.out";
    my $rc = system("perl $script --radius=10 $fa > $out_f 2>/dev/null");
    my $content = do { open my $r, '<', $out_f or die $!; local $/; <$r> };
    my ($row) = grep { /^ctgA\t/ } split /\n/, $content;
    if ($row) {
        my @cols = split /\t/, $row;
        is($cols[3] . '..' . $cols[4], '1224..1244',
            "case 2: --radius=10 produces start=1224, end=1244 for pos 1234");
    } else {
        fail("case 2: no ctgA row found");
    }
}
