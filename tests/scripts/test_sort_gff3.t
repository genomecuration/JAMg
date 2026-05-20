#!/usr/bin/env perl
# Unit test for bin/sort_gff3.pl
# Script: sort_gff3.pl <GFF3> [delimiter]
# Behaviour: writes <GFF3>.sorted with gene records sorted by seqid then by
# gene start. Features within a gene record reordered by sort_order
# (gene<mRNA<exon/intron/UTR<CDS<nucleotide_to_protein_match).

use strict;
use warnings;
use Test::More tests => 4;
use File::Temp qw(tempdir);
use FindBin qw($Bin);

my $script = "$Bin/../../bin/sort_gff3.pl";
my $tmp = tempdir(CLEANUP => 1);

# ---------------------------------------------------------------------------
# Case 1: two genes on different seqids in reverse order -> sorted by seqid.
# Each gene block separated by blank lines (record delimiter).
# ---------------------------------------------------------------------------
{
    my $infile = "$tmp/case1.gff3";
    open my $fh, '>', $infile or die $!;
    print $fh "##gff-version 3\n\n";  # header line for delimiter discovery
    # second gene first (on ctgB)
    print $fh "ctgB\tmaker\tgene\t500\t900\t.\t+\t.\tID=gB\n";
    print $fh "ctgB\tmaker\tmRNA\t500\t900\t.\t+\t.\tID=mB;Parent=gB\n";
    print $fh "ctgB\tmaker\texon\t500\t900\t.\t+\t.\tID=eB.1;Parent=mB\n\n";
    # then first gene (on ctgA)
    print $fh "ctgA\tmaker\tgene\t100\t300\t.\t+\t.\tID=gA\n";
    print $fh "ctgA\tmaker\tmRNA\t100\t300\t.\t+\t.\tID=mA;Parent=gA\n";
    print $fh "ctgA\tmaker\texon\t100\t300\t.\t+\t.\tID=eA.1;Parent=mA\n\n";
    close $fh;

    system("perl $script $infile") == 0 or die "case1 script failed: $?";
    my $sorted = "$infile.sorted";
    ok(-s $sorted, "case 1: sorted file produced");
    my $content = do { open my $rh, '<', $sorted or die $!; local $/; <$rh> };
    my @gene_lines = grep { /\tgene\t/ } split /\n/, $content;
    is(scalar(@gene_lines), 2, "case 1: two gene lines present");
    like(
        $content,
        qr/ctgA.+ID=gA.+ctgB.+ID=gB/s,
        "case 1: ctgA gene appears before ctgB gene after sort"
    );
}

# ---------------------------------------------------------------------------
# Case 2: feature reordering within a record - CDS should come after exon.
# We feed a record with CDS before exon and check the order flips.
# ---------------------------------------------------------------------------
{
    my $infile = "$tmp/case2.gff3";
    open my $fh, '>', $infile or die $!;
    print $fh "##gff-version 3\n\n";
    print $fh "ctgA\tmaker\tgene\t100\t300\t.\t+\t.\tID=gA\n";
    print $fh "ctgA\tmaker\tmRNA\t100\t300\t.\t+\t.\tID=mA;Parent=gA\n";
    print $fh "ctgA\tmaker\tCDS\t100\t300\t.\t+\t0\tID=cA.1;Parent=mA\n";
    print $fh "ctgA\tmaker\texon\t100\t300\t.\t+\t.\tID=eA.1;Parent=mA\n\n";
    close $fh;

    system("perl $script $infile") == 0 or die "case2 script failed";
    my $sorted = "$infile.sorted";
    my $content = do { open my $rh, '<', $sorted or die $!; local $/; <$rh> };
    # sort_order puts exon=3, CDS=5, so exon must precede CDS within record.
    like(
        $content,
        qr/\texon\t.+\tCDS\t/s,
        "case 2: exon row precedes CDS row after feature-level sort"
    );
}
