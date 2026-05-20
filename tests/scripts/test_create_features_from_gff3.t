#!/usr/bin/env perl
# Unit test for bin/create_features_from_gff3.pl.
# Signature: create_features_from_gff3.pl -gff <GFF3> -genome <FASTA> [-rename] [-strip_name] [-fix_first_phase]
# Behaviour: indexes the GFF, builds Gene_obj's, emits <GFF3>.gff3 + .pep.fasta
# + .cds.fasta + .mRNA.fasta etc.

use strict;
use warnings;
use Test::More tests => 11;
use File::Temp qw(tempdir);
use FindBin qw($Bin);

my $script = "$Bin/../../bin/create_features_from_gff3.pl";
my $tmp = tempdir(CLEANUP => 1);

# ---------------------------------------------------------------------------
# Build a 600-nt genome with a single-CDS gene at positions 101..199 on +
# strand. ATG at 101..103, in-frame TAA at 197..199. Sequence outside the CDS
# is filler.
# ---------------------------------------------------------------------------
my $genome = "$tmp/genome.fa";
my $gff    = "$tmp/case.gff3";
{
    # Build a 96-codon ORF: ATG + 30 codons of arbitrary nt + TAA  (gives M..*).
    # Length = 3 + 90 + 3 = 96 nt at positions 101..196.
    my $cds_seq = 'ATG' . ('GCT' x 30) . 'TAA';  # 96 nt
    is(length($cds_seq), 96, "internal: CDS sequence is 96 nt");
    my $genome_seq = ('A' x 100) . $cds_seq . ('A' x 100);
    open my $fh, '>', $genome or die $!;
    print $fh ">ctgA\n$genome_seq\n";
    close $fh;

    open my $g, '>', $gff or die $!;
    print $g "##gff-version 3\n";
    print $g "ctgA\tmaker\tgene\t101\t196\t.\t+\t.\tID=g1;Name=gene_name_1\n";
    print $g "ctgA\tmaker\tmRNA\t101\t196\t.\t+\t.\tID=m1;Parent=g1;Name=gene_name_1\n";
    print $g "ctgA\tmaker\texon\t101\t196\t.\t+\t.\tID=e1.1;Parent=m1\n";
    print $g "ctgA\tmaker\tCDS\t101\t196\t.\t+\t0\tID=cds1.1;Parent=m1\n";
    close $g;
}

# ---------------------------------------------------------------------------
# Case 1: -rename rewrites gene IDs to JAMg_model_N.
# ---------------------------------------------------------------------------
{
    my $gff1 = "$tmp/case1.gff3";
    system("cp", $gff, $gff1) == 0 or die $!;
    my $rc = system("cd $tmp && perl $script -gff $gff1 -genome $genome -rename >case1.log 2>&1");
    is($rc, 0, "case 1 (-rename): exits 0");
    ok(-s "$gff1.gff3", "case 1: <gff>.gff3 output written");
    my $content = do { open my $r, '<', "$gff1.gff3" or die $!; local $/; <$r> };
    like($content, qr/ID=JAMg_model_/, "case 1: -rename rewrites gene IDs as JAMg_model_*");
}

# ---------------------------------------------------------------------------
# Case 2: -strip_name keeps the script happy on a Name-bearing GFF (vanilla
# run produces gff3 + pep with M..*).
# ---------------------------------------------------------------------------
{
    my $gff2 = "$tmp/case2.gff3";
    system("cp", $gff, $gff2) == 0 or die $!;
    my $rc = system("cd $tmp && perl $script -gff $gff2 -genome $genome -strip_name >case2.log 2>&1");
    is($rc, 0, "case 2 (-strip_name): exits 0");
    ok(-s "$gff2.gff3", "case 2: <gff>.gff3 output written");
    my $pep = "$gff2.pep.fasta";
    ok(-s $pep, "case 2: <gff>.pep.fasta written");
    my $pep_content = do { open my $r, '<', $pep or die $!; local $/; <$r> };
    like($pep_content, qr/^M.*\*$/m, "case 2: peptide starts with M and ends with *");
}

# ---------------------------------------------------------------------------
# Case 3: -fix_first_phase exercises the phase-correction branch; verify the
# resulting GFF still emits a coherent gene record.
# ---------------------------------------------------------------------------
{
    my $gff3f = "$tmp/case3.gff3";
    system("cp", $gff, $gff3f) == 0 or die $!;
    my $rc = system("cd $tmp && perl $script -gff $gff3f -genome $genome -fix_first_phase >case3.log 2>&1");
    is($rc, 0, "case 3 (-fix_first_phase): exits 0");
    ok(-s "$gff3f.gff3", "case 3: <gff>.gff3 output written");
    my $content = do { open my $r, '<', "$gff3f.gff3" or die $!; local $/; <$r> };
    like($content, qr/\tCDS\t/, "case 3: -fix_first_phase still emits a CDS record");
}
