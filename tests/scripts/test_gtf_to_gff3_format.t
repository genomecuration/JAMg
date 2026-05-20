#!/usr/bin/env perl
# Unit test for bin/gtf_to_gff3_format.pl.
# Signature: gtf_to_gff3_format.pl <GTF> [genome.fasta] [source] > <GFF3>
# Reads GTF, indexes gene objects via GTF_utils, emits GFF3 to STDOUT.

use strict;
use warnings;
use Test::More tests => 6;
use File::Temp qw(tempdir);
use FindBin qw($Bin);

my $script = "$Bin/../../bin/gtf_to_gff3_format.pl";
my $tmp = tempdir(CLEANUP => 1);

# ---------------------------------------------------------------------------
# Build a minimal GeneMark/Augustus-style GTF with one gene/one transcript and
# two CDS features (mid-intron split). This is enough for the GTF_utils indexer
# to produce a gene_obj with two CDS exons.
# ---------------------------------------------------------------------------
my $gtf = "$tmp/case.gtf";
{
    open my $fh, '>', $gtf or die $!;
    print $fh qq(chr1\tAUGUSTUS\tCDS\t1100\t1500\t0.9\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n);
    print $fh qq(chr1\tAUGUSTUS\tCDS\t1700\t1900\t0.9\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n);
    print $fh qq(chr1\tAUGUSTUS\tstop_codon\t1898\t1900\t.\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n);
    close $fh;
}

# ---------------------------------------------------------------------------
# Case 1: GTF -> GFF3 emits gene, mRNA, CDS records to STDOUT.
# ---------------------------------------------------------------------------
{
    my $out_f = "$tmp/case1.gff3";
    my $rc = system("perl $script $gtf > $out_f 2>/dev/null");
    is($rc, 0, "case 1: gtf_to_gff3 exits 0");
    ok(-s $out_f, "case 1: emitted GFF3 is non-empty");
    my $content = do { open my $r, '<', $out_f or die $!; local $/; <$r> };
    like($content, qr/\tgene\t/, "case 1: emitted at least one gene feature");
    like($content, qr/\tmRNA\t|\tCDS\t/, "case 1: emitted at least one mRNA or CDS feature");
}

# ---------------------------------------------------------------------------
# Case 2: passing a user-supplied source as the 3rd positional arg rewrites
# the GFF3 source column. We pass dummy.fa as genome to satisfy positional
# order; we use a 600bp scaffold so GTF coords (1100..1900) still resolve.
# ---------------------------------------------------------------------------
{
    my $genome = "$tmp/genome.fa";
    open my $fh, '>', $genome or die $!;
    print $fh ">chr1\n", ('A' x 2500), "\n";
    close $fh;
    my $out_f = "$tmp/case2.gff3";
    my $rc = system("perl $script $gtf $genome MYSOURCE > $out_f 2>/dev/null");
    is($rc, 0, "case 2: exits 0 with user-source");
    my $content = do { open my $r, '<', $out_f or die $!; local $/; <$r> };
    like($content, qr/\tMYSOURCE\t/, "case 2: user-supplied source token appears in column 2");
}
