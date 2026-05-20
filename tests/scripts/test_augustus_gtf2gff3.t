#!/usr/bin/env perl
# Unit test for bin/augustus_gtf2gff3.pl.
# Signature: augustus_gtf2gff3.pl --out <FILE> [--gff3] < input.gtf
# Reads GTF on STDIN, writes either GTF or GFF3 to --out depending on --gff3 flag.

use strict;
use warnings;
use Test::More tests => 8;
use File::Temp qw(tempdir);
use FindBin qw($Bin);

my $script = "$Bin/../../bin/augustus_gtf2gff3.pl";
my $tmp = tempdir(CLEANUP => 1);

# ---------------------------------------------------------------------------
# Build an Augustus-style GTF with gene/transcript/CDS/intron lines.
# ---------------------------------------------------------------------------
my $gtf = "$tmp/aug.gtf";
{
    open my $fh, '>', $gtf or die $!;
    print $fh qq(chr1\tAUGUSTUS\tgene\t1000\t2000\t0.5\t+\t.\tg1\n);
    print $fh qq(chr1\tAUGUSTUS\ttranscript\t1000\t2000\t0.5\t+\t.\tg1.t1\n);
    print $fh qq(chr1\tAUGUSTUS\tCDS\t1100\t1500\t0.9\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n);
    print $fh qq(chr1\tAUGUSTUS\tCDS\t1700\t1900\t0.9\t+\t0\ttranscript_id "g1.t1"; gene_id "g1";\n);
    close $fh;
}

# ---------------------------------------------------------------------------
# Case 1: --gff3 produces ID=...;Parent=... attributes and 'mRNA' line type.
# ---------------------------------------------------------------------------
{
    my $out = "$tmp/case1.gff3";
    my $rc = system("perl $script --out $out --gff3 < $gtf");
    is($rc, 0, "case 1: --gff3 exits 0");
    ok(-s $out, "case 1: GFF3 output written");
    my $content = do { open my $r, '<', $out or die $!; local $/; <$r> };
    like($content, qr/ID=g1/, "case 1: gene ID=g1 propagates into GFF3");
    like($content, qr/Parent=g1\.t1|ID=g1\.t1/, "case 1: transcript gets ID/Parent of g1.t1");
    like($content, qr/\tmRNA\t/, "case 1: 'transcript' is rewritten as 'mRNA' in GFF3 mode");
}

# ---------------------------------------------------------------------------
# Case 2: --printExon causes synthetic exon lines to be emitted from the CDS+UTR
# input. Without --printExon (and no input exon lines) the output should
# contain no '\texon\t' rows; with --printExon they appear.
# ---------------------------------------------------------------------------
{
    my $out_no_exon = "$tmp/case2_noexon.gff3";
    my $rc1 = system("perl $script --out $out_no_exon --gff3 < $gtf");
    is($rc1, 0, "case 2 (no --printExon): exits 0");

    my $out_with_exon = "$tmp/case2_withexon.gff3";
    my $rc2 = system("perl $script --out $out_with_exon --gff3 --printExon < $gtf");
    is($rc2, 0, "case 2 (--printExon): exits 0");

    my $with = do { open my $r, '<', $out_with_exon or die $!; local $/; <$r> };
    like($with, qr/\texon\t/, "case 2: --printExon produces 'exon' feature rows");
}
