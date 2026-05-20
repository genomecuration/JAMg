#!/usr/bin/env perl
# Unit test for bin/gff2hints.pl
# Script signature: gff2hints.pl <GFF> [is_golden]
# Behaviour: emits <GFF>.hints in Augustus hints format.
#   - second arg truthy => src=GLD pri=7, no suffix on exon/CDS/intron
#   - second arg falsy  => src=XNT pri=5, "part" suffix appended

use strict;
use warnings;
use Test::More tests => 7;
use File::Temp qw(tempdir);
use FindBin qw($Bin);

my $script = "$Bin/../../bin/gff2hints.pl";

# Build a minimal GFF3 with gene/mRNA/exon/CDS/intron, separated by blank line
# (the script reads records using `\n\n` as $/).
my $gff_template = <<'EOF';
ctgA	maker	gene	100	300	.	+	.	ID=g1
ctgA	maker	mRNA	100	300	.	+	.	ID=m1;Parent=g1
ctgA	maker	exon	100	200	.	+	.	ID=e1;Parent=m1
ctgA	maker	intron	201	250	.	+	.	ID=i1;Parent=m1
ctgA	maker	CDS	100	200	.	+	0	ID=c1;Parent=m1

EOF

my $tmp = tempdir(CLEANUP => 1);

# ---------------------------------------------------------------------------
# Case 1: non-golden mode -> src=XNT, exonpart/CDSpart/intronpart
# ---------------------------------------------------------------------------
{
    my $infile = "$tmp/case1.gff3";
    open my $fh, '>', $infile or die $!;
    print $fh $gff_template;
    close $fh;
    system("perl $script $infile") == 0 or die "case1 script failed";
    my $hints = "$infile.hints";
    ok(-s $hints, "case 1: hints file produced");
    my $content = do { open my $rh, '<', $hints or die $!; local $/; <$rh> };
    like($content, qr/\tsrc=XNT;pri=5/, "case 1: src=XNT and priority 5");
    like($content, qr/\b(exonpart|CDSpart|intronpart)\b/,
        "case 1: feature types get 'part' suffix in non-golden mode");
}

# ---------------------------------------------------------------------------
# Case 2: golden mode -> src=GLD, no 'part' suffix on exon/CDS/intron
# ---------------------------------------------------------------------------
{
    my $infile = "$tmp/case2.gff3";
    open my $fh, '>', $infile or die $!;
    print $fh $gff_template;
    close $fh;
    system("perl $script $infile 1") == 0 or die "case2 script failed";
    my $hints = "$infile.hints";
    my $content = do { open my $rh, '<', $hints or die $!; local $/; <$rh> };
    like($content, qr/\tsrc=GLD;pri=7/, "case 2: src=GLD and priority 7");
    unlike($content, qr/\b(exonpart|CDSpart|intronpart)\b/,
        "case 2: no 'part' suffix in golden mode");
}

# ---------------------------------------------------------------------------
# Case 3: cDNA_match records take the exonpart branch with src=E;pri=3 and a
# grp= attribute pulled from the cDNA_match ID. This is a distinct code path
# from the regular gene/mRNA/exon flow.
# ---------------------------------------------------------------------------
{
    my $infile = "$tmp/case3.gff3";
    open my $fh, '>', $infile or die $!;
    print $fh "ctgZ\tgmap\tcDNA_match\t1000\t1500\t.\t+\t.\tID=cdna_42;Name=cdna42\n";
    print $fh "ctgZ\tgmap\tcDNA_match\t1700\t2000\t.\t+\t.\tID=cdna_42;Name=cdna42\n\n";
    close $fh;
    system("perl $script $infile") == 0 or die "case3 script failed";
    my $content = do { open my $r, '<', "$infile.hints" or die $!; local $/; <$r> };
    like($content, qr/\texonpart\t/, "case 3: cDNA_match emits exonpart");
    like($content, qr/src=E;pri=3/, "case 3: cDNA_match attribute src=E;pri=3");
}
