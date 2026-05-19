#!/usr/bin/env perl
# Pure-Perl validator unit test for Golden::Filter::validate_gene_structure.
# No SIF needed — runs on the host perl.
#
# 14 cases per plan §4.1b (12 original + ambiguous-codon + frame-shift-stop-
# at-exon-boundary). Synthetic genomes are built from deterministic
# building blocks so each assertion is independently checkable by reading
# the local block.

use strict;
use warnings;
use Test::More tests => 14;
use FindBin qw($Bin);
use lib "$Bin/../../PerlLib";
use Golden::Filter qw(validate_gene_structure _splice_cds _revcomp);

# ---------------------------------------------------------------------------
# Canonical building blocks. Reused across cases so the math stays auditable.
#   CDS1   = ATG AAA CGG       (9 nt; M-K-R in frame)
#   intron = GT (CCCCCCCCCCCCCCC) AG  (19 nt; canonical GT..AG splice)
#   CDS2   = GGG TAA           (6 nt; G-stop in continued frame)
# Plus-strand spliced product = ATGAAACGGGGGTAA (15 nt; M-K-R-G-*).
# ---------------------------------------------------------------------------
my $CDS1   = 'ATGAAACGG';
my $CDS2   = 'GGGTAA';
my $intron = 'GT' . ('C' x 15) . 'AG';
my $pre    = 'N' x 100;   # filler so CDS1 starts at genome position 101
my $post   = 'N' x 10;

# 1. plus-strand canonical → OK
{
    my $genome = $pre . $CDS1 . $intron . $CDS2 . $post;
    my $cds = [
        { type => 'CDS', start => 101, end => 109, strand => '+' },
        { type => 'CDS', start => 129, end => 134, strand => '+' },
    ];
    my ($s, $r) = validate_gene_structure($cds, $genome, { complete => 1 });
    is($s, 'OK', "1: plus-strand canonical two-exon gene validates (reason='$r')");
}

# 2. minus-strand canonical → OK
# Genome carries rc(transcript). Transcript pos T → genome pos (101+34-T).
#   CDS1 (transcript 1-9) → genome 126-134
#   CDS2 (transcript 29-34) → genome 101-106
{
    my $tr = $CDS1 . $intron . $CDS2;
    my $genome = $pre . _revcomp($tr) . $post;
    my $cds = [
        { type => 'CDS', start => 101, end => 106, strand => '-' },
        { type => 'CDS', start => 126, end => 134, strand => '-' },
    ];
    my ($s, $r) = validate_gene_structure($cds, $genome, { complete => 1 });
    is($s, 'OK', "2: minus-strand canonical two-exon gene validates (reason='$r')");
}

# 3. missing start (GTG, not ATG) → FAIL: missing_start
{
    my $cds1_bad = 'GTGAAACGG';
    my $genome = $pre . $cds1_bad . $intron . $CDS2 . $post;
    my $cds = [
        { type => 'CDS', start => 101, end => 109, strand => '+' },
        { type => 'CDS', start => 129, end => 134, strand => '+' },
    ];
    my ($s, $r) = validate_gene_structure($cds, $genome, { complete => 1 });
    is($r, 'missing_start', "3: missing-start rejected (status=$s)");
}

# 4. missing stop (GGG CCC instead of GGG TAA) → FAIL: missing_stop
{
    my $cds2_bad = 'GGGCCC';
    my $genome = $pre . $CDS1 . $intron . $cds2_bad . $post;
    my $cds = [
        { type => 'CDS', start => 101, end => 109, strand => '+' },
        { type => 'CDS', start => 129, end => 134, strand => '+' },
    ];
    my ($s, $r) = validate_gene_structure($cds, $genome, { complete => 1 });
    is($r, 'missing_stop', "4: missing-stop rejected (status=$s)");
}

# 5. non-canonical splice (AT..AC instead of GT..AG) → FAIL: non_canonical_splice:AT-AC
{
    my $intron_bad = 'AT' . ('C' x 15) . 'AC';
    my $genome = $pre . $CDS1 . $intron_bad . $CDS2 . $post;
    my $cds = [
        { type => 'CDS', start => 101, end => 109, strand => '+' },
        { type => 'CDS', start => 129, end => 134, strand => '+' },
    ];
    my ($s, $r) = validate_gene_structure($cds, $genome);
    like($r, qr/^non_canonical_splice/, "5: non-canonical splice rejected (reason='$r')");
}

# 6. internal stop in single-exon CDS → FAIL: internal_stop_count_1
{
    my $cds_seq = 'ATGTAAGGGAAATAA';   # 15 nt; M-*-G-K-*; 1 internal stop
    my $genome = $pre . $cds_seq . $post;
    my $cds = [
        { type => 'CDS', start => 101, end => 115, strand => '+' },
    ];
    my ($s, $r) = validate_gene_structure($cds, $genome);
    is($r, 'internal_stop_count_1', "6: single-exon internal stop rejected (status=$s)");
}

# 7. length not divisible by 3 → FAIL: cds_length_not_mod3
{
    my $cds_seq = 'ATGAAAGGGT';        # 10 nt — not mod 3
    my $genome = $pre . $cds_seq . $post;
    my $cds = [
        { type => 'CDS', start => 101, end => 110, strand => '+' },
    ];
    my ($s, $r) = validate_gene_structure($cds, $genome);
    is($r, 'cds_length_not_mod3', "7: non-mod3 CDS length rejected (status=$s)");
}

# 8. empty CDS array → FAIL: no_CDS
{
    my ($s, $r) = validate_gene_structure([], $pre);
    is($r, 'no_CDS', "8: empty CDS rejected (status=$s)");
}

# 9. single-exon plus-strand → OK
{
    my $cds_seq = 'ATGAAATAA';         # 9 nt; M-K-*
    my $genome = $pre . $cds_seq . $post;
    my $cds = [
        { type => 'CDS', start => 101, end => 109, strand => '+' },
    ];
    my ($s, $r) = validate_gene_structure($cds, $genome, { complete => 1 });
    is($s, 'OK', "9: single-exon plus-strand validates (reason='$r')");
}

# 10. single-exon minus-strand → OK
{
    my $cds_seq = 'ATGAAATAA';
    my $genome = $pre . _revcomp($cds_seq) . $post;
    my $cds = [
        { type => 'CDS', start => 101, end => 109, strand => '-' },
    ];
    my ($s, $r) = validate_gene_structure($cds, $genome, { complete => 1 });
    is($s, 'OK', "10: single-exon minus-strand validates (reason='$r')");
}

# 11. adjacent CDS pieces (no real intron, intron_end < intron_start) → OK
# Splice check sees intron_start=106, intron_end=105; skipped via `next if`.
{
    my $cds_seq = 'ATGAAACGGGGGTAA';   # 15 nt; M-K-R-G-*
    my $genome = $pre . $cds_seq . $post;
    my $cds = [
        { type => 'CDS', start => 101, end => 105, strand => '+' },   # ATGAA
        { type => 'CDS', start => 106, end => 115, strand => '+' },   # ACGGGGGTAA
    ];
    my ($s, $r) = validate_gene_structure($cds, $genome, { complete => 1 });
    is($s, 'OK', "11: adjacent CDS pieces validate (reason='$r')");
}

# 12. mixed-strand within one gene → die (caller bug)
{
    my $cds_seq = 'ATGAAACGGGGGTAA';
    my $genome = $pre . $cds_seq . $post;
    my $cds = [
        { type => 'CDS', start => 101, end => 105, strand => '+' },
        { type => 'CDS', start => 106, end => 115, strand => '-' },
    ];
    eval { validate_gene_structure($cds, $genome) };
    like($@, qr/mixed strand/i,
        "12: mixed-strand input dies with informative message ($@)");
}

# 13. ambiguous codon (NNN in CDS) → FAIL: ambiguous_codon_count_1
{
    my $cds_seq = 'ATGNNNAAATAA';      # 12 nt; M-X-K-* (1 X codon)
    my $genome = $pre . $cds_seq . $post;
    my $cds = [
        { type => 'CDS', start => 101, end => 112, strand => '+' },
    ];
    my ($s, $r) = validate_gene_structure($cds, $genome);
    is($r, 'ambiguous_codon_count_1', "13: ambiguous codon rejected (status=$s)");
}

# 14. frame-shift-induced stop at exon boundary → FAIL: internal_stop_count_>=1
# CDS1='ATGAAA' (M-K, 6 nt), CDS2='TAAGGGTAA' (*-G-*, 9 nt) in continued frame.
# Spliced='ATGAAATAAGGGTAA' → M-K-*-G-*; internal stop arises only after splice.
{
    my $cds1 = 'ATGAAA';
    my $cds2 = 'TAAGGGTAA';
    my $genome = $pre . $cds1 . $intron . $cds2 . $post;
    my $cds = [
        { type => 'CDS', start => 101, end => 106, strand => '+' },
        { type => 'CDS', start => 126, end => 134, strand => '+' },
    ];
    my ($s, $r) = validate_gene_structure($cds, $genome);
    like($r, qr/^internal_stop_count_\d+$/,
        "14: frame-shift-induced stop at exon boundary rejected (reason='$r')");
}
