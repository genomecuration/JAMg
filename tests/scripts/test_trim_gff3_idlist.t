#!/usr/bin/env perl
# Unit test for bin/trim_gff3.pl (existing v1 ID-list filter — distinct from
# the new bin/trim_overlap_gff3.py).
# Signature: trim_gff3.pl <LIST> <GFF>
#   - LIST: one ID per line (FASTA headers w/ '>' OK; leading whitespace OK)
#   - GFF: records separated by blank lines; column 0 (seqid) is matched against
#     IDs in LIST. Matches go to stdout, non-matches to stderr.

use strict;
use warnings;
use Test::More tests => 7;
use File::Temp qw(tempdir);
use FindBin qw($Bin);

my $script = "$Bin/../../bin/trim_gff3.pl";
my $tmp = tempdir(CLEANUP => 1);

# ---------------------------------------------------------------------------
# Case 1: LIST contains ctgA -> ctgA records to stdout, ctgB records to stderr.
# ---------------------------------------------------------------------------
{
    my $list = "$tmp/case1.list";
    open my $fh, '>', $list or die $!;
    print $fh "ctgA\n";
    close $fh;

    my $gff = "$tmp/case1.gff3";
    open my $fh2, '>', $gff or die $!;
    print $fh2 "ctgA\tmaker\tgene\t1\t100\t.\t+\t.\tID=gA\n";
    print $fh2 "ctgA\tmaker\tmRNA\t1\t100\t.\t+\t.\tID=mA;Parent=gA\n\n";
    print $fh2 "ctgB\tmaker\tgene\t1\t100\t.\t+\t.\tID=gB\n";
    print $fh2 "ctgB\tmaker\tmRNA\t1\t100\t.\t+\t.\tID=mB;Parent=gB\n\n";
    close $fh2;

    my $stdout_f = "$tmp/case1.stdout";
    my $stderr_f = "$tmp/case1.stderr";
    my $rc = system("perl $script $list $gff > $stdout_f 2> $stderr_f");
    is($rc, 0, "case 1: exit 0");

    my $so = do { open my $r, '<', $stdout_f or die $!; local $/; <$r> };
    my $se = do { open my $r, '<', $stderr_f or die $!; local $/; <$r> };
    like($so, qr/ctgA\tmaker\tgene/, "case 1: ctgA gene record on stdout");
    unlike($so, qr/ctgB\t/, "case 1: ctgB does NOT appear on stdout");
    like($se, qr/ctgB\tmaker\tgene/, "case 1: ctgB gene record on stderr");
}

# ---------------------------------------------------------------------------
# Case 2: LIST file uses FASTA-header form ('>id') — leading '>' is stripped
# during parsing. ID with leading whitespace also accepted.
# ---------------------------------------------------------------------------
{
    my $list = "$tmp/case2.list";
    open my $fh, '>', $list or die $!;
    print $fh ">ctgKEEP\n";       # FASTA-style line
    print $fh "  ctgKEEP2\n";     # leading whitespace
    close $fh;

    my $gff = "$tmp/case2.gff3";
    open my $fh2, '>', $gff or die $!;
    print $fh2 "ctgKEEP\tmaker\tgene\t1\t100\t.\t+\t.\tID=gK\n\n";
    print $fh2 "ctgKEEP2\tmaker\tgene\t200\t300\t.\t+\t.\tID=gK2\n\n";
    print $fh2 "ctgDROP\tmaker\tgene\t1\t100\t.\t+\t.\tID=gD\n\n";
    close $fh2;

    my $stdout_f = "$tmp/case2.stdout";
    my $stderr_f = "$tmp/case2.stderr";
    my $rc = system("perl $script $list $gff > $stdout_f 2> $stderr_f");
    is($rc, 0, "case 2: exit 0");
    my $so = do { open my $r, '<', $stdout_f or die $!; local $/; <$r> };
    like($so, qr/ctgKEEP\t/,  "case 2: '>ctgKEEP' list line matches ctgKEEP records");
    like($so, qr/ctgKEEP2\t/, "case 2: '  ctgKEEP2' (leading whitespace) line matches ctgKEEP2");
}
