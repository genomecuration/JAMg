#!/usr/bin/env perl
# Unit test for bin/fix_pasa_evm_gff.pl.
# Signature: fix_pasa_evm_gff.pl <infile> [do_alias]
# Behaviour: splits a PASA-EVM GFF3 into per-status outfiles named
# <infile>.<state>.gff3 (state = 'undefined' if no status= present, plus
# 'weird' bucket for malformed records, plus any state captured from a
# status= attribute on a gene/mRNA line).

use strict;
use warnings;
use Test::More tests => 6;
use File::Temp qw(tempdir);
use FindBin qw($Bin);

my $script = "$Bin/../../bin/fix_pasa_evm_gff.pl";
my $tmp = tempdir(CLEANUP => 1);

# ---------------------------------------------------------------------------
# Build a 2-gene PASA-EVM-style snippet. The script reads ### / blank-line
# as the gene-record delimiter and expects the first line of each record to
# carry "\tgene\t" in column 3.
# ---------------------------------------------------------------------------
my $gff = "$tmp/case.gff";
{
    open my $fh, '>', $gff or die $!;
    print $fh "##gff-version 3\n";
    print $fh "ctgA\tEVM\tgene\t100\t300\t.\t+\t.\tID=g1;Name=g1;status=Finished\n";
    print $fh "ctgA\tEVM\tmRNA\t100\t300\t.\t+\t.\tID=m1;Parent=g1;status=Finished\n";
    print $fh "ctgA\tEVM\texon\t100\t300\t.\t+\t.\tID=e1;Parent=m1\n";
    print $fh "ctgA\tEVM\tCDS\t100\t300\t.\t+\t0\tID=c1;Parent=m1\n";
    print $fh "###\n";
    print $fh "ctgB\tEVM\tgene\t500\t800\t.\t-\t.\tID=g2;Name=g2;status=Inferred\n";
    print $fh "ctgB\tEVM\tmRNA\t500\t800\t.\t-\t.\tID=m2;Parent=g2;status=Inferred\n";
    print $fh "ctgB\tEVM\texon\t500\t800\t.\t-\t.\tID=e2;Parent=m2\n";
    print $fh "ctgB\tEVM\tCDS\t500\t800\t.\t-\t0\tID=c2;Parent=m2\n";
    print $fh "###\n";
    close $fh;
}

# ---------------------------------------------------------------------------
# Case 1: script splits into one or more <infile>.<state>.gff3 files.
# Specifically, status=Finished (g1) and status=Inferred (g2) should be in
# distinct files.
# ---------------------------------------------------------------------------
{
    my $rc = system("perl $script $gff >/dev/null 2>&1");
    is($rc, 0, "case 1: exits 0");
    my $finished_f = "$gff.Finished.gff3";
    my $inferred_f = "$gff.Inferred.gff3";
    ok(-s $finished_f, "case 1: Finished state file produced");
    ok(-s $inferred_f, "case 1: Inferred state file produced");

    my $f = do { open my $r, '<', $finished_f or die $!; local $/; <$r> };
    like($f, qr/ID=g1/, "case 1: Finished file contains g1");
}

# ---------------------------------------------------------------------------
# Case 2: a gene with no status= attribute falls into the 'undefined' bucket.
# ---------------------------------------------------------------------------
{
    my $gff2 = "$tmp/case2.gff";
    open my $fh, '>', $gff2 or die $!;
    print $fh "##gff-version 3\n";
    print $fh "ctgX\tEVM\tgene\t10\t90\t.\t+\t.\tID=gNoStatus;Name=gNoStatus\n";
    print $fh "ctgX\tEVM\tmRNA\t10\t90\t.\t+\t.\tID=mNoStatus;Parent=gNoStatus\n";
    print $fh "ctgX\tEVM\texon\t10\t90\t.\t+\t.\tID=eNoStatus;Parent=mNoStatus\n";
    print $fh "ctgX\tEVM\tCDS\t10\t90\t.\t+\t0\tID=cNoStatus;Parent=mNoStatus\n";
    print $fh "###\n";
    close $fh;

    my $rc = system("perl $script $gff2 >/dev/null 2>&1");
    is($rc, 0, "case 2: exits 0");
    my $undef_f = "$gff2.undefined.gff3";
    ok(-s $undef_f, "case 2: status-less gene routed to 'undefined' bucket file");
}
