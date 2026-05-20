#!/usr/bin/env perl
# Unit test for bin/add_source_gff.pl
# Script signature: add_source_gff.pl <GFF> <SOURCE>
# Behaviour: rewrites column 2 (source) to SOURCE for every line that has
# tab-separated columns; writes result to <GFF>.out. Comment-only lines
# (no tab-split second column) pass through unchanged.

use strict;
use warnings;
use Test::More tests => 6;
use File::Temp qw(tempdir);
use FindBin qw($Bin);

my $script = "$Bin/../../bin/add_source_gff.pl";
ok(-x $script || -f $script, "script exists at $script");

my $tmp = tempdir(CLEANUP => 1);

# ---------------------------------------------------------------------------
# Case 1: every non-comment tabbed line gets column 2 rewritten to "FOO".
# ---------------------------------------------------------------------------
{
    my $infile = "$tmp/case1.gff3";
    open my $fh, '>', $infile or die $!;
    print $fh "##gff-version 3\n";
    print $fh "ctgA\tmaker\tgene\t1\t100\t.\t+\t.\tID=g1\n";
    print $fh "ctgA\tmaker\tmRNA\t1\t100\t.\t+\t.\tID=m1;Parent=g1\n";
    print $fh "ctgA\tmaker\texon\t1\t100\t.\t+\t.\tID=e1;Parent=m1\n";
    print $fh "ctgB\tother\tgene\t200\t300\t.\t-\t.\tID=g2\n";
    close $fh;

    system("perl $script $infile FOO") == 0 or die "script failed";
    my $out = "$infile.out";

    open my $rh, '<', $out or die "no output: $!";
    my @lines = <$rh>;
    close $rh;

    my @noncomment = grep { !/^##/ && /\t/ } @lines;
    is(scalar(@noncomment), 4, "case 1: four data lines emitted");
    my $col2_all_foo = 1;
    for my $ln (@noncomment) {
        my @data = split("\t", $ln);
        $col2_all_foo = 0 if $data[1] ne 'FOO';
    }
    ok($col2_all_foo, "case 1: column 2 rewritten to FOO on every data line");
}

# ---------------------------------------------------------------------------
# Case 2: comment-only line (no tab) passes through unchanged.
# ---------------------------------------------------------------------------
{
    my $infile = "$tmp/case2.gff3";
    open my $fh, '>', $infile or die $!;
    print $fh "##gff-version 3\n";
    print $fh "# a free-form comment with no tabs\n";
    print $fh "ctgX\told_src\tgene\t1\t99\t.\t+\t.\tID=gX\n";
    close $fh;

    system("perl $script $infile NEWSRC") == 0 or die "script failed";
    my $out = "$infile.out";
    open my $rh, '<', $out or die "no output: $!";
    my $content = do { local $/; <$rh> };
    close $rh;

    like($content, qr/^# a free-form comment with no tabs$/m,
        "case 2: untabbed comment line passes through unchanged");
}

# ---------------------------------------------------------------------------
# Case 3: empty input produces an empty (or near-empty) output and exits 0.
# ---------------------------------------------------------------------------
{
    my $infile = "$tmp/case3.gff3";
    open my $fh, '>', $infile or die $!;
    close $fh;  # empty

    my $rc = system("perl $script $infile EMPTYSRC");
    is($rc, 0, "case 3: empty input does not crash the script");
    my $out = "$infile.out";
    ok(-e $out && -z $out, "case 3: empty input produces an empty output file");
}
