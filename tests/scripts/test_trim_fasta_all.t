#!/usr/bin/env perl
# Unit test for bin/trim_fasta_all.pl
# Script uses Bio::SeqIO via -fa <FASTA> -length <CUTOFF>. Sequences shorter
# than CUTOFF (in nt; the script does not do aa conversion) go to <FASTA>.discard;
# those >= CUTOFF land in <FASTA>.trim.
#
# (Plan note: the plan-text said "-le 50 keeps seqs >= 50aa" but the script
# is purely nucleic-length; we assert on the script's actual behavior.)

use strict;
use warnings;
use Test::More tests => 8;
use File::Temp qw(tempdir);
use FindBin qw($Bin);

my $script = "$Bin/../../bin/trim_fasta_all.pl";
my $tmp = tempdir(CLEANUP => 1);

# ---------------------------------------------------------------------------
# Build a fasta with 3 short (10 nt) and 2 long (100 nt) sequences.
# ---------------------------------------------------------------------------
my $fa = "$tmp/case.fa";
{
    open my $fh, '>', $fa or die $!;
    for my $name (qw(short1 short2 short3)) {
        print $fh ">$name\nACGTACGTAC\n";  # 10 nt
    }
    for my $name (qw(long1 long2)) {
        print $fh ">$name\n", "ACGT" x 25, "\n";  # 100 nt
    }
    close $fh;
}

# ---------------------------------------------------------------------------
# Case 1: -length 50 trims sequences < 50 nt.
#   .trim file contains the 2 long seqs.
#   .discard file contains the 3 short seqs.
# ---------------------------------------------------------------------------
{
    my $rc = system("cd $tmp && perl $script -fa $fa -length 50 >/dev/null 2>&1");
    is($rc, 0, "case 1: trim_fasta_all -length 50 exits 0");
    my $trim    = "$fa.trim";
    my $discard = "$fa.discard";
    ok(-s $trim, "case 1: .trim file written");

    my $trim_content = do { open my $r, '<', $trim or die $!; local $/; <$r> };
    my @trim_hdrs = ($trim_content =~ /^>(\S+)/gm);
    is(scalar(@trim_hdrs), 2, "case 1: 2 sequences kept (long1, long2)");
    is_deeply(
        [sort @trim_hdrs],
        [sort qw(long1 long2)],
        "case 1: the 2 long sequences survive"
    );

    if (-s $discard) {
        my $disc = do { open my $r, '<', $discard or die $!; local $/; <$r> };
        my @disc_hdrs = ($disc =~ /^>(\S+)/gm);
        is(scalar(@disc_hdrs), 3, "case 1: 3 sequences discarded (short1..3)");
    } else {
        fail("case 1: discard file expected but missing");
    }
}

# ---------------------------------------------------------------------------
# Case 2: -list <FILE> filters by ID. List contains 'long1' -> long1 is
# discarded (the list is the set of IDs to REMOVE per script semantics);
# the surviving .trim file holds long2 plus the 3 shorts that did not match
# any -length cutoff (we omit -length here).
# ---------------------------------------------------------------------------
{
    my $fa2 = "$tmp/case2.fa";
    open my $fh, '>', $fa2 or die $!;
    for my $name (qw(keepA keepB removeMe)) {
        print $fh ">$name\n", "ACGT" x 25, "\n";  # 100 nt each
    }
    close $fh;

    my $list = "$tmp/case2.list";
    open my $lh, '>', $list or die $!;
    print $lh "removeMe\n";
    close $lh;

    my $rc = system("cd $tmp && perl $script -fa $fa2 -list $list >/dev/null 2>&1");
    is($rc, 0, "case 2: trim_fasta_all -list exits 0");
    my $trim = "$fa2.trim";
    my $disc = "$fa2.discard";
    ok(-s $trim, "case 2: .trim file written");
    my $trim_content = do { open my $r, '<', $trim or die $!; local $/; <$r> };
    my @trim_hdrs = ($trim_content =~ /^>(\S+)/gm);
    is_deeply([sort @trim_hdrs], [sort qw(keepA keepB)],
        "case 2: keepA/keepB survive; removeMe is filtered out");
}
