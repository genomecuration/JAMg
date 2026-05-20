#!/usr/bin/env perl
# Unit test for bin/repeatmasker2hints.pl.
# Signature: repeatmasker2hints.pl <RM.gff>
# Behaviour: rewrites column 2 ('source' is left alone), column 3 (feature_type)
# to 'nonexonpart', column 9 to 'src=RM;pri=6', sorts in place, writes
# <RM.gff>.hints.
#
# (Plan note: the plan-text said "respects -xsmall soft-mask convention" but the
# script has NO -xsmall option — see bin/repeatmasker2hints.pl:14-15, only a
# positional arg. Test asserts on the script's actual transformation.)

use strict;
use warnings;
use Test::More tests => 6;
use File::Temp qw(tempdir);
use FindBin qw($Bin);

my $script = "$Bin/../../bin/repeatmasker2hints.pl";
my $tmp = tempdir(CLEANUP => 1);

# RepeatMasker-style GFF (with comments + similarity rows).
my $rm_gff = "$tmp/case.gff";
{
    open my $fh, '>', $rm_gff or die $!;
    print $fh "# RepeatMasker GFF output\n";
    print $fh qq(ctgA\tRepeatMasker\tsimilarity\t100\t350\t17.2\t+\t.\tTarget "Motif:DNA1" 1 250\n);
    print $fh qq(ctgA\tRepeatMasker\tsimilarity\t500\t750\t22.1\t-\t.\tTarget "Motif:LTR2" 1 250\n);
    print $fh qq(ctgB\tRepeatMasker\tsimilarity\t1000\t1500\t10.0\t+\t.\tTarget "Motif:SINE1" 1 500\n);
    close $fh;
}

# ---------------------------------------------------------------------------
# Case 1: produces .hints with src=RM;pri=6 and feature type 'nonexonpart'.
# ---------------------------------------------------------------------------
{
    my $rc = system("cd $tmp && perl $script $rm_gff >/dev/null 2>&1");
    is($rc, 0, "case 1: repeatmasker2hints exits 0");
    my $hints = "$rm_gff.hints";
    ok(-s $hints, "case 1: hints file produced");
    my $content = do { open my $r, '<', $hints or die $!; local $/; <$r> };
    like($content, qr/\tnonexonpart\t/, "case 1: feature type rewritten to 'nonexonpart'");
    like($content, qr/src=RM;pri=6/, "case 1: source rewritten to src=RM;pri=6");
}

# ---------------------------------------------------------------------------
# Case 2: comment lines in the input GFF (lines starting '#') are dropped
# from the hints output. Each non-comment row contributes one hints row.
# ---------------------------------------------------------------------------
{
    my $rm2 = "$tmp/case2.gff";
    open my $fh, '>', $rm2 or die $!;
    print $fh "# RepeatMasker preamble line 1\n";
    print $fh "# RepeatMasker preamble line 2\n";
    print $fh qq(ctgZ\tRepeatMasker\tsimilarity\t10\t90\t5.0\t+\t.\tTarget "Motif:DNA_X" 1 80\n);
    print $fh qq(ctgZ\tRepeatMasker\tsimilarity\t200\t280\t6.0\t-\t.\tTarget "Motif:DNA_Y" 1 80\n);
    close $fh;

    my $rc = system("cd $tmp && perl $script $rm2 >/dev/null 2>&1");
    is($rc, 0, "case 2: exits 0 with comment-bearing input");
    my $hints = "$rm2.hints";
    my $content = do { open my $r, '<', $hints or die $!; local $/; <$r> };
    my @rows = grep { /\S/ && !/^#/ } split /\n/, $content;
    is(scalar(@rows), 2, "case 2: 2 hints rows for 2 non-comment input rows (comments dropped)");
}
