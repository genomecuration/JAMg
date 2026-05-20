#!/usr/bin/env perl
# Unit test for bin/splitfasta.pl
# Per POD + GetOptions: actual flag is -depth N (sequences per file),
# -dir <DIR>, -suffix <S>, -i|-fa|-fasta <FILE>.
# (The plan's earlier spec said -size; -depth is what the script actually
# parses. See bin/splitfasta.pl:27-32.)

use strict;
use warnings;
use Test::More tests => 6;
use File::Temp qw(tempdir);
use FindBin qw($Bin);

my $script = "$Bin/../../bin/splitfasta.pl";
my $tmp = tempdir(CLEANUP => 1);

# Build a 5-sequence FASTA fixture.
sub make_fasta {
    my $path = shift;
    open my $fh, '>', $path or die $!;
    for my $name (qw(seqA seqB seqC seqD seqE)) {
        print $fh ">$name short header\nACGTACGT\n";
    }
    close $fh;
}

# ---------------------------------------------------------------------------
# Case 1: -depth 1 -> one sequence per file, exactly 5 files written.
# ---------------------------------------------------------------------------
{
    my $fa  = "$tmp/case1.fa";
    my $dir = "$tmp/case1_dir";
    make_fasta($fa);
    system("perl $script -i $fa -depth 1 -dir $dir") == 0 or die "case1: $?";
    ok(-d $dir, "case 1: output directory created");
    opendir my $dh, $dir or die $!;
    my @files = grep { !/^\.\.?$/ } readdir $dh;
    closedir $dh;
    is(scalar(@files), 5, "case 1: -depth 1 produces 5 per-seq files");
}

# ---------------------------------------------------------------------------
# Case 2: -depth 3 with 5 seqs -> 2 files (seqs 1-3 in file 1, seqs 4-5 in
# file 2 plus the new-file rollover at seq=4). Total headers across all
# output files MUST equal input (5).
# ---------------------------------------------------------------------------
{
    my $fa  = "$tmp/case2.fa";
    my $dir = "$tmp/case2_dir";
    make_fasta($fa);
    system("perl $script -i $fa -depth 3 -dir $dir") == 0 or die "case2: $?";
    opendir my $dh, $dir or die $!;
    my @files = sort grep { !/^\.\.?$/ } readdir $dh;
    closedir $dh;
    ok(scalar(@files) >= 1, "case 2: -depth 3 over 5 seqs produces >= 1 file");
    my $total_headers = 0;
    for my $f (@files) {
        my $content = do { open my $rh, '<', "$dir/$f" or die $!; local $/; <$rh> };
        $total_headers += () = ($content =~ /^>/gm);
    }
    is($total_headers, 5, "case 2: all 5 input sequences emitted across output files");
}

# ---------------------------------------------------------------------------
# Case 3: -suffix .pep -> output filename respects suffix.
# ---------------------------------------------------------------------------
{
    my $fa  = "$tmp/case3.fa";
    my $dir = "$tmp/case3_dir";
    make_fasta($fa);
    system("perl $script -i $fa -depth 1 -suffix .pep -dir $dir") == 0 or die "case3: $?";
    opendir my $dh, $dir or die $!;
    my @files = grep { /\.pep$/ } readdir $dh;
    closedir $dh;
    ok(scalar(@files) >= 1, "case 3: -suffix .pep produces files ending in .pep");
    is(scalar(@files), 5, "case 3: 5 .pep files for 5 sequences");
}
