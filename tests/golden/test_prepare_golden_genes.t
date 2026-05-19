#!/usr/bin/env perl
# Driver-level test for bin/prepare_golden_genes.pl (Phase 4 rewrite).
# Asserts the predictor-output cleanup commitments from plan §4.1:
# v2 emits ONLY Augustus/EVM-shaped artefacts; no SNAP/zff/geneid/glimmer/fathom.

use strict;
use warnings;
use Test::More tests => 9;
use File::Temp qw(tempdir);
use File::Spec;
use Cwd qw(abs_path);
use FindBin qw($Bin);

my $repo     = abs_path("$Bin/../..");
my $sif      = "$repo/containers/jamg.sif";
my $genome   = "$repo/test_suite/mini-genome.fasta";
my $mrna     = "$repo/test_suite/mini-transcripts.fasta";
my $softmask = "$repo/test_suite/output/repeats/mini-genome.fasta.softmasked";

if (! -s $sif) {
    plan skip_all => "containers/jamg.sif missing; run 'make sifs RM_LIB_HOST=...' first";
}
if (! -s $softmask) {
    plan skip_all => "softmasked fixture missing at $softmask; run 'bash tests/rules/test_repeats.sh' first";
}

my $tmp = tempdir(CLEANUP => 0);

# Two execution modes:
#   1. JAMG_GOLDEN_PRESTAGED=<dir> — re-use a prior driver-run output dir.
#      Set this when running the test off the Slurm node where the driver
#      has already been dispatched via `bash tests/golden/run_driver_sbatch.sh`.
#   2. SLURM_JOB_ID set (we're already inside an allocation) — run inline.
# Otherwise the test is skipped: invoking apptainer-exec from `prove` on
# the lazebnik head-node violates the "all heavy compute via Slurm" rule
# and is also painfully slow (~5 min for what takes 55 s on a compute node).
my ($rc, $out);
if (my $pre = $ENV{JAMG_GOLDEN_PRESTAGED}) {
    if (! -d $pre) {
        plan skip_all => "JAMG_GOLDEN_PRESTAGED='$pre' is not a directory";
    }
    $tmp = $pre;
    $rc  = 0;
    $out = '';
}
elsif ($ENV{SLURM_JOB_ID}) {
    my $driver = "$repo/bin/prepare_golden_genes.pl";
    my $cmd = qq{apptainer exec --bind $repo:$repo --bind $tmp:$tmp $sif }
            . qq{perl $driver }
            . qq{-genome $genome -mrna $mrna -softmasked $softmask }
            . qq{-outdir $tmp -threads $ENV{SLURM_CPUS_PER_TASK} 2>&1};
    $out = `$cmd`;
    $rc  = $?;
}
else {
    plan skip_all => "set JAMG_GOLDEN_PRESTAGED=<outdir> (e.g. from "
        . "'bash tests/golden/run_driver_sbatch.sh <outdir>' via sbatch) "
        . "or run inside a SLURM_JOB_ID allocation; refusing to invoke "
        . "apptainer-exec on the head node";
}

is($rc, 0, 'driver exits 0');

ok(-s "$tmp/final_golden_genes.gff3.nr.golden.gff3",
    'golden gff3 produced');

ok(-s "$tmp/final_golden_genes.gff3.nr.golden.train.good.gb",
    'GenBank train.good.gb produced');

unlike($out, qr/\b(snap|zff|geneid|glimmer|fathom)\b/i,
    'no SNAP/zff/geneid/glimmer/fathom references in driver output');

ok(! -e "$tmp/final_golden_genes.gff3.nr.golden.train.good.gb.geneid",
    'no geneid output file emitted');

ok(! -e "$tmp/final_golden_genes.gff3.nr.golden.train.good.gb.glimmer",
    'no GlimmerHMM output file emitted');

ok(! -e "$tmp/final_golden_genes.gff3.nr.golden.zff",
    'no SNAP zff output file emitted');

ok(! -e "$tmp/final_golden_genes.gff3.nr.golden.train.zff",
    'no SNAP train.zff output file emitted');

ok(-s "$tmp/gene_validation.log",
    'gene_validation.log records the pure-Perl validator decisions');

diag("driver output captured at: $tmp") unless $rc == 0;
