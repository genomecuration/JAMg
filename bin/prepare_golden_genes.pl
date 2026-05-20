#!/usr/bin/env perl

=pod

=head1 NAME

 prepare_golden_genes.pl

 Use an exonerate result, PASA assemblies, peptides, or mRNAs to prepare a
 high-confidence "golden" gene set for training Augustus and seeding
 EvidenceModeler. The script aligns the input to the genome (exonerate +
 GMAP), filters for canonical splicing + completeness, partitions into
 training / test / optimisation subsets, and emits the Augustus-shaped
 artefacts (GenBank files + extrinsic hints).

 Emits Augustus-shaped artefacts only (GenBank training files + extrinsic
 hints). The canonical-splice + CDS-completeness check is implemented in
 pure Perl by PerlLib/Golden/Filter.pm.

=head1 USAGE

Mandatory:

  -genome       :s    Genome FASTA. Repeat-masked with Ns if -softmasked
                      is also given.

And one of:

1. PASA outputs from pasa_asmbls_to_training_set.dbi:

    -pasa_assembly :s   PASA contig FASTA
    -pasa_gff      :s   PASA-emitted GFF
    -pasa_peptides :s   PASA-emitted protein FASTA
    -pasa_cds      :s   PASA-emitted CDS FASTA
    -pasa_genome   :s   PASA genome-alignment GFF

2. -peptides   :s   Protein FASTA

3. -mrna       :s   cDNA FASTA

4. -exonerate  :s   Existing exonerate result file

Other options:

    -outdir          :s   Output directory (default: current working dir).
                          chdir into it before all subsequent work.
    -softmasked      :s   Softmasked genome (recommended)
    -gmap_dir        :s   Where the GMAP databases live
                          (default: <JAMG_PATH>/databases/gmap)
    -training        :i   Number of genes for the random training set
                          (default: min(4000, 40% of total))
    -complete             Only accept full-length CDS
    -nosingle             Drop single-exon genes
    -flanks          :i   Bases left+right when building the GenBank file
    -identical       :i   Identity cutoff for an exonerate alignment
                          (default: 95)
    -similar         :i   Similarity cutoff (default: 98)
    -intron          :i   Max intron length (default: 70000)
    -threads|cpu     :i   CPU threads (default: 4)
    -minorf          :i   Minimum ORF length in bp (default: 290)
    -mismatch_cutoff :i   Max mismatches allowed (default: 10)
    -same_species         Input is from the same species as the genome
    -liberal              Use looser cut-offs
    -norefine             Skip exonerate's --refine pass
    -exhaustive           Run exonerate exhaustively (slow)
    -norerun              Trust an existing exonerate result file
    -augustus        :s   Path to the Augustus install dir
                          (used to find the auxiliary scripts)
    -extra_gff       :s   Extra GFF lines to merge in
    -build_only           Build genome indexes and exit
    -help                 Show this manual

=head1 DESCRIPTION

The exonerate post-processing uses this --ryo:

 --ryo "RYOAP_START\nRYOAP_STATS_D\tAlignment length\tidentical\tsimilar\tmismatch\tidentical%%\tsimilar%%\nRYOAP_STATS\t%s\t%et\t%ei\t%es\t%em\t%pi\t%ps\nRYOAP_CODING_QUERY\t%qi\t%qab\t%qae\n>%qi\n%qas\nRYOAP_CODING_GENOME\t%ti\t%tcb\t%tce\n>%ti\n%tcs\nRYOAP_END\n"

Acceptance criteria (canonical mode):

    Methionine + * in protein sequence
    Start + stop codon in genome
    Identity > -identical%
    Similarity > -similar%
    Mismatches <= -mismatch_cutoff
    No overlapping genes (best score wins)
    Canonical splice sites (GT..AG or GC..AG); no isoforms

Golden::Filter::validate_gene_structure (pure Perl) writes one line per
rejected gene to gene_validation.log in -outdir.

=head1 AUTHORS

 Alexie Papanicolaou (2012-) — CSIRO Ecosystem Sciences; phase-4 modular
 rewrite 2026-05.

=head1 LICENSE

 Copyright 2012- the Commonwealth Scientific and Industrial Research
 Organization. See LICENSE.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use POSIX qw(ceil);
use FindBin qw($RealBin);
use File::Basename;
use Cwd qw(abs_path getcwd);

use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";

use threads;
use Thread_helper;
use List::Util 'shuffle';

# JAMg PerlLib modules used directly by the driver
use Gene_obj;
use Gene_obj_indexer;
use CdbTools;
use GFF3_utils;
use GTF_utils;
use Nuc_translator;
use Fasta_reader;
use Carp;

# Phase-4 modules. All three reference the driver's globals via $main::var.
use Golden::Filter    qw(validate_gene_structure filter_gff sort_gff3
                         get_gff_delimiter order_fasta process_cmd
                         check_program check_program_optional check_sort_version
                         read_fasta gff3_fix_phase remove_overlapping_gff
                         parse_genome_gff bed_to_gff3 partition_genes);
use Golden::Alignment qw(prepare_pasa_output run_gmap run_blast run_exonerate
                         correct_exonerate_gff
                         create_golden_gffs_from_exonerate
                         create_golden_gffs_from_gff);
use Golden::Augustus  qw(check_augustus parse_gb gb2gff3 gff2hints
                         shuffle_fasta split_fasta wrap_text count_seq
                         gff_to_gtf);

# ---------------------------------------------------------------------------
# Globals (declared `our` so Golden::* modules can reference them as
# $main::var). Set by GetOptions below.
# ---------------------------------------------------------------------------
our (
    $debug, $verbose, $contig_file, $exonerate_file,
    $genome_file, $is_cdna, $no_single_exon,
    $overwrite, $only_complete, $peptide_file,
    $pasa_gff, $pasa_cds,
    $pasa_assembly_file, $pasa_peptides, $mrna_file,
    $softmasked_genome, $stop_after_correction, $norefine,
    $nodataprint, $no_gmap, $no_exonerate,
    $pasa_genome_gff, $extra_gff_file, $show_help,
    $liberal_cutoffs, $aat_dir, $parafly_exec, $do_exhaustive,
    $filter_exec, $build_only, $outdir,
);

our $no_rerun_exonerate;
our $force_blast;
our $same_species             = '';
our $intron_size              = 70000;
our $training_set_size        = 4000;
our $aug_optimization_geneset = 450;
our $augustus_flank_region    = 3000;
our $minorf                   = 290;
our $identical_fraction_cutoff = 95;
our $similar_fraction_cutoff   = 98;
our $mismatch_cutoff           = 10;
our $threads                   = 4;
our $cwd                       = getcwd() . '/';

our $failed_cutoff = 0;
our ($cdbfasta_exec, $cdbyank_exec);
our (%get_id_seq_from_fasta_hash, $augustus_dir);
our $gmap_dir = "$RealBin/../databases/gmap/";

pod2usage $! unless GetOptions(
    'gmap_dir:s'          => \$gmap_dir,
    'force_blast'         => \$force_blast,
    'help'                => \$show_help,
    'debug'               => \$debug,
    'verbose'             => \$verbose,
    'exonerate:s'         => \$exonerate_file,
    'genome:s'            => \$genome_file,
    'softmasked:s'        => \$softmasked_genome,
    'pasa_gff:s'          => \$pasa_gff,
    'pasa_peptides:s'     => \$pasa_peptides,
    'pasa_cds:s'          => \$pasa_cds,
    'pasa_assembly:s'     => \$pasa_assembly_file,
    'pasa_genome:s'       => \$pasa_genome_gff,
    'peptides:s'          => \$peptide_file,
    'mrnas|mrna:s'        => \$mrna_file,
    'cdna'                => \$is_cdna,
    'training:i'          => \$training_set_size,
    'flanks:i'            => \$augustus_flank_region,
    'identical:i'         => \$identical_fraction_cutoff,
    'similar:i'           => \$similar_fraction_cutoff,
    'overwrite'           => \$overwrite,
    'intron:i'            => \$intron_size,
    'threads|cpus:i'      => \$threads,
    'nosingle|no_single'  => \$no_single_exon,
    'complete'            => \$only_complete,
    'stop_exonerate'      => \$stop_after_correction,
    'minorf:i'            => \$minorf,
    'mismatch_cutoff:i'   => \$mismatch_cutoff,
    'same_species'        => \$same_species,
    'norefine'            => \$norefine,
    'exhaustive'          => \$do_exhaustive,
    'norerun'             => \$no_rerun_exonerate,
    'nodataprint'         => \$nodataprint,
    'augustus|augustus_dir:s' => \$augustus_dir,
    'no_exonerate'        => \$no_exonerate,
    'no_gmap'             => \$no_gmap,
    'extra_gff:s'         => \$extra_gff_file,
    'liberal'             => \$liberal_cutoffs,
    'build_only'          => \$build_only,
    'outdir:s'            => \$outdir,
);
pod2usage if $show_help;

# chdir into the output directory before all subsequent work so every
# emitted artefact lands under $outdir (Phase 4 -outdir contract).
$outdir //= '.';
unless (-d $outdir) {
    mkdir $outdir or die "Cannot create -outdir '$outdir': $!\n";
}
# Resolve inputs against the original cwd BEFORE we chdir; absolutise so
# downstream sub calls work regardless of cwd.
for my $ref ( \$genome_file, \$softmasked_genome, \$pasa_gff, \$pasa_cds,
              \$pasa_assembly_file, \$pasa_peptides, \$pasa_genome_gff,
              \$peptide_file, \$mrna_file, \$exonerate_file, \$extra_gff_file ) {
    $$ref = abs_path($$ref) if defined $$ref && length $$ref && -e $$ref;
}
$gmap_dir = abs_path($gmap_dir) if defined $gmap_dir && -d $gmap_dir;
chdir $outdir or die "chdir $outdir: $!\n";
$cwd = getcwd() . '/';

our $sort_buffer = '5G';
our $tmpdir = $ENV{TMP} // $ENV{TMPDIR}
    // die "TMP or TMPDIR env var must be set; /tmp is forbidden on this host\n";
our $sort_exec = Golden::Filter::check_sort_version();

# Required programs
($cdbfasta_exec, $cdbyank_exec) = Golden::Filter::check_program('cdbfasta', 'cdbyank');
our ($makeblastdb_exec, $tblastn_exec, $tblastx_exec)
    = Golden::Filter::check_program('makeblastdb', 'tblastn', 'tblastx');
our ($gmap_build_exec, $gmap_exec, $gff3_introns_exec)
    = Golden::Filter::check_program('gmap_build', 'gmap', 'gff3_introns');
our ($samtools_exec) = Golden::Filter::check_program('samtools');

die "No genome found in '$genome_file'\n" if $genome_file && !-s $genome_file;
pod2usage "No genome found\n" unless ( $genome_file && -s $genome_file );

# Auto-pick gmapl for big genomes (>~4 GB).
if ( $genome_file && -s $genome_file > 4306887543 ) {
    ($gmap_exec) = Golden::Filter::check_program('gmapl');
}

our ( $gff2gb_exec, $augustus_exec, $augustus_train_exec, $augustus_filterGenes_exec )
    = Golden::Filter::check_program_optional(
        'gff2gbSmallDNA.pl', 'augustus', 'etraining', 'filterGenes.pl'
    );

# Plain genome-sequence file is the (possibly softmasked) one we feed BLAST.
our $genome_sequence_file =
    $softmasked_genome ? $softmasked_genome : $genome_file;
our $genome_sequence_file_dir  = dirname($genome_sequence_file);
our $genome_sequence_file_base = basename($genome_sequence_file);
our $genome_dir                = basename($genome_file) . '_dir';
# run_blast / run_aat / parse_genome_gff each write per-scaffold filter
# files into $genome_dir. mkdir up-front rather than per-caller.
mkdir $genome_dir unless -d $genome_dir;
our $genome_dbname             = $genome_sequence_file_base . '.gmap';
print "Processing genome $genome_sequence_file_dir/$genome_sequence_file_base\n";

our ($scaffold_seq_hashref, $scaffold_seq_length)
    = Golden::Filter::read_fasta($genome_sequence_file);

print "Building BLAST directory at $genome_sequence_file\n";
Golden::Filter::process_cmd(
    "$makeblastdb_exec -in $genome_sequence_file -out $genome_sequence_file -hash_index -parse_seqids -dbtype nucl"
) unless ( -s "$genome_sequence_file.nin" || -s "$genome_sequence_file.nal" );

print "Building GMAP directory at $gmap_dir/$genome_dbname\n";
Golden::Filter::process_cmd(
    "$gmap_build_exec -e 0 -D $gmap_dir -d $genome_dbname $genome_file > /dev/null"
) unless ( -d "$gmap_dir/$genome_dbname" || ($peptide_file && !$build_only) );

print "Indexing...\n";
Golden::Filter::process_cmd("$samtools_exec faidx $genome_sequence_file")
    unless -s "$genome_sequence_file.fai";
Golden::Filter::process_cmd("$cdbfasta_exec $genome_sequence_file")
    unless -s "$genome_sequence_file.cidx";

exit 0 if $build_only;

check_for_options();

if ($liberal_cutoffs) {
    print "\nLiberal cut-offs requested\n";
    $identical_fraction_cutoff = 40   if $identical_fraction_cutoff == 95;
    $similar_fraction_cutoff   = 50   if $similar_fraction_cutoff   == 98;
    $mismatch_cutoff           = 1000 if $mismatch_cutoff           == 10;
    print "\nMethionine, stop codons, splice sites and overlapping genes will not be imposed\n";
}
else {
    print "\nMethionine + * in protein sequence; canonical splice sites (GT..AG or GC..AG); no overlapping genes\n";
}
print "These cutoffs will be used:
 -identical_fraction : $identical_fraction_cutoff
 -similar_fraction   : $similar_fraction_cutoff
 -mismatch           : $mismatch_cutoff
\n";

# ---------------------------------------------------------------------------
# Main flow
# ---------------------------------------------------------------------------
my ( $gmap_gff, $gmap_passed, $exonerate_gff, $exonerate_passed,
     %passed_check, @to_evaluate );

unless ($no_gmap) {
    if ( $mrna_file && -s $mrna_file ) {
        ($gmap_gff, $gmap_passed) = Golden::Alignment::run_gmap($mrna_file);
    }
    elsif ( $pasa_cds && -s $pasa_cds ) {
        ($gmap_gff, $gmap_passed) = Golden::Alignment::run_gmap($pasa_cds);
    }
}
($exonerate_gff, $exonerate_passed) = Golden::Alignment::run_exonerate()
    unless $no_exonerate;

my $final_gff = 'final_golden_genes.gff3';
open( my $out1, '>', $final_gff )         or die "Cannot write $final_gff: $!\n";
open( my $out2, '>', "$final_gff.passed" ) or die "Cannot write $final_gff.passed: $!\n";

if ( $exonerate_gff && -s $exonerate_gff ) {
    open( my $in, '<', $exonerate_gff ) or die "Cannot open $exonerate_gff: $!";
    push @to_evaluate, $exonerate_gff;
    while (<$in>) { print $out1 $_ }
    close $in;
    open( $in, '<', $exonerate_passed ) or die "Cannot open $exonerate_passed: $!";
    while (my $ln = <$in>) {
        if ($ln =~ /^(\S+)/) {
            next if $passed_check{$1};
            $passed_check{$1}++;
            print $out2 $ln;
        }
    }
    close $in;
}
if ( $gmap_gff && -s $gmap_gff ) {
    open( my $in, '<', $gmap_gff ) or die "Cannot open $gmap_gff: $!";
    push @to_evaluate, $gmap_gff;
    while (<$in>) { print $out1 $_ }
    close $in;
    open( $in, '<', $gmap_passed ) or die "Cannot open $gmap_passed: $!";
    while (my $ln = <$in>) {
        if ($ln =~ /^(\S+)/) {
            next if $passed_check{$1};
            $passed_check{$1}++;
            print $out2 $ln;
        }
    }
    close $in;
}
if ( $extra_gff_file && -s $extra_gff_file ) {
    print "Processing user-provided GFF3 $extra_gff_file\n";
    Golden::Alignment::create_golden_gffs_from_gff($extra_gff_file, "$extra_gff_file.n");
    Golden::Filter::remove_overlapping_gff("$extra_gff_file.n", "$extra_gff_file.nr");
    $extra_gff_file .= '.nr';
    Golden::Alignment::create_golden_gffs_from_gff($extra_gff_file);
    push @to_evaluate, "$extra_gff_file.golden";
    open( my $in, '<', "$extra_gff_file.golden" )
        or die "Cannot open $extra_gff_file.golden: $!";
    while (<$in>) { print $out1 $_ }
    close $in;
    open( $in, '<', "${extra_gff_file}passed" )
        or die "Cannot open ${extra_gff_file}passed: $!";
    while (my $ln = <$in>) {
        if ($ln =~ /^(\S+)/) {
            next if $passed_check{$1};
            $passed_check{$1}++;
            print $out2 $ln;
        }
    }
    close $in;
}
close $out1;
close $out2;

if ( $final_gff && -s $final_gff ) {
    Golden::Filter::remove_overlapping_gff(
        $final_gff, "$final_gff.nr", undef, "$final_gff.passed"
    );
    Golden::Filter::partition_genes(
        "$final_gff.nr", "$final_gff.passed", $training_set_size,
        $augustus_flank_region, $augustus_dir, $genome_sequence_file,
        $cdbyank_exec, $gff3_introns_exec, '.', $only_complete,
    );
    print "Done!\n";
}
else {
    warn "There was an issue producing the final GFF file. Exiting.\n";
    exit 1;
}

exit 0;

# ---------------------------------------------------------------------------
# Driver-local helper: check_for_options
# ---------------------------------------------------------------------------
sub check_for_options {
    pod2usage "No genome found\n" unless ( $genome_file && -s $genome_file );
    pod2usage "Provide at least one set of input files\n"
        unless ( $exonerate_file
            || ($pasa_gff && $pasa_assembly_file && $pasa_peptides
                && $pasa_cds && $pasa_genome_gff)
            || $peptide_file || $mrna_file || $build_only );

    unless ( ($exonerate_file && -s $exonerate_file)
        || ($pasa_gff && -s $pasa_gff
            && $pasa_assembly_file && -s $pasa_assembly_file
            && $pasa_peptides && -s $pasa_peptides
            && $pasa_cds && -s $pasa_cds
            && $pasa_genome_gff && -s $pasa_genome_gff)
        || ($peptide_file && -s $peptide_file)
        || ($mrna_file && -s $mrna_file)
        || $build_only )
    {
        warn "A required input file is missing\n";
        for my $f ( [$exonerate_file, '-exonerate'], [$pasa_gff, '-pasa_gff'],
                    [$pasa_assembly_file, '-pasa_assembly'],
                    [$pasa_peptides, '-pasa_peptides'],
                    [$pasa_cds, '-pasa_cds'],
                    [$pasa_genome_gff, '-pasa_genome'],
                    [$peptide_file, '-peptides'],
                    [$mrna_file, '-mrnas'] )
        {
            warn "\t$f->[1]: $f->[0] not found\n"
                if $f->[0] && !-s $f->[0];
        }
        die "\n";
    }

    if ($pasa_genome_gff && $pasa_genome_gff =~ /\.bed$/) {
        die "Please provide the GFF file for -pasa_genome (got a .bed)\n";
    }

    die "Max intron size (-intron) cannot be 0\n" unless $intron_size && $intron_size > 0;
    die "CPUs (-cpu or -thread) cannot be 0\n"    unless $threads && $threads > 0;
    die "Cannot have both PASA and -peptide\n"   if $pasa_gff && $peptide_file;
    die "Cannot have both PASA and -mrna\n"      if $pasa_gff && $mrna_file;
    die "Cannot have both -mRNA and -peptide\n"  if $mrna_file && $peptide_file;
    die "cDNA mode is only used without peptides\n" if $peptide_file && $is_cdna;

    $is_cdna = 1 if $pasa_gff || $mrna_file;

    die "Softmasked genome file does not exist\n"
        if $softmasked_genome && !-s $softmasked_genome;

    $augustus_dir = readlink($augustus_dir)
        if $augustus_dir && -l $augustus_dir;
    $augustus_dir //= dirname(dirname($augustus_exec)) if $augustus_exec;
    Golden::Augustus::check_augustus() if !$stop_after_correction;

    $same_species = '-same_species' if $same_species;
}
