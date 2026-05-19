package Golden::Alignment;

use strict;
use warnings;
use Carp qw(confess);
use File::Basename;
use Exporter 'import';
use Pod::Usage;
use POSIX qw(ceil);
use List::Util qw(shuffle);
use threads;
use Thread_helper;

# BioPerl + JAMg PerlLib modules used by lifted subs
use Gene_obj;
use Gene_obj_indexer;
use CdbTools;
use GFF3_utils;
use GTF_utils;
use Nuc_translator;
use Fasta_reader;

# Pull cross-module utilities from sibling Golden modules.
# TEMPORARY HACK: sibling modules Golden::Filter and Golden::Augustus do not
# yet export the symbols used here (Filter.pm only exports validate_gene_structure
# + a partial helper set; Augustus.pm does not exist yet). Body calls remain
# fully-qualified (Golden::Filter::xxx, Golden::Augustus::xxx) which compile fine
use Golden::Filter qw(process_cmd filter_gff sort_gff3 get_gff_delimiter
                      order_fasta read_fasta get_id_seq_from_fasta
                      parse_genome_gff gff3_fix_phase remove_overlapping_gff
                      check_program);
use Golden::Augustus qw(wrap_text count_seq split_fasta);

no warnings 'once';                # $main::var globals are populated by the driver

our @EXPORT_OK = qw(
    prepare_pasa_output
    run_aligner
    do_blast_cmd
    do_gmap_cmd
    run_gmap
    run_blast
    run_exonerate
    correct_exonerate_gff
    create_golden_gffs_from_exonerate
    create_golden_gffs_from_gff
    predict_orfs
);

# All globals are owned by the driver (bin/prepare_golden_genes.pl) in package
# `main`. Lifted subs reference them via fully-qualified $main::var syntax.
# This is the simplest mechanism for sharing state without refactoring every
# sub signature. Driver declares them as `our $var;` at script top, and
# GetOptions writes to them.

###############################################################################
# prepare_pasa_output  (v1 lines 2522-2624)
###############################################################################
sub prepare_pasa_output {

# the point of this exercise is to get the cDNA from pasa (not just the ORF) and an
# exonerate compatible annotations file.

 my $fasta_contigs = "$main::pasa_gff.contigs";
 return $fasta_contigs if -s $fasta_contigs;
 print "Processing pasa output to produce exonerate-compatible input\n";

 # will sort and remove overlapping genes. returns valid genes
 my $allowed_data_hashref = Golden::Filter::parse_genome_gff( $main::pasa_genome_gff, 'PASA' );
 print "Indexing PASA assembly file...\n";
 my ($contig_seq_hashref,$contig_length_hashref) = Golden::Filter::read_fasta($main::pasa_assembly_file);
 my $index_file = "$main::pasa_gff.inx";
 my $gene_obj_indexer = new Gene_obj_indexer( { "create" => $index_file } );
 my $asmbl_id_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs( $main::pasa_gff, $gene_obj_indexer );
 die "There was an error indexing the data!" unless $asmbl_id_to_gene_list_href && scalar(keys %{$asmbl_id_to_gene_list_href})>0;

 print "Finding contigs from pasa...\n";
 my ( $contig_counter, $gene_counter, $mrna_counter ) =
   ( int(0), int(0), int(0) );
 my %contigs_used;

 open( ANNOT, ">$fasta_contigs.annotations" );
 open( TROUT, ">$fasta_contigs" );

 foreach my $asmbl_id ( keys %{$asmbl_id_to_gene_list_href} ) {
  my @gene_ids = @{ $asmbl_id_to_gene_list_href->{$asmbl_id} };
  $contig_counter++;
  foreach my $gene_id (@gene_ids) {
   my %params;
   my $gene_obj_ref = $gene_obj_indexer->get_gene($gene_id);
   die "Cannot find gene $gene_id" unless $gene_obj_ref;
   $gene_obj_ref->create_all_sequence_types( \$contig_seq_hashref->{$asmbl_id}, %params );
   $gene_counter++;

   foreach my $isoform ( $gene_obj_ref, $gene_obj_ref->get_additional_isoforms() ) {
    my $com_name   = $isoform->{com_name};
    my $isoform_id = $isoform->{Model_feat_name};
    next unless $allowed_data_hashref->{$isoform_id};

    if ($main::only_complete) {
     next
       if $isoform->{is_5prime_partial}
        || $isoform->{is_3prime_partial}
        || $isoform->{is_pseudogene}
        || ( $com_name && $com_name =~ /prime_partial/ );
    }
    #$mrna_seq is now always given in + orientation
    #my $mrna_seq = $isoform->get_cDNA_sequence();
    # my $orientation = '+';
    my $orientation = $isoform->get_orientation();
    my $mrna_seq    = $contig_seq_hashref->{$asmbl_id};

    # set it to contig as there is no other way to get co-ords currently.

    my ( $model_lend, $model_rend ) =
      sort { $a <=> $b } $isoform->get_model_span();

    # die "Already exists: $isoform_id " if $contigs_used{$isoform_id};
    my $orf_seq_length = int( abs( $model_rend - $model_lend ) + 1 );
    next if $orf_seq_length < $main::minorf;

    if ( $orientation eq '+' ) {
     print ANNOT join(
                       ' ',
                       (
                         $isoform_id, $orientation,
                         $model_lend, $orf_seq_length
                       )
     ) . "\n";
    }
    else {
     print ANNOT join(
      ' ',
      (
       $isoform_id,
       $orientation,

       #the co-ordinates are revcomplemented
       ( length($mrna_seq) - $model_rend + 1 ),
       $orf_seq_length
      )
     ) . "\n";
    }
    print TROUT ">$isoform_id\n" . Golden::Augustus::wrap_text($mrna_seq);
    $contigs_used{$isoform_id} = 1;

   }

  }

 }

 close ANNOT;
 close TROUT;
 die "Could not create $fasta_contigs\n"
   unless -s $fasta_contigs && -s $fasta_contigs . '.annotations';
 print "Prepared $fasta_contigs with "
   . scalar( keys %contigs_used )
   . " sequences larger than $main::minorf b.p. from $mrna_counter transcripts,$gene_counter genes and $contig_counter contigs\n";
 return $fasta_contigs;
}

###############################################################################
# run_aligner  (v1 lines 2804-2854)
###############################################################################
sub run_aligner {
 my $fasta_in = shift;
 my $aligner  = shift;
 my $type     = shift;
 my $seq_count = Golden::Augustus::count_seq($fasta_in);
 my $splits    = ceil( $seq_count / $main::threads );
 print "$aligner: Preparing input $seq_count sequences for up to $main::threads processes\n";
 my $output_directory = $main::cwd.basename($main::genome_file).".vs.".basename($fasta_in).".".$aligner."_dir";

 if ( -d $output_directory ) {
  warn
"$output_directory already exists. Will NOT overwrite and will skip existing output. Stop, delete it and restart otherwise\n";
  sleep(1);
 }
 else {
  if (    ( $main::pasa_cds && $fasta_in eq $main::pasa_cds )
       || ( $main::pasa_peptides && $fasta_in eq $main::pasa_peptides ) )
  {
   Golden::Augustus::split_fasta( $fasta_in, $output_directory, $splits, 'type:complete', 1, 300 );
  }
  else {
   Golden::Augustus::split_fasta( $fasta_in, $output_directory, $splits, undef, 1, 100 );
  }
 }
 my @fastas = glob("$output_directory/*");
 return unless @fastas && scalar(@fastas) > 0;
 print "$aligner: Running ".scalar(@fastas)." alignments\n";

 my $thread_helper = new Thread_helper($main::threads);

 foreach my $fasta ( sort @fastas ) {
  next unless -f $fasta && -s $fasta;
  next unless $fasta =~ /\d$/;
  my $thread = threads->create( 'Golden::Alignment::do_' . $aligner . '_cmd', $fasta,
                                $output_directory,         $type );
  $thread_helper->add_thread($thread);
 }

 $thread_helper->wait_for_all_threads_to_complete();

 my @failed_threads = $thread_helper->get_failed_threads();
 if (@failed_threads) {
  die "Error, " . scalar(@failed_threads) . " threads failed.\n";
  exit(1);
 }
 Golden::Filter::process_cmd("find $output_directory -empty -delete");
 my @aligner_out = glob("$output_directory/*");
 die "$aligner failed" unless @aligner_out && scalar(@aligner_out) > 0;
 print "$aligner: Completed\n";
 return \@aligner_out;
}

###############################################################################
# do_blast_cmd  (v1 lines 2856-2867)
###############################################################################
sub do_blast_cmd {
 my $fasta            = shift;
 my $output_directory = shift;    # currently ignored.
 my $type             = shift;
 my $blast_opt = $type eq 'protein' ? $main::tblastn_exec : $main::tblastx_exec;
 $blast_opt .=
" -num_threads 1 -evalue 1e-20 -max_target_seqs 1 -outfmt 6 -db $main::genome_sequence_file";
 $blast_opt .= " -max_intron_length $main::intron_size" if $type eq 'protein';
 $blast_opt .= " -lcase_masking" if $main::softmasked_genome;
 Golden::Filter::process_cmd("$blast_opt -query $fasta -out $fasta.blast") unless -s "$fasta.blast";
 unlink($fasta);
}

###############################################################################
# do_gmap_cmd  (v1 lines 2869-2881)
###############################################################################
sub do_gmap_cmd {
 my $fasta            = shift;
 my $output_directory = shift; # currently ignored.
 my $type             = shift;
 my $id_fraction = $main::identical_fraction_cutoff / 100;
 my $gmap_opt = "-D $main::gmap_dir -d $main::genome_dbname -f gff3_gene -n 1 --nofails --use-shared-memory=1 --no-chimeras --min-identity=$id_fraction --min-trimmed-coverage=0.80";
 if ($type ne 'protein'){
  Golden::Filter::process_cmd("$main::gmap_exec $gmap_opt --split-output=$fasta.gmap $fasta 2>/dev/null") unless -s "$fasta.gmap.uniq";
 }else{
  die "GMAP only supported for nucleotides not $type\n";
 }
 unlink($fasta);
}

###############################################################################
# run_gmap  (v1 lines 2883-2970)
###############################################################################
sub run_gmap {
 my $fasta       = shift;
 my $gmap_output = shift;
 $gmap_output = $main::cwd.basename($main::genome_file).".vs.".basename($fasta).".gmap.gff3" if !$gmap_output;

 my $output_files_ref = run_aligner( $fasta, 'gmap', 'nuc' );
 my $orig_sep = $/;

 # our gmap output has the gene as the main ID but we need it to
 # the mrna
 print "Post-processing GMAP output...\n";
 unless ( -s $gmap_output ) {
  $/ = "###\n";
  open( OUT, ">$gmap_output" );
  foreach my $file (@$output_files_ref) {
   next unless $file =~ /uniq$/;
   open( IN, $file );
   while ( my $record = <IN> ) {
    chomp($record);
    next if !$record || $record =~ /^\s*$/;
    my @record_data = split( "\n", $record );
    my @corrected_record_data;
    for ( my $i = 0 ; $i < scalar(@record_data) ; $i++ ) {
     next if $record_data[$i] =~ /^#/;
     my @data = split( "\t", $record_data[$i] );
     $data[8] =~ s/\.path\d+//g;
     $data[8] =~ s/\.mrna\d+//;
     $data[1] = 'GMAP';
     push( @corrected_record_data, join( "\t", @data ) );
    }

    my $gene_line = shift @corrected_record_data;
    my $mrna_line = shift @corrected_record_data;

    my @gene_data = split( "\t", $gene_line );
    my @mrna_data = split( "\t", $mrna_line );
    my ( $gene_id, $mrna_id );
    if ( $gene_data[8] =~ /ID=([^;]+)/ ) {
     $gene_id = $1;
    }
    else {
     confess "Unkown GFF format for $record\n$gene_data[8]";
    }
    if ( $mrna_data[8] =~ /ID=([^;]+)/ ) {
     $mrna_id = $1;
    }
    else {
     confess "Unkown GFF format for $record\n$mrna_data[8]";
    }

    # if PASA, convert to pasa format
    if ( $gene_id =~ s/^(asmbl_\d+\|)m(\.\d+)$/$1g$2/ ) {
     $gene_line =~ s/(asmbl_\d+\|)m(\.\d+)/$1g$2/;
    }
    else {
     $gene_id .= '.gene';
     $gene_line =~ s/(ID=[^;]+)/$1.gene/;
    }
    $mrna_line =~ s/Parent=[^;]+/Parent=$gene_id/;
    print OUT $gene_line . "\n";
    print OUT $mrna_line . "\n";

    foreach my $other_line (@corrected_record_data) {
     $other_line =~ s/\.path\d+//g;
     $other_line =~ s/Parent=[^;]+/Parent=$mrna_id/;
     print OUT $other_line . "\n";
    }
    print OUT "##\n";
   }
   close IN;
  }
  close OUT;
  $/ = $orig_sep;
 }

# consider searching for overlaps after we process golden. it will
# much slower ( a lot more genes) but will not discard golden overlaps
# i tried it and found 30 more genes (i.e 0.4%) and costs 10 minutes for 100,000 mrnas

 create_golden_gffs_from_gff( $gmap_output, $gmap_output . '.n' )
   unless ( -s "$gmap_output.passed" );
 Golden::Filter::remove_overlapping_gff( $gmap_output . '.n', $gmap_output . '.nr' )
   unless -s $gmap_output . '.nr';
 $gmap_output .= '.nr';
 create_golden_gffs_from_gff($gmap_output) unless ( -s "$gmap_output.passed" );

 return ( $gmap_output . '.golden', $gmap_output . ".passed" );
}

###############################################################################
# run_blast  (v1 lines 2972-3026)
# Produces a <genome_dir>/<hit>.blast.filter file for each hit.
###############################################################################
sub run_blast {
 my $fasta        = shift;
 my $type         = shift;

 my $blast_output = $main::cwd.basename($main::genome_file).".vs.".basename($fasta).".blast.output";

 print "Preparing for BLAST...\n";
 if ( !-s $blast_output ) {
  my $output_files_ref = run_aligner( $fasta, 'blast', $type );
  foreach my $file (@$output_files_ref) {
   Golden::Filter::process_cmd(
         "sort -nk9,10 $file | sort -s -k1,1 | sort -s -k2,2 >> $blast_output");
  }
 }

 #produce a hit.filter file in $genome_dir
 my %print_hash;

 open( BLAST, $blast_output ) || die("Cannot find $blast_output $!");
 while ( my $ln = <BLAST> ) {
  next if $ln =~ /^\s*$/;
  chomp($ln);
  my @data = split( "\t", $ln );
  next unless $data[10];
  next if $data[10] > 1e-20;

# AAT ( $dstart, $dend,  $score, $astart, $aend,$orient, $zero1, $zero2, $query ) = ( $1, $2, $3, $4, $5, $6, $7, $8, $9 );
# BLAST qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
  my $g_start  = $data[8] < $data[9] ? $data[8] : $data[9];
  my $g_end    = $data[8] < $data[9] ? $data[9] : $data[8];
  my $g_strand = $data[8] < $data[9] ? '1'      : '0';
  my $score =
      $data[11] > 100
    ? $data[11]
    : 101
    ; # force to keep all alignments, even the very short ones if evalue is high enough
  $print_hash{ $data[1] } .=
      "$g_start\t$g_end\t"
    . $score . "\t"
    . $data[6] . "\t"
    . $data[7]
    . "\t$g_strand\t0\t0\t"
    . $data[0] . "\n";

 }
 foreach my $hit ( keys %print_hash ) {
  open( OUT, ">$main::genome_dir/$hit.blast.filter" )
    || die("Cannot write to $main::genome_dir/$hit.blast.filter $!");
  print OUT "\t$hit\tProduced by BLAST\n";
  print OUT $print_hash{$hit};
  close OUT;
 }
 close BLAST;
 return $blast_output;
}

###############################################################################
# run_exonerate  (v1 lines 2626-2781)
# DEVIATION FROM VERBATIM LIFT: the PASA branch (line ~2719 "we have pasa data")
# previously called &run_aat($fasta_contigs,'nucl') guarded by `if (!$same_species)`.
# Per v2 plan §4: AAT.pl + filter binary are not in jamg.sif; BLAST is. So both
# the guard and the run_aat call are replaced with an unconditional
# Golden::Alignment::run_blast($fasta_contigs, 'nucl').
###############################################################################
sub run_exonerate {
 my $is_blast;
 ($main::parafly_exec) = Golden::Filter::check_program('ParaFly');

 if ($main::pasa_gff){
	$main::exonerate_file = basename($main::pasa_gff) . '.exonerate.results';
 }elsif($main::peptide_file){
	$main::exonerate_file = basename($main::peptide_file) . '.exonerate.results';
 }elsif($main::mrna_file){
	$main::exonerate_file = basename($main::mrna_file) . '.exonerate.results';
 }
 my $exonerate_command_file = "run_exonerate_commands.cmd";

 # check if it already has been processed
 $main::no_rerun_exonerate = 1
   if ( -s $main::exonerate_file
  && -s $exonerate_command_file
  && ( -s $exonerate_command_file == -s $exonerate_command_file . '.completed' )
   );
 unless ($main::no_rerun_exonerate) {
  print "Finding rough co-ordinates...\n";
  if (
   -s $exonerate_command_file
   && ( !-s $exonerate_command_file . '.completed'
    || (
     -s $exonerate_command_file != -s $exonerate_command_file . '.completed' ) )
    )
  {
   print "Re-processing with exonerate\n";
   Golden::Filter::process_cmd("$main::parafly_exec -shuffle -CPU $main::threads -c $exonerate_command_file -failed_cmds $exonerate_command_file.failed -v "   );
  }
  else {

   #TODO mrna_file only?
   if ($main::peptide_file) {

    # v2: AAT removed from jamg.sif; the peptide path now always seeds
    # exonerate from BLAST (was: AAT when !same_species, BLAST otherwise).
    $is_blast = 1;
    my $blast_output = run_blast( $main::peptide_file, 'protein' );
    print "Finding accurate co-ordinates using exonerate...\n";
    unlink($main::exonerate_file);
    my @filter_files = glob("$main::genome_dir/*filter");
    die "No filter files found" unless @filter_files;
    my $exonerate_options = " -minorf $main::minorf -protein -in $main::peptide_file -separate -aat_dir $main::genome_dir -aat_suffix filter -threads $main::threads "
	."-intron_max $main::intron_size  $main::same_species ";
    $exonerate_options .= " -softmask -ref $main::softmasked_genome "
      if ($main::softmasked_genome);
    $exonerate_options .= " -ref $main::genome_file "
      if ( !$main::softmasked_genome );
    $exonerate_options .= " -norefine "   if $main::norefine;
    $exonerate_options .= " -local_protein "   if !$main::do_exhaustive;
    $exonerate_options .= " -from_blast " if $is_blast;
    my $aat_score = $main::same_species ? 100 : 20;
    $exonerate_options .= " -score_dps $aat_score";

    Golden::Filter::process_cmd( 'run_exonerate.pl' . $exonerate_options );
    die "Exonerate run failed for some reason....\n" if ( !-d basename($main::peptide_file) . "_queries" );
    my @files_to_cat = glob(basename($main::peptide_file)."_queries/*exonerate_results");
    die "No *exonerate_results files found!" unless @files_to_cat && $files_to_cat[0];
    foreach my $f (@files_to_cat){
	 system("cat $f >> $main::exonerate_file");
    }
   }
   elsif ($main::mrna_file) {
    #die "Not implemented yet - see GMAP output\n";
     my $gmap_output = $main::cwd.basename($main::genome_file).".vs.".basename($main::mrna_file).".gmap.gff3";
     die "Cannot find GMAP output $gmap_output\n" unless -s $gmap_output;

    print "Finding accurate co-ordinates using exonerate...\n";
    unlink($main::exonerate_file);
    my $exonerate_options = " -minorf $main::minorf -in $main::mrna_file -separate -threads $main::threads "
	."-intron_max $main::intron_size  $main::same_species -gmap_file $gmap_output ";
    $exonerate_options .= " -softmask -ref $main::softmasked_genome "
      if ($main::softmasked_genome);
    $exonerate_options .= " -ref $main::genome_file "
      if ( !$main::softmasked_genome );
    $exonerate_options .= " -norefine "   if $main::norefine;
    $exonerate_options .= " -from_blast " if $is_blast;
    my $aat_score = $main::same_species ? 100 : 20;
    $exonerate_options .= " -score_dps $aat_score";
    Golden::Filter::process_cmd( 'run_exonerate.pl' . $exonerate_options );
    die "Exonerate run failed for some reason....\n" if ( !-d basename($main::mrna_file) . "_queries" );
    my @files_to_cat = glob(basename($main::mrna_file)."_queries/*exonerate_results");
    die "No *exonerate_results files found!" unless @files_to_cat && $files_to_cat[0];
    foreach my $f (@files_to_cat){
	 system("cat $f >> $main::exonerate_file");
    }

   }
   else {

    # we have pasa data. no peptide file
    my $fasta_contigs = prepare_pasa_output();
    # v2 deviation: was `&run_aat($fasta_contigs, 'nucl') if !$same_species`.
    # AAT is removed in v2; BLAST replaces it. Guard dropped, call unconditional.
    Golden::Alignment::run_blast( $fasta_contigs, 'nucl' );
    $is_blast = 1;
    unlink($main::exonerate_file);
    print "Finding accurate co-ordinates using exonerate...\n";
    my $exonerate_options =" -minorf $main::minorf -annotation $fasta_contigs.annotations -in $fasta_contigs -separate -aat_dir $main::genome_dir -aat_suffix .filter "
      . " -threads $main::threads -intron_max $main::intron_size $main::same_species ";
    $exonerate_options .= " -softmask -ref $main::softmasked_genome "
      if ($main::softmasked_genome);
    $exonerate_options .= " -ref $main::genome_file "
      if ( !$main::softmasked_genome );
    $exonerate_options .= " -norefine "   if $main::norefine;
    $exonerate_options .= " -local_protein "   if !$main::do_exhaustive;
    $exonerate_options .= " -from_blast " if $is_blast;
    my $aat_score = $main::same_species ? 100 : 20;
    $exonerate_options .= " -score_dps $aat_score";

    Golden::Filter::process_cmd( 'run_exonerate.pl' . $exonerate_options );
    die "Exonerate run failed for some reason....\n"
      if ( !-d basename($fasta_contigs) . "_queries" );
    my @files_to_cat = glob(basename($fasta_contigs)."_queries/*exonerate_results");
    die "No *exonerate_results files found!" unless @files_to_cat && $files_to_cat[0];
    foreach my $f (@files_to_cat){
	 system("cat $f >> $main::exonerate_file");
    }
   }

  }
  $main::no_rerun_exonerate = 1
    if ( -s $main::exonerate_file
   && -s $exonerate_command_file
   && -s $exonerate_command_file . '.completed'
   && (
    -s $exonerate_command_file == -s $exonerate_command_file . '.completed' ) );

 }
  $main::no_rerun_exonerate = 1
    if ( -s $main::exonerate_file
   && -s $exonerate_command_file
   && -s $exonerate_command_file . '.completed'
   && (
    -s $exonerate_command_file != -s $exonerate_command_file . '.completed' ) );
 die"Have not completed with the current exonerate run. Rerun this program once more. If it still fails, you can force this message to be ignored with -norerun \n"
   unless $main::no_rerun_exonerate;
 print "Completed processing exonerate file as $main::exonerate_file\n";

 correct_exonerate_gff($main::exonerate_file)
   unless -s $main::exonerate_file . '.corrected.gff3'
    && $main::exonerate_file . '.corrected.passed';
 die "Cannot find $main::exonerate_file.corrected.gff3\n"
   unless -s $main::exonerate_file . '.corrected.gff3';
 if ($main::stop_after_correction) {
  print "User asked to stop after exonerate correction\n";
  exit();
 }

 return ( $main::exonerate_file . '.corrected.golden',
          $main::exonerate_file . '.corrected.passed' );
}

###############################################################################
# correct_exonerate_gff  (v1 lines 398-740)
###############################################################################
sub correct_exonerate_gff {
 my $exonerate_file = shift;
 print "Processing $exonerate_file using these cut-offs: identical:$main::identical_fraction_cutoff similar:$main::similar_fraction_cutoff;  No more than $main::mismatch_cutoff mismatches in total; Methionine and * in protein sequence; Start codon and stop codon in genome; No overlapping genes; Splice sites (GT..AG or GC..AG)\n" unless $main::liberal_cutoffs;
 print "Processing $exonerate_file using few cut-offs (liberal): identical:$main::identical_fraction_cutoff similar:$main::similar_fraction_cutoff;  No more than $main::mismatch_cutoff mismatches in total\n" if $main::liberal_cutoffs;

 #pod2usage unless $exonerate_file && -s $exonerate_file;
 # IF WE NEED TO GET WHOLE SCAFFOLD SEQUENCE:
 pod2usage "Cannot find a file $exonerate_file or $main::genome_sequence_file\n"
   unless $exonerate_file
    && $main::genome_sequence_file
    && -s $exonerate_file
    && -s $main::genome_sequence_file;
 open( EXONERATE,               $exonerate_file );
 open( CORRECTED_EXONERATE_GFF, ">$exonerate_file.corrected.gff3" );

 my ( %details, $flag, %vulgar_data, %already_printed );
 my ( $gene_counter, $scounter, $ecounter ) = ( int(0), int(0), int(0) );

###########################################################################################
# EXONERATE protein has no UTR, the gene is actually the coding part of the mRNA
# get gene data
###########################################################################################
 my $header = <EXONERATE>;
 die
"It seems that the exonerate was run using protein2genome but -cdna was requested. Aborting\n"
   if ( $header =~ /protein2genome/ && $main::is_cdna );
 die
"It seems that the exonerate was run using cdna2genome but -cdna was not requested. Aborting\n"
   if ( ( $header =~ /cdna2genome/ || $header =~ /coding2genome/ )
        && !$main::is_cdna );
 while ( my $ln = <EXONERATE> ) {
  next if $ln =~ /^\s*$/;
  if ( $ln =~ /^vulgar:\s+(\S+)/ ) {
   $vulgar_data{$1} = $ln;
  }
  if ( $ln =~ /^# --- START OF GFF DUMP/ ) {
   my @exonerate_data;
   my @overlap_check;
   while ( my $ln2 = <EXONERATE> ) {
    last if ( $ln2 =~ /^# --- END OF GFF DUMP/ );
    next if ( $ln2 =~ /^#/ );
    my @data = split( "\t", $ln2 );
    if ( $data[8] ) {
     next if $data[2] eq 'similarity';
     my $offset = int(0);

# this is where the offset is taken from. hope that people don't run with references with this coords...
     if ( $data[0] =~ s/:(\d+)-(\d+)$// ) {
      $offset = $1 - 1;
      $data[3] += $offset;
      $data[4] += $offset;
     }
     push( @exonerate_data, join( "\t", @data ) );
     if ( $data[2] eq 'cds' ) {
      my %h = ( 'start' => $data[3], 'end' => $data[4] );
      push( @overlap_check, \%h );
     }
    }
   }
###########################################################################################
   # making other corrections. and fixing exonerate, snap etc bugs
###########################################################################################
   my ( $gene_start, $gene_end, $mRNA, $mRNA_id );
   my $exon_counter   = int(0);
   my $cds_counter    = int(0);
   my $splice_counter = int(0);
   $exonerate_data[0] =~ /sequence (\S+)/;
   my $gene_id = $1 if $1;

# if it has already been printed, skip it. someitmes there is a bug and this creates non-unique IDs
   if ( $already_printed{$gene_id} ) {
    warn "Warning: $gene_id was already printed before. skipping...\n";
    next;
   }
   $already_printed{$gene_id} = 1;

   for ( my $i = 0 ; $i < scalar(@exonerate_data) ; $i++ ) {
    my $gff_line = $exonerate_data[$i] || next;
    my @data = split( "\t", $gff_line );
    if ( $data[2] eq 'gene' ) {
     $data[8] =~ /sequence\s+(\S+)/;
     $gene_id    = $1;
     $gene_start = $data[3];
     $gene_end   = $data[4];
     undef($mRNA);    #reset
     undef($mRNA_id);

     #PASA HACK - review if there are problems in the future
     if ( $gene_id =~ s/^m\.(\d+)$/g.$1/ ) {
      $mRNA_id = "m.$1";
      $data[8] = "ID=$gene_id\n";
      $mRNA =
          $data[0] . "\t"
        . $data[1]
        . "\tmRNA\t"
        . $data[3] . "\t"
        . $data[4] . "\t.\t"
        . $data[6] . "\t.\t"
        . "ID=$mRNA_id;Parent=$gene_id\n";
     }

     # not PASA
     else {
      $mRNA_id = $gene_id . ".mRNA";
      $mRNA =
          $data[0] . "\t"
        . $data[1]
        . "\tmRNA\t"
        . $data[3] . "\t"
        . $data[4] . "\t.\t"
        . $data[6] . "\t.\t"
        . "ID=$mRNA_id;Parent=$gene_id\n";

      $data[8] = "ID=$gene_id\n";
     }
     $details{$gene_id}{'source'} = $data[1];
    }
    elsif ( $data[2] =~ /utr/ ) {
     if ( $data[2] =~ /5/ ) {
      $data[2] = 'five_prime_utr';
     }
     elsif ( $data[2] =~ /3/ ) {
      $data[2] = 'three_prime_utr';
     }
     else {
      $data[2] = 'UTR';
     }

     $data[8] = "ID=$gene_id." . $data[2] . ";Parent=$mRNA_id\n";
    }
    elsif ( $data[2] =~ /splice/ ) {
     $splice_counter++;
     if ( $data[8] =~ /splice_site "([A-Z]+)"/ ) {
      my $site = $1;
      if ( $data[2] eq 'splice5' ) {
       unless ( $site eq 'GT' || $site eq 'GC' ) {
        $details{$gene_id}{'non_canonical_splice_sites'}++;
       }
      }
      elsif ( $data[2] eq 'splice3' && $site ne 'AG' ) {
       $details{$gene_id}{'non_canonical_splice_sites'}++;
      }
     }
     $data[8] =
       "ID=$gene_id." . $data[2] . ":$splice_counter;Parent=$mRNA_id\n";
     $data[2] = 'splice_junction';

    }
    elsif ( $data[2] eq 'cds' ) {
     if ( $cds_counter > 0 ) {
      ### check if the next line is an exon. if it is not, check if there is an overlapping previous CDS (exonerate cdna2genome bug)
      # if there is we have to delete $exonerate_data[$i]; and move on
      my $previous = $overlap_check[ $cds_counter - 1 ];
      my @next_data = split( "\t", $exonerate_data[ $i + 1 ] )
        if $exonerate_data[ $i + 1 ];
      if (    $previous
           && $next_data[2]
           && $next_data[2] ne 'exon' )
      {
       if ( $data[6] eq '+' ) {

        # start must be higher than previous end
        if ( $data[3] < $previous->{'end'} ) {

#warn $data[8].":Start is smaller than previous end:\n".$data[3]." vs ".$previous->{'end'}."\n";
         delete( $exonerate_data[$i] );
         next;
        }
       }
       else {

        # end must be lower than previous start
        if ( $data[4] > $previous->{'start'} ) {
         delete( $exonerate_data[$i] );
         next;
        }
       }
      }
     }
     elsif ( scalar(@overlap_check) == 1 ) {

# if this is the first and only CDS, then it is a single exon gene. exonerate cdna2genome bug includes UTR co-ordidates. fix it
# by getting previous UTR line and adjust co-ordidates to by utr[end]+1 if + or utr[start]-1 if -.
      my @previous_data =
        split( "\t", $exonerate_data[ $i - 1 ] );
      if ( @previous_data && $previous_data[2] =~ /utr/i ) {
       if ( $data[6] eq '+' ) {
        if ( $data[3] < $previous_data[4] ) {
         $data[3] = $previous_data[4] + 1;
        }
       }
       else {
        if ( $data[4] > $previous_data[3] ) {
         $data[4] = $previous_data[3] - 1;
        }
       }
      }
     }
     $cds_counter++;
     $data[2] = 'CDS';
     $data[8] = "ID=cds.$gene_id;Parent=$mRNA_id\n";
    }
    elsif ( $data[2] eq 'exon' ) {
     $exon_counter++;
     $data[8] = "ID=$gene_id.exon.$exon_counter;Parent=$mRNA_id\n";
    }
    elsif ( $data[2] eq 'intron' ) {

     # we capture this from the GFF now;
     next;

     #$intron_counter++;
     #$data[8] = "ID=$gene_id.intron.$intron_counter;Parent=$mRNA_id\n";
    }
    $data[8] =~ s/\s+\;\s+/\;/g;
    $exonerate_data[$i] = join( "\t", @data );
   }
#################################################################################
   # now printing data
################################################################################
   for ( my $i = 0 ; $i < scalar(@exonerate_data) ; $i++ ) {
    my $gff_line = $exonerate_data[$i] || next;
    my @data = split( "\t", $gff_line );
    if ( $data[2] eq 'gene' ) {
     print CORRECTED_EXONERATE_GFF join( "\t", @data );
     print CORRECTED_EXONERATE_GFF $mRNA if $mRNA;
     $gene_counter++;
    }
    else {
     print CORRECTED_EXONERATE_GFF join( "\t", @data );
    }
   }

   # delimit end of gene model
   print CORRECTED_EXONERATE_GFF "##\n";

  }
###########################################################################################
  # for golden set:
  if ( $ln =~ /^RYOAP_END/ ) {
   $flag = 0;
   $ecounter++;
  }
  if ( $ln =~ /^RYOAP_START/ ) {
   $flag = 1;
   $scounter++;
   my $stats_header = <EXONERATE>;
   my $stats_str    = <EXONERATE>;
   chomp($stats_str);
   my @stats_data = split( "\t", $stats_str );
   my $query_data_str = <EXONERATE>;
   chomp($query_data_str);
   my @query_data = split( "\t", $query_data_str );
   my $query_id = $query_data[1];
   my $query_start =
       $query_data[2] <= $query_data[3]
     ? $query_data[2]
     : $query_data[3];
   my $query_end =
       $query_data[3] >= $query_data[2]
     ? $query_data[3]
     : $query_data[2];
   my $query_fasta = <EXONERATE>;
   my $query_seq;
   $ln = <EXONERATE>;

   while ( $ln !~ /^RYOAP/ ) {
    chomp($ln);
    $query_seq .= $ln;
    $ln = <EXONERATE>;
   }

   #       $query_fasta .= $query_seq;
   #TARGET
   my $target_data_str = $ln;
   chomp($target_data_str);
   my @target_data = split( "\t", $target_data_str );
   my $target_id = $target_data[1];
   my $target_start =
       $target_data[2] <= $target_data[3]
     ? $target_data[2]
     : $target_data[3];
   my $target_end =
       $target_data[3] >= $target_data[2]
     ? $target_data[3]
     : $target_data[2];
   my $alignment_strand = $target_data[2] <= $target_data[3] ? int(1) : int(-1);
   my $target_fasta = <EXONERATE>;
   my $target_seq;
   $ln = <EXONERATE>;

   while ( $ln !~ /^RYOAP/ ) {
    chomp($ln);
    $target_seq .= $ln;
    $ln = <EXONERATE>;
   }
   if ( $ln =~ /^RYOAP_END/ ) {
    $flag = 0;
    $ecounter++;
   }

   #       $target_fasta .= $target_seq;
   my $qlength = $query_data[3] - $query_data[2];
   my $tlength = abs( $target_data[3] - $target_data[2] );
   my $offset  = int(0);
   if ( $target_id =~ s/:(\d+)\-\d+$// ) {
    $offset = $1 - 1;
   }
   $details{$query_id}{'counter'}                    = $scounter;
   $details{$query_id}{'qlength'}                    = $qlength;
   $details{$query_id}{'tlength'}                    = $tlength + $offset;
   $details{$query_id}{'score'}                      = $stats_data[1];
   $details{$query_id}{'mismatch_number'}            = $stats_data[5];
   $details{$query_id}{'identical_frac'}             = $stats_data[6];
   $details{$query_id}{'similar_frac'}               = $stats_data[7];
   $details{$query_id}{'query_id'}                   = $query_id;
   $details{$query_id}{'query_seq'}                  = uc($query_seq);
   $details{$query_id}{'target_seq'}                 = uc($target_seq);
   $details{$query_id}{'query_start'}                = $query_start;
   $details{$query_id}{'query_end'}                  = $query_end;
   $details{$query_id}{'target_id'}                  = $target_id;
   $details{$query_id}{'target_start'}               = $target_start + $offset;
   $details{$query_id}{'target_end'}                 = $target_end + $offset;
   $details{$query_id}{'alignment_strand'}           = $alignment_strand;
   $details{$query_id}{'non_canonical_splice_sites'} = int(0)
     if !$details{$query_id}{'non_canonical_splice_sites'};
  }
 }
 print "Processed $gene_counter gene models.\n";
 close(EXONERATE);
 close CORRECTED_EXONERATE_GFF;
 die "Report seems incomplete!\n"
   unless $scounter == $ecounter;
 die "No genes have been found!\n" if $scounter == 0;

 if ($main::stop_after_correction) {
  print "User asked to stop after correction\n";
  exit();
 }
 create_golden_gffs_from_exonerate( $exonerate_file . '.corrected.gff3',
                                     \%details, \%vulgar_data )
   unless -s "$exonerate_file.corrected.passed";
}

###############################################################################
# create_golden_gffs_from_exonerate  (v1 lines 742-1032)
###############################################################################
sub create_golden_gffs_from_exonerate {
 my $exonerate_file      = shift;
 my $details_hashref     = shift;
 my $vulgar_data_hashref = shift;

 Golden::Filter::gff3_fix_phase($exonerate_file) unless -s $exonerate_file . '.gtf';

 # it's getting a bit long...
 my $basename_exonerate = $exonerate_file;
 $basename_exonerate =~ s/\.gff3$//;
 print "Finding golden subset...\n";

 my $query_outfile =
   $main::is_cdna
   ? "$basename_exonerate.passed.cDNA"
   : "$basename_exonerate.passed.protein";
 open( OUT, ">$basename_exonerate.data" ) unless $main::nodataprint;
 open( VULGAR, ">$basename_exonerate.passed.vulgar" );
 open( LOG,    ">$basename_exonerate.golden.log" );
 open( PASSP,  ">$query_outfile" );

 #  open( PASSG,  ">$basename_exonerate.passed.genome" )  unless $main::nodataprint;
 open( PASSGC, ">$basename_exonerate.passed.genome.cds" ) unless $main::nodataprint;
 print OUT
"#Hit number\tQuery length\tTarget length\tScore\tNumber of mismatches\tFraction of identical\tFraction of similar\tNon-canonical splice_sites\tquery name\tquery start\tquery end\tquery sequence\tAlignment direction\ttarget name\ttarget start\ttarget end\ttarget sequence\n"
   unless $main::nodataprint;
 my $query_number = keys %{$details_hashref};
 my $pass_counter = int(0);
 my %pass_check;
 my %accepted;
 my $region_exists = int(0);
 my $counter       = int(0);

 foreach my $query_id (
  sort {
   $details_hashref->{$b}{'score'} <=> $details_hashref->{$a}{'score'}
  }
  keys %{$details_hashref}
   )
 {
  $counter++;
  print "\r$counter/$query_number    $pass_counter passed           ";
  print OUT "Hit "
    . $details_hashref->{$query_id}{'counter'} . "\t"
    . $details_hashref->{$query_id}{'qlength'} . "\t"
    . $details_hashref->{$query_id}{'tlength'} . "\t"
    . $details_hashref->{$query_id}{'score'} . "\t"
    . $details_hashref->{$query_id}{'mismatch_number'} . "\t"
    . $details_hashref->{$query_id}{'identical_frac'} . "\t"
    . $details_hashref->{$query_id}{'similar_frac'} . "\t"
    . $details_hashref->{$query_id}{'non_canonical_splice_sites'} . "\t"
    . $query_id . "\t"
    . $details_hashref->{$query_id}{'query_start'} . "\t"
    . $details_hashref->{$query_id}{'query_end'} . "\t"
    . $details_hashref->{$query_id}{'query_seq'} . "\t"
    . $details_hashref->{$query_id}{'alignment_strand'} . "\t"
    . $details_hashref->{$query_id}{'target_id'} . "\t"
    . $details_hashref->{$query_id}{'target_start'} . "\t"
    . $details_hashref->{$query_id}{'target_end'} . "\t"
    . $details_hashref->{$query_id}{'target_seq'} . "\n"
    unless $main::nodataprint;
  if ( !$main::liberal_cutoffs &&
       $pass_check{
            's'
          . $details_hashref->{$query_id}{'target_start'} . 'e'
          . $details_hashref->{$query_id}{'target_end'}
       }
    )
  {
   print LOG
     "$query_id rejected: higher scoring alignment with same co-ordinates ("
     . 's'
     . $details_hashref->{$query_id}{'target_start'} . 'e'
     . $details_hashref->{$query_id}{'target_end'} . "):"
     . $pass_check{ 's'
      . $details_hashref->{$query_id}{'target_start'} . 'e'
      . $details_hashref->{$query_id}{'target_end'} }
     . "\n";
   $region_exists++;
   next;
  }
  unless (
        $details_hashref->{$query_id}{'identical_frac'} >= $main::identical_fraction_cutoff )
  {
   print LOG
     "$query_id rejected: identical_frac less than $main::identical_fraction_cutoff: "
     . $details_hashref->{$query_id}{'identical_frac'} . "\n";
   $main::failed_cutoff++;
   next;
  }
  unless ( $details_hashref->{$query_id}{'similar_frac'} >= $main::similar_fraction_cutoff )
  {
   print LOG "$query_id rejected: similar_frac less than $main::similar_fraction_cutoff: "
     . $details_hashref->{$query_id}{'similar_frac'} . "\n";
   $main::failed_cutoff++;
   next;
  }
  unless (
          $details_hashref->{$query_id}{'mismatch_number'} <= $main::mismatch_cutoff )
  {
   print LOG "$query_id rejected: mismatches more than $main::mismatch_cutoff: "
     . $details_hashref->{$query_id}{'mismatch_number'} . "\n";
   $main::failed_cutoff++;
   next;
  }
  if ( !$main::liberal_cutoffs && $details_hashref->{$query_id}{'non_canonical_splice_sites'} > 0 ) {
   print LOG "$query_id rejected: there are non-canonical splice sites : "
     . $details_hashref->{$query_id}{'non_canonical_splice_sites'} . "\n";
   $main::failed_cutoff++;
   next;
  }

  unless ( $main::liberal_cutoffs || !$main::only_complete || $main::is_cdna
          || substr( $details_hashref->{$query_id}{'query_seq'}, 0, 1 ) eq 'M' )
  {
   print LOG "$query_id rejected: query does not start with M: "
     . substr( $details_hashref->{$query_id}{'query_seq'}, 0, 1 ) . "\n";
   $main::failed_cutoff++;
   next;
  }
  unless ( $main::liberal_cutoffs || !$main::only_complete || $main::is_cdna
         || substr( $details_hashref->{$query_id}{'query_seq'}, -1, 1 ) eq '*' )
  {
   print LOG "$query_id rejected: query does not end with stop codon (*): "
     . substr( $details_hashref->{$query_id}{'query_seq'}, -1, 1 ) . "\n";
   $main::failed_cutoff++;
   next;
  }
  unless ( $main::liberal_cutoffs || !$main::only_complete ||
          substr( $details_hashref->{$query_id}{'target_seq'}, 0, 3 ) eq 'ATG' )
  {
   print LOG "$query_id rejected: target does not start with ATG: "
     . substr( $details_hashref->{$query_id}{'target_seq'}, 0, 3 ) . "\n";
   $main::failed_cutoff++;
   next;
  }
  unless ( $main::liberal_cutoffs || !$main::only_complete ||
       (
           substr( $details_hashref->{$query_id}{'target_seq'}, -3, 3 ) eq 'TAG'
        || substr( $details_hashref->{$query_id}{'target_seq'}, -3, 3 ) eq 'TAA'
        || substr( $details_hashref->{$query_id}{'target_seq'}, -3, 3 ) eq 'TGA'
       )
    )
  {
   print LOG
"$query_id rejected: target does not end with stop codon (TAG, TAA or TGA): "
     . substr( $details_hashref->{$query_id}{'target_seq'}, -3, 3 ) . "\n";
   $main::failed_cutoff++;
   next;
  }
  my $scaffold_seq =
    $main::scaffold_seq_hashref->{ $details_hashref->{$query_id}{'target_id'} };
  if ( !$scaffold_seq ) {
   warn "Cannot find genome entry "
     . $details_hashref->{$query_id}{'target_id'} . "\n";
   next;
  }

  #bug check
  unless (int( $details_hashref->{$query_id}{'target_start'} )
       && int( $details_hashref->{$query_id}{'target_end'} )
       && $details_hashref->{$query_id}{'target_start'} > 0
       && $details_hashref->{$query_id}{'target_end'} > 0
       && $details_hashref->{$query_id}{'target_end'} <= length($scaffold_seq) )
  {
   print LOG "$query_id had to be skipped: has issues with target end ("
     . $details_hashref->{$query_id}{'target_id'}
     . ") being negative or larger than scaffold size - possible corruption of exonerate file! ("
     . $details_hashref->{$query_id}{'target_end'}
     . " is <0 or > "
     . length($scaffold_seq) . ")\n";
   next;
  }
  my $scaffold_subseq;

  #consider: adding stop codon as it is not reported at end??
  if ( $details_hashref->{$query_id}{'alignment_strand'} == -1 ) {
   $scaffold_subseq = Golden::Filter::_revcomp(
    substr(
            $scaffold_seq,
            $details_hashref->{$query_id}{'target_start'} - 2 - 1,
            $details_hashref->{$query_id}{'target_end'} - 1
      )    # starts from 0
   );
  }
  else {
   $scaffold_subseq = substr(
                              $scaffold_seq,
                              $details_hashref->{$query_id}{'target_start'},
                              $details_hashref->{$query_id}{'target_end'} + 2
   );
  }
  if ( $main::no_single_exon
       && length($scaffold_subseq) <=
       length( $details_hashref->{$query_id}{'target_seq'} ) )
  {
   $main::failed_cutoff++;
   print LOG "$query_id rejected: user requested no single coding exon\n";
   next;
  }

  my $number_of_ns = $scaffold_subseq =~ tr/N/N/;
  if ( $number_of_ns > ( length($scaffold_subseq) * 0.40 ) ) {
   $main::failed_cutoff++;
   print LOG
"$query_id rejected: genome sequence is more than 40% Ns ($number_of_ns)\n";
   next;
  }
  print LOG "$query_id accepted\n";
  print VULGAR $vulgar_data_hashref->{$query_id}
    if $vulgar_data_hashref->{$query_id};
  $pass_check{ 's'
     . $details_hashref->{$query_id}{'target_start'} . 'e'
     . $details_hashref->{$query_id}{'target_end'} } = $query_id;
  print PASSP ">$query_id\n"
    . $details_hashref->{$query_id}{'query_seq'} . "\n";

#    print PASSG ">$query_id "     . $details_hashref->{$query_id}{'target_id'} . " genome sequence\n$scaffold_subseq\n" unless $main::nodataprint;
  print PASSGC ">$query_id "
    . $details_hashref->{$query_id}{'target_id'}
    . " genome CDS\n"
    . $details_hashref->{$query_id}{'target_seq'} . "\n"
    unless $main::nodataprint;
  $accepted{$query_id} = $details_hashref->{$query_id}{'source'};
  $pass_counter++;
 }
 print "\r$counter/$query_number    $pass_counter passed           \n";
 close OUT   unless $main::nodataprint;
 close PASSP unless $main::nodataprint;

 #  close PASSG unless $main::nodataprint;
 close VULGAR;

 my $number_of_passing_genes = scalar( keys %accepted );
 open( PASS, ">$basename_exonerate.passed" );
 my @shuffled_genes = shuffle( keys %accepted );
 foreach my $gene (@shuffled_genes) {
  print PASS $gene . "\t" . $accepted{$gene} . "\n";
 }
 close PASS;

 my $orig_sep = $/;
 $/ = Golden::Filter::get_gff_delimiter($exonerate_file);
 open( IN,  $exonerate_file );
 open( OUT, ">$basename_exonerate.golden" );

GENE: while ( my $record = <IN> ) {
  chomp($record);
  next if !$record || $record =~ /^\s*$/;
  my @record_data = split( "\n", $record );
  my $gene_line   = shift @record_data;
  my $mrna_line   = shift @record_data;
  my @gene_data   = split( "\t", $gene_line );
  my @mrna_data   = split( "\t", $mrna_line );
  my ( $gene_id, $mrna_id );

  if ( $gene_data[8] =~ /ID=([^;]+)/ ) {
   $gene_id = $1;
  }
  else {
   confess "Unkown GFF format for $record\n$gene_data[8]";
  }
  if ( $mrna_data[8] =~ /ID=([^;]+)/ ) {
   $mrna_id = $1;
  }
  else {
   confess "Unkown GFF format for $record\n$mrna_data[8]";
  }
  print OUT $record . $/ if $accepted{$gene_id} || $accepted{$mrna_id};
 }
 close OUT;
 $/ = $orig_sep;

 print LOG
"Processed $counter transcripts.\nFound $number_of_passing_genes sequences passing criteria.\n\n";
 close LOG;
 print
"Processed $counter transcripts.\nFound $number_of_passing_genes sequences passing criteria.\n";
 print "See $basename_exonerate.golden.log for details what happened to each gene\n\n";

 #cleanup
 unlink( $exonerate_file . '.pep' );
 unlink( $exonerate_file . '.cds' );
 unlink( $exonerate_file . '.gene' );
 unlink( $exonerate_file . '.original' );
}

###############################################################################
# create_golden_gffs_from_gff  (v1 lines 3106-3260)
###############################################################################
sub create_golden_gffs_from_gff {
 my $gff_file = shift;
 my $output   = shift;
 $output = $gff_file . '.golden' if !$output;
 my $delimiter = shift;
 my $orig_sep  = $/;

 Golden::Filter::gff3_fix_phase($gff_file);

 $delimiter = Golden::Filter::get_gff_delimiter($gff_file) if !$delimiter;
 print "Finding golden subset...\n";

 #TODO :
 # if thomas changes the gmap format:
 # $mismatch_cutoff
 #DONE:
 # $identical_fraction_cutoff
 # mrna: coverage=100.0;identity=100.0
 # pep: M and *
 # genome my $number_of_ns = $scaffold_subseq =~ tr/N/N/; (40%)
 # splice sites
 # acceptor_site=AG
 # donor_site=GC
 # donor_site=GT

 my ( %tocheck, %accepted );
 my ( $counter, $failed_cutoff ) = ( int(0), int(0), int(0) );
 my ($gene_seqs_hashref,$gene_length_hashref) = Golden::Filter::read_fasta( $gff_file . '.gene' );
 my ($pep_seqs_hashref,$pep_length_hashref)  = Golden::Filter::read_fasta( $gff_file . '.pep' );

 $/ = $delimiter;

 open( LOG, ">$output.log" );
 open( OUT, ">$output" );
 open( GFF, $gff_file );

GENE: while ( my $record = <GFF> ) {
  chomp($record);
  next if !$record || $record =~ /^\s*$/;
  my @record_data = split( "\n", $record );
  my $gene_line   = shift @record_data;
  my $mrna_line   = shift @record_data;
  my @gene_data   = split( "\t", $gene_line );
  my @mrna_data   = split( "\t", $mrna_line );
  my ( $gene_id, $mrna_id );

  if ( $gene_data[8] =~ /ID=([^;]+)/ ) {
   $gene_id = $1;
  }
  else {
   confess "Unkown GFF format for $record\n$gene_data[8]";
  }
  if ( $mrna_data[8] =~ /ID=([^;]+)/ ) {
   $mrna_id = $1;
  }
  else {
   confess "Unkown GFF format for $record\n$mrna_data[8]";
  }
  $counter++;

  if ( $mrna_line =~ /coverage=([\d\.]+)/ ) {
   my $coverage = $1;
   if ( $coverage < 90 ) {
    print LOG "$mrna_id rejected: coverage ($coverage) below 90%.\n";
    $failed_cutoff++;
    next;
   }
  }
  if ( $mrna_line =~ /identity=([\d\.]+)/ ) {
   my $identity = $1;
   if ( $identity < $main::identical_fraction_cutoff ) {
    print LOG
      "$mrna_id rejected: identity ($identity) below $main::identical_fraction_cutoff.\n";
    $failed_cutoff++;
    next;
   }
  }
  my $gene_seq = $gene_seqs_hashref->{$mrna_id};
  my $pep_seq  = $pep_seqs_hashref->{$mrna_id};
  unless ( $gene_seq && $pep_seq ) {
   print LOG
"$mrna_id rejected: cannot find protein or genome sequence. Possible CDS is unavailable.\n";
   $failed_cutoff++;
   next;
  }

  if ( $pep_seq !~ /^M/ ) {
   print LOG "$mrna_id rejected: pep query does not start with M:"
   . substr( $pep_seq, 0, 1 ) . "\n";;
   $failed_cutoff++;
   next;
  }
  if ( $pep_seq !~ /\*$/ ) {
   print LOG "$mrna_id rejected: pep query does not end with a star (stop codon): "
   . substr( $pep_seq, -1, 1 ) . "\n";;
   $failed_cutoff++;
   next;
  }

  my $number_of_ns = $gene_seq =~ tr/N/N/;
  if ( $number_of_ns > ( length($gene_seq) * 0.40 ) ) {
   $failed_cutoff++;
   print LOG
     "$mrna_id rejected: genome sequence is more than 40% Ns ($number_of_ns)\n";
   next;
  }

  foreach my $other_lines (@record_data) {
   next unless $other_lines =~ /splice/;
   if ( $other_lines =~ /acceptor_site=([A-Z]+)/ ) {
    my $site = $1;
    unless ( $site eq 'AG' ) {
     print LOG
       "$mrna_id rejected: non-canonical acceptor splice site ($site)\n";
     $failed_cutoff++;
     next GENE;
    }
   }
   elsif ( $other_lines =~ /donor_site=([A-Z]+)/ ) {
    my $site = $1;
    unless ( $site eq 'GC' || $site eq 'GT' ) {
     print LOG "$mrna_id rejected: non-canonical donor splice site ($site)\n";
     $failed_cutoff++;
     next GENE;
    }
   }
  }
  $accepted{$mrna_id} = $gene_data[1];
  print OUT $record . $delimiter;
 }
 close OUT;
 close GFF;
 $/ = $orig_sep;

 open( PASS, ">$gff_file.passed" );
 my @shuffled_genes = shuffle( keys %accepted );
 foreach my $gene (@shuffled_genes) {
  print PASS $gene . "\t" . $accepted{$gene} . "\n";
 }
 close PASS;

 my $number_of_passing_genes = scalar(@shuffled_genes);
 print LOG
"Processed $counter transcripts.\nFound $number_of_passing_genes sequences passing criteria.\n\n";
 close LOG;

 print
"Processed $counter transcripts.\nFound $number_of_passing_genes sequences passing criteria.\n\n";

 #cleanup
 unlink( $gff_file . '.pep' );
 unlink( $gff_file . '.gene' );
 unlink( $gff_file . '.cds' );

}

###############################################################################
# predict_orfs  (v1 lines 3375-3388)
###############################################################################
sub predict_orfs {

 # will not use PFAM
 $ENV{PATH} .= ":$main::RealBin/../3rd_party/transdecoder";
 my $fasta               = $main::mrna_file;
 my ($transdecoder_exec) = Golden::Filter::check_program('TransDecoder.LongOrfs');
 my $cmd                 = $transdecoder_exec
   . " -t $fasta --workdir $fasta.transdecoder --CPU $main::threads ";
 Golden::Filter::process_cmd($cmd) unless -s "$fasta.transdecoder.pep";
 die "Failed to produce $fasta.transdecoder.pep\n"
   unless -s "$fasta.transdecoder.pep";
 return "$fasta.transdecoder.pep";

}

1;
