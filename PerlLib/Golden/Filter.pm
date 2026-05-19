package Golden::Filter;

use strict;
use warnings;
use Carp qw(confess);
use Exporter 'import';
use List::Util qw(shuffle);
use File::Basename;

use Nuc_translator;
use Gene_obj;
use Gene_obj_indexer;
use GFF3_utils;

our @EXPORT_OK = qw(
    validate_gene_structure
    _splice_cds
    _revcomp
    _translate
    filter_gff
    sort_gff3
    get_gff_delimiter
    order_fasta
    process_cmd
    gff3_fix_phase
    remove_overlapping_gff
    parse_genome_gff
    get_id_seq_from_fasta
    read_fasta
    bed_to_gff3
    check_program
    check_program_optional
    check_sort_version
    partition_genes
);

# ---------------------------------------------------------------------------
# Pure-Perl gene-structure validator. Replaces the fathom-based validator at
# lines 1054-1098 of v1 prepare_golden_genes_for_predictors.pl.
#
# Returns ('OK', '') if the gene is valid; ('FAIL', <reason>) otherwise.
# Caller is expected to remove rejected genes from the accepted set and
# append the reason to gene_validation.log.
# ---------------------------------------------------------------------------
sub validate_gene_structure {
    my ($gene_features, $genome_seq, $opts) = @_;
    my $require_complete = $opts->{complete}         // 0;
    my $reject_ambiguous = $opts->{reject_ambiguous} // 1;

    # 1. Collect CDS pieces in genomic-ascending order.
    my @cds = sort { $a->{start} <=> $b->{start} }
              grep { $_->{type} eq 'CDS' } @$gene_features;
    return ('FAIL', 'no_CDS') unless @cds;

    my $liberal          = $opts->{liberal}          // 0;

    # 2. Canonical splice sites: every intron between consecutive CDS pieces
    #    must be GT..AG (or GC..AG on the donor side for U12) read in
    #    transcript orientation. Skipped under liberal mode.
    unless ($liberal) {
        for (my $i = 0; $i < $#cds; $i++) {
            my $intron_start = $cds[$i]->{end} + 1;
            my $intron_end   = $cds[$i + 1]->{start} - 1;
            next if $intron_end < $intron_start;            # adjacent CDS = no real intron
            my $donor    = substr($genome_seq, $intron_start - 1, 2);
            my $acceptor = substr($genome_seq, $intron_end - 2,   2);
            if ($cds[0]->{strand} eq '-') {
                ($donor, $acceptor) = (_revcomp($acceptor), _revcomp($donor));
            }
            unless (($donor eq 'GT' || $donor eq 'GC') && $acceptor eq 'AG') {
                return ('FAIL', "non_canonical_splice:${donor}-${acceptor}");
            }
        }
    }

    # 3. CDS length divisible by 3.
    my $cds_seq = _splice_cds(\@cds, $genome_seq);
    return ('FAIL', 'cds_length_not_mod3') if length($cds_seq) % 3 != 0;

    # 4. No internal stops (skipped under liberal); no ambiguous codons
    #    unless caller opts out.
    my $aa = _translate($cds_seq);
    unless ($liberal) {
        my $internal_stops = () = $aa =~ /\*(?=.)/g;
        return ('FAIL', "internal_stop_count_$internal_stops") if $internal_stops > 0;
    }
    if ($reject_ambiguous) {
        my $x_count = () = $aa =~ /X/g;
        return ('FAIL', "ambiguous_codon_count_$x_count") if $x_count > 0;
    }

    # 5. Completeness check only when caller requested it (skipped under liberal).
    if ($require_complete && !$liberal) {
        return ('FAIL', 'missing_start') unless substr($cds_seq, 0, 3)  =~ /^ATG$/i;
        return ('FAIL', 'missing_stop')  unless substr($cds_seq, -3)    =~ /^(TAA|TAG|TGA)$/i;
    }

    return ('OK', '');
}

sub _revcomp {
    my $s = shift;
    $s = reverse $s;
    $s =~ tr/ACGTacgtNn/TGCAtgcaNn/;
    return $s;
}

sub _splice_cds {
    my ($cds_arr, $genome_seq) = @_;
    confess "_splice_cds: empty cds_arr" unless @$cds_arr;
    my $strand = $cds_arr->[0]->{strand};
    confess "_splice_cds: mixed strands within one gene"
        if grep { $_->{strand} ne $strand } @$cds_arr;

    my @ordered = $strand eq '-'
        ? sort { $b->{start} <=> $a->{start} } @$cds_arr
        : sort { $a->{start} <=> $b->{start} } @$cds_arr;

    my $out = '';
    for my $piece (@ordered) {
        my $seg = substr($genome_seq, $piece->{start} - 1,
                                       $piece->{end} - $piece->{start} + 1);
        $seg = _revcomp($seg) if $strand eq '-';
        $out .= $seg;
    }
    return $out;
}

sub _translate {
    return Nuc_translator::translate_sequence(shift, 1);
}

# ---------------------------------------------------------------------------
# GFF utility subs lifted from v1. Used by the driver and by Augustus.pm
# during golden / training / test partitioning.
# ---------------------------------------------------------------------------

sub get_gff_delimiter {
    my $gff = shift;
    confess "GFF file $gff not found or empty" unless -s $gff;
    my $orig_sep = $/;
    open( my $gfh, '<', $gff ) or confess "Cannot open $gff: $!";
    my $delimiter;
    while (my $line = <$gfh>) {
        if ($line =~ /^\s*\n$/) { $delimiter = "\n\n"; last; }
        elsif ($line =~ /^\#\#\#\n$/) { $delimiter = "\n\#\#\#\n"; last; }
    }
    close $gfh;
    $/ = $orig_sep;
    return $delimiter // "\n\n";
}

sub sort_gff3 {
    my $gff = shift;
    my $delimiter = shift // get_gff_delimiter($gff);
    my $orig_sep = $/;
    open( my $gfh, '<', $gff ) or confess "Cannot open $gff: $!";
    open( my $out, '>', "$gff.sorted" ) or confess "Cannot write $gff.sorted: $!";
    $/ = $delimiter;
    my @records;
    while (my $rec = <$gfh>) {
        chomp $rec;
        next unless $rec =~ /\S/;
        my ($first) = split /\n/, $rec, 2;
        my @f = split /\t/, $first;
        push @records, [ $f[0] // '', $f[3] // 0, $rec ];
    }
    close $gfh;
    for my $r (sort { $a->[0] cmp $b->[0] || $a->[1] <=> $b->[1] } @records) {
        print $out $r->[2] . $delimiter;
    }
    close $out;
    rename( "$gff.sorted", $gff ) or confess "rename $gff.sorted -> $gff: $!";
    $/ = $orig_sep;
    return $gff;
}

# filter_gff: writes accepted records to <filter_out>; rejects to <filter_out2>.
# Records are matched on (gene_id OR mrna_id) AND source — v1 keeps a
# hashref of $id => $source from the .passed file, so an alignment from
# one source can't be substituted for another with the same id.
# Returns the number of rejected records.
sub filter_gff {
    my ($gff, $hash_ref, $filter_out, $filter_out2) = @_;
    return if !$gff || !$hash_ref || !-s $gff;
    $filter_out  //= "$gff.filtered";
    $filter_out2 //= "$filter_out.rest";

    my $delimiter = get_gff_delimiter($gff);
    my $orig_sep  = $/;
    $/ = $delimiter;

    open( my $gfh, '<', $gff ) or confess "Cannot open $gff: $!";
    open( my $ok,  '>', $filter_out )  or confess "Cannot write $filter_out: $!";
    open( my $rej, '>', $filter_out2 ) or confess "Cannot write $filter_out2: $!";
    my $rejects = 0;

GENE: while ( my $record = <$gfh> ) {
        chomp $record;
        next if !$record || $record =~ /^\s*$/;
        my ( $gene_id, $mrna_id );
        my @lines = split /\n/, $record;
        my @gf    = split /\t/, ($lines[0] // '');
        if ( $gf[8] && $gf[8] =~ /ID=([^;]+)/ ) {
            $gene_id = $1;
            my @mf = split /\t/, ($lines[1] // '');
            if ( $mf[8] && $mf[8] =~ /ID=([^;]+)/ ) {
                $mrna_id = $1;
            }
        }
        else {
            $gf[8] =~ s/\.mRNA$// if defined $gf[8];
            $gene_id = $gf[8];
        }
        confess "No gene ID for this record!\n$record\n" unless $gene_id;
        my $source = $gf[1] // '';
        if (   ($gene_id && $hash_ref->{$gene_id} && $hash_ref->{$gene_id} eq $source)
            || ($mrna_id && $hash_ref->{$mrna_id} && $hash_ref->{$mrna_id} eq $source) )
        {
            print $ok $record . $delimiter;
        }
        else {
            print $rej $record . $delimiter;
            $rejects++;
        }
    }
    close $gfh; close $ok; close $rej;
    warn "Warning: filtered file $filter_out is empty\n"  if !-s $filter_out;
    warn "Warning: rejected file $filter_out2 is empty\n" if !-s $filter_out2;
    $/ = $orig_sep;
    return $rejects;
}

# order_fasta: extract scaffolds named in $gff_file out of $genome_seq_file,
# producing $gff_file.fasta. Lifted from v1 sub order_fasta (line 1524).
sub order_fasta {
    my ($genome_seq_file, $gff_file) = @_;
    confess "order_fasta: genome file '$genome_seq_file' missing" unless -s $genome_seq_file;
    confess "order_fasta: gff file '$gff_file' missing"           unless -s $gff_file;
    open( my $gfh, '<', $gff_file ) or confess "Cannot open $gff_file: $!";
    my %seen;
    while (my $line = <$gfh>) {
        next if $line =~ /^\s*#/ || $line !~ /\S/;
        my @f = split /\t/, $line;
        $seen{$f[0]} = 1 if defined $f[0] && $f[0] ne '';
    }
    close $gfh;

    open( my $ifh, '<', $genome_seq_file ) or confess "Cannot open $genome_seq_file: $!";
    open( my $ofh, '>', "$gff_file.fasta" ) or confess "Cannot write $gff_file.fasta: $!";
    local $/ = "\n>";
    while (my $rec = <$ifh>) {
        chomp $rec;
        $rec =~ s/^>//;
        next unless $rec =~ /\S/;
        my ($hdr, @body) = split /\n/, $rec;
        my ($id) = $hdr =~ /^(\S+)/;
        next unless defined $id && $seen{$id};
        print $ofh ">$hdr\n", join("\n", @body), "\n";
    }
    close $ifh; close $ofh;
    return "$gff_file.fasta";
}

# ---------------------------------------------------------------------------
# process_cmd: command runner lifted verbatim from v1 (line 2280).
# ---------------------------------------------------------------------------
sub process_cmd {
 my ($cmd) = @_;
 print "CMD: $cmd\n" if $main::verbose;
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  die "Error, cmd died with ret $ret\n";
 }
 return $ret;
}

# ---------------------------------------------------------------------------
# gff3_fix_phase: index a GFF3, recompute CDS phases via Gene_obj, emit fixed
# gff/pep/cds/gene files. Lifted from v1 (line 2290).
# ---------------------------------------------------------------------------
sub gff3_fix_phase {
 my $gff3_file = shift;
 open( IN, $gff3_file ) || confess( "Cannot find $gff3_file " . $! );
 my $index_file = "$gff3_file.inx";
 my $gene_obj_indexer = new Gene_obj_indexer( { "create" => $index_file } );
 my $asmbl_id_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs( $gff3_file, $gene_obj_indexer );
 die "There was an error indexing the gff file $gff3_file" unless $asmbl_id_to_gene_list_href && scalar(keys %{$asmbl_id_to_gene_list_href})>0;

 open( OUT,  ">$gff3_file.gff3" );
 open( PEP,  ">$gff3_file.pep" );
 open( CDS,  ">$gff3_file.cds" );
 open( GENE, ">$gff3_file.gene" );

 foreach my $asmbl_id ( sort keys %$asmbl_id_to_gene_list_href ) {
  my $genome_seq = $main::scaffold_seq_hashref->{$asmbl_id};
  if ( !$genome_seq ) {
   warn "Cannot find sequence $asmbl_id from genome\n";
   next;
  }
  my @gene_ids = @{ $asmbl_id_to_gene_list_href->{$asmbl_id} };
  foreach my $gene_id (@gene_ids) {
   my %params;
   my %preferences;
   $preferences{'sequence_ref'} = \$genome_seq;
   $params{unspliced_transcript} = 1;

   my $gene_obj_ref = $gene_obj_indexer->get_gene($gene_id);

   $gene_obj_ref->create_all_sequence_types( \$genome_seq, %params );

   foreach
     my $isoform ( $gene_obj_ref, $gene_obj_ref->get_additional_isoforms() )
   {

    $isoform->delete_isoforms();
    my $isoform_id = $isoform->{Model_feat_name};
    my @model_span = $isoform->get_CDS_span();
    next
      if ( !$isoform->get_CDS_span()
           || abs( $model_span[0] - $model_span[1] ) < 3 );

    eval { $isoform->set_CDS_phases( \$genome_seq ); };

    # get sequences
    # CDS
    my $seq = $isoform->get_CDS_sequence();
    $seq =~ s/(\S{60})/$1\n/g if $seq;
    chomp $seq if $seq;
    if ( $seq && length($seq) >= $main::minorf ) {
     print CDS ">$isoform_id $gene_id\n$seq\n";
    }
    else {
     next;
    }

    # proteins
    $seq = $isoform->get_protein_sequence();
    $seq =~ s/(\S{60})/$1\n/g;
    chomp $seq;
    print PEP ">$isoform_id $gene_id\n$seq\n";

    # gene
    $seq = $isoform->get_gene_sequence();
    $seq =~ s/(\S{60})/$1\n/g;
    chomp $seq;
    print GENE ">$isoform_id $gene_id\n$seq\n";

    # GFF3
    print OUT $isoform->to_GFF3_format_extended(%preferences) . "\n";

   }
  }
 }
 close OUT;
 close PEP;
 close CDS;
 close GENE;
 rename( $gff3_file, "$gff3_file.original" );
 unlink $index_file;
 rename( "$gff3_file.gff3", $gff3_file );
 &sort_gff3("$gff3_file");
 return;
}

# ---------------------------------------------------------------------------
# remove_overlapping_gff: keep the longest-CDS gene per overlap cluster on
# each ref/strand. Lifted from v1 (line 1618).
# ---------------------------------------------------------------------------
sub remove_overlapping_gff {

 # we allow genes on opposite strands
 my $file               = shift;
 my $out                = shift;
 my $delimiter          = shift;
 my $optional_pass_file = shift;

 $delimiter = &get_gff_delimiter($file) if !$delimiter;

 &sort_gff3( $file, $delimiter );

 my ( %overlap_found, %master_gene_list, %cds_sizes );
 my $skipped = int(0);
 my ( $previous_ref, $previous_mrna_id, %last_position_check );

 open( LOG2, ">$out.log" );
 my $orig_sep = $/;
 $/ = $delimiter;

 open( IN, $file );

GENE: while ( my $record = <IN> ) {
  chomp($record);
  my ( $gene_id, $mrna_id );
  next if !$record || $record =~ /^\s*$/;
  my @record_data = split( "\n", $record );
  my @gene_data   = split( "\t", $record_data[0] );
  confess "Weird GFF for $record" unless $gene_data[8];
  my $ref_id = $gene_data[0];
  my $strand = $gene_data[6];

  if ( $gene_data[8] =~ /ID=([^;]+)/ ) {
   $gene_id = $1;
  }
  else {

   # gene_id
   $gene_data[8] =~ s/\.mRNA$//;
   $gene_id = $gene_data[8];
  }
  my @mrna_data = split( "\t", $record_data[1] );
  if ( $mrna_data[8] && $mrna_data[8] =~ /ID=([^;]+)/ ) {
   $mrna_id = $1;
   $mrna_id =~ s/\.mRNA$//;
  }
  else {
   $mrna_id = $gene_id;
  }
  confess "No gene ID for this record!\n$record\n" unless $gene_id;

  unless ($mrna_id) {
   print LOG2 "No mRNA found for gene $mrna_id skipped\n";
   $skipped++;
   next GENE;
  }
  if ( $master_gene_list{$mrna_id} ) {
   print LOG2 "Gene $mrna_id found more than once. Skipping new data and keeping what I found first\n";
   $skipped++;
   next GENE;
  }

  foreach my $d (@record_data) {
   my @t = split( "\t", $d );
   if ( $t[2] && $t[2] eq 'CDS' ) {
    $cds_sizes{$mrna_id} += abs( $t[4] - $t[3] );
   }
  }
  if ( !$cds_sizes{$mrna_id} ) {
   print LOG2 "Non-coding transcript $mrna_id skipped\n";
   $skipped++;
   next GENE;
  }

  my ( $smallest_coord, $largest_coord ) =
    sort { $a <=> $b } ( $gene_data[3], $gene_data[4] );

  if (    $previous_mrna_id
       && $previous_ref
       && $last_position_check{$ref_id}{$strand}
       && $previous_ref eq $ref_id
       && $smallest_coord <= $last_position_check{$ref_id}{$strand}
       )
  {
   print LOG2 "Gene overlapping $previous_mrna_id found: $mrna_id\n";
   $overlap_found{$previous_mrna_id}{$mrna_id} = $record . $delimiter;
   next;
  }
  $master_gene_list{$mrna_id}            = $record . $delimiter;
  $previous_mrna_id                      = $mrna_id;
  $previous_ref                          = $ref_id;
  $last_position_check{$ref_id}{$strand} = $largest_coord;
 }
 $/ = $orig_sep;
 close IN;
 print LOG2 "\nDeciding on overlaps:\n";

 my %kept;

 open( OUT, ">$out" );
 foreach my $gene ( sort keys %master_gene_list ) {
  if ( $overlap_found{$gene} ) {
   my $master_size = $cds_sizes{$gene};
   my $gene_longest;
   foreach my $overlap_gene ( keys %{ $overlap_found{$gene} } ) {
    my $size = $cds_sizes{$overlap_gene};
    $gene_longest = $overlap_gene if $size > $master_size;
   }
   if ( !$gene_longest ) {
    print LOG2 "$gene overlaps: $gene kept as longest CDS\n";
    print OUT $master_gene_list{$gene};
    my @lines = split( "\n", $master_gene_list{$gene} );
    my @data  = split( "\t", $lines[0] );
    $kept{$gene} = $data[1];
   }
   else {
    print LOG2 "$gene overlaps: longer CDS $gene_longest found\n";
    print OUT $overlap_found{$gene}{$gene_longest};
    my @lines = split( "\n", $overlap_found{$gene}{$gene_longest} );
    my @data  = split( "\t", $lines[0] );
    $kept{$gene_longest} = $data[1];
   }
   $skipped += scalar( keys %{ $overlap_found{$gene} } );
  }
  else {

   # no overlaps
   print OUT $master_gene_list{$gene};
   my @lines = split( "\n", $master_gene_list{$gene} );
   my @data  = split( "\t", $lines[0] );
   $kept{$gene} = $data[1] if !$kept{$gene};    # only the first one
  }

 }
 close OUT;
 close LOG2;
 &sort_gff3($out);
 my %passed;
 if ( $optional_pass_file && -s $optional_pass_file ) {
  open( IN, $optional_pass_file );
  while ( my $ln = <IN> ) {
   chomp($ln);
   my @data = split( "\t", $ln );
   if ( $kept{ $data[0] } && $kept{ $data[0] } eq $data[1] ) {
    $passed{ $data[0] } = $data[1];
   }
  }
  close IN;

  open( OUT, ">$optional_pass_file" );
  foreach my $id ( keys %passed ) {
   print OUT $id . "\t" . $passed{$id} . "\n";
  }
  close OUT;

 }

 print "\tOverlaps checked. Kept ".scalar(keys %kept)." genes. See $out.log for details\n";
}

# ---------------------------------------------------------------------------
# parse_genome_gff: build a per-scaffold filter file for PASA-style genome
# GFFs. Lifted from v1 (line 2467).
# ---------------------------------------------------------------------------
sub parse_genome_gff {
 my %hash;
 my ( $skipped, %allowed_data );
 my $gff_file = shift;
 die unless $gff_file && -s $gff_file;
 my $program = shift;
 $program = 'PASA' unless $program;
 &remove_overlapping_gff( $gff_file, $gff_file . ".nr" );
 $gff_file .= '.nr';

 open( IN, $gff_file ) || die("Cannot find $gff_file\n");
 while ( my $ln = <IN> ) {
  next if $ln =~ /^#/ || $ln =~ /^\s*$/ || $ln !~ /\tmRNA\t/;
  chomp($ln);
  my @data = split( "\t", $ln );
  next unless $data[8];
  my $mrna_id;
  if ( $data[8] =~ /ID=([^;]+)/ ) {
   $mrna_id = $1;
   $hash{ $data[0] }{$mrna_id} = $ln;
  }
  if ($main::no_single_exon) {
   my $exon_counting = int(0);
   while ( my $ln2 = <IN> ) {
    last if $ln2 =~ /^#/ || $ln2 =~ /^\s*$/;
    $exon_counting++ if $ln2 =~ /\texon\t/;
   }
   if ( $exon_counting <= 1 ) {
    delete( $hash{ $data[0] }{$mrna_id} );
    $skipped++;
   }
  }
 }
 close IN;

 mkdir $main::genome_dir unless -d $main::genome_dir;
 foreach my $scaffold ( keys %hash ) {
  open( OUT, ">$main::genome_dir/$scaffold.gff.filter" )
    or confess "Cannot write $main::genome_dir/$scaffold.gff.filter: $!";
  print OUT "\t$scaffold\tProduced by $gff_file and $program\n";
  foreach my $mrna_id ( keys %{ $hash{$scaffold} } ) {
   my @data = split( "\t", $hash{$scaffold}{$mrna_id} );
   my $orient = $data[6] eq '+' ? 1 : 0;
   print OUT $data[3] . "\t"
     . $data[4]
     . "\t999\t0\t0\t$orient\t0\t0\t$mrna_id\n";
   $allowed_data{$mrna_id} = $data[1];
  }
  close OUT;
 }
 print "\tSkipped $skipped single exon genes\n" if $skipped;

 # update to non-overlapping file in case we need to use it again.
 $main::pasa_genome_gff = $gff_file;
 return \%allowed_data;
}

# ---------------------------------------------------------------------------
# get_id_seq_from_fasta: cdbyank-backed single-record FASTA fetcher with a
# memoising hash. Lifted from v1 (line 2443).
# ---------------------------------------------------------------------------
sub get_id_seq_from_fasta {
 my ( $acc, $fasta_db ) = @_;
 my $seq;
 if ( !$main::get_id_seq_from_fasta_hash{$fasta_db}{$acc} ) {
  $seq = `$main::cdbyank_exec $fasta_db.cidx -a '$acc'`;
  chomp($seq);
  $main::get_id_seq_from_fasta_hash{$fasta_db}{$acc} = $seq;
 }
 else {
  $seq = $main::get_id_seq_from_fasta_hash{$fasta_db}{$acc};
 }
 unless ($seq) {

  #warn "WARNING: couldn't retrieve seq for $acc from $fasta_db\n";
  return;
 }
 my @x = split( /\n/, $seq );
 my $id = shift @x;
 $seq = join( "", @x );
 $seq = uc($seq);
 $seq =~ s/\s//g;
 return $seq;
}

# ---------------------------------------------------------------------------
# read_fasta: slurp a FASTA into ( \%seq, \%length ). Lifted from v1 (3390).
# ---------------------------------------------------------------------------
sub read_fasta {
 my $fasta = shift;
 my (%hash,%hash_length);
 my $orig_sep = $/;
 $/ = '>';
 open( IN, $fasta ) || confess( "Cannot open $fasta : " . $! );
 while ( my $record = <IN> ) {
  chomp($record);
  next unless $record;
  my @lines = split( "\n", $record );
  my $id    = shift(@lines);
  my $seq   = join( '', @lines );
  $seq =~ s/\s+//g;
  if ( $id && $seq && $id =~ /^(\S+)/ ) {
   $hash{$1} = $seq;
   $hash_length{$1} = length($seq);
  }
 }
 close IN;
 $/ = $orig_sep;
 return (\%hash,\%hash_length);
}

# ---------------------------------------------------------------------------
# bed_to_gff3: BED12 -> GFF3 (gene_obj). Lifted from v1 (3514).
# ---------------------------------------------------------------------------
sub bed_to_gff3 {
	# from $JAMG_PATH/3rd_party/PASA/misc_utilities/bed_to_gene_gff3.pl
	my $input_bed = shift;
	my $outfile = $input_bed;
	$outfile=~s/\.bed$//;$outfile.='.gff3';
	print "Converting $input_bed BED to $outfile GFF\n";
	if (-s $outfile){
		print "GFF $outfile already exists. Will use it. Stop and delete it otherwise\n";
		sleep(2);
		return $outfile;
	}
	open (OUT,">$outfile");
	my $counter = 0;
	open (my $fh, $input_bed) or die "Error, cannot open file $input_bed";
	while (<$fh>) {
		my @x = split(/\t/);
		next unless $x[5];
		my $scaff = $x[0];
		my $gene_lend = $x[1] + 1;
		my $gene_rend = $x[2];
		my $com_name = $x[3];
		my $score = $x[4];
		my $orient = $x[5];
		$orient = '+' if ($orient eq '*');
		my $coding_lend = $x[6] + 1;
		my $coding_rend = $x[7];
		my $rgb_color = $x[8];

		my $num_exons = $x[9];

		my $lengths_text = $x[10];
		my $exon_relative_starts_text = $x[11];

		my @lengths = split(/,/, $lengths_text);
		my @exon_relative_starts = split(/,/, $exon_relative_starts_text);

		my @exons;

		my $sum_len = 0;

		while (@lengths) {
			my $len = shift @lengths;
			my $start = shift @exon_relative_starts;
			my $exon_lend = $gene_lend + $start;
			my $exon_rend = $exon_lend + $len - 1;
			push (@exons, [$exon_lend, $exon_rend]);
			$sum_len += $len;
		}

		if ($sum_len < 3) {
			next;
		}

		eval {

			my $gene_obj = new Gene_obj();

            if ($coding_lend == $coding_rend +1) { ## not coding
                $coding_lend = 0;
                $coding_rend = 0;
            }

			$gene_obj->build_gene_obj_exons_n_cds_range(\@exons, $coding_lend, $coding_rend, $orient);
			if ($com_name && $com_name=~/ID=([^;]+)/){
				#'ID=asmbl_40814|m.42780;Name=ORF_(asmbl_40814|g.42780,_asmbl_40814|m.42780);'
				$gene_obj->{Model_feat_name} = $1;
				if ($com_name=~/Name=([^;]+)/){
					$gene_obj->{com_name} = $1;
					if ($gene_obj->{com_name}=~/(asmbl_\d+\|g\.\d+)/){
						$gene_obj->{TU_feat_name} = $1;
					}
				}
			}else{
				$gene_obj->{Model_feat_name} = "model.$counter";
				$gene_obj->{TU_feat_name} = "gene.$counter";
				$gene_obj->{com_name} = $com_name;
			}
			$gene_obj->{asmbl_id} = $scaff;
			$counter++;
			print OUT $gene_obj->to_GFF3_format() . "\n";

		};
	}
	close OUT;
	print "Error when conveting to GFF $outfile from BED $input_bed:\n$@\n" unless -s $outfile;
	return $outfile;
}

# ---------------------------------------------------------------------------
# check_program / check_program_optional: which-based binary locators.
# Lifted from v1 (2416 / 2429).
# ---------------------------------------------------------------------------
sub check_program {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  die "Error, path to required $prog cannot be found\n"
    unless $path =~ /^\//;
  chomp($path);
  push( @paths, $path );
 }
 return @paths;
}

sub check_program_optional {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  warn
"Warning: path to optional $prog cannot be found in your path environment.\n"
    unless $path =~ /^\//;
  chomp($path);
  push( @paths, $path );
 }
 return @paths;
}

# ---------------------------------------------------------------------------
# check_sort_version: probe GNU sort for --parallel support. Lifted (3819).
# ---------------------------------------------------------------------------
sub check_sort_version {
	my ($sort_exec) = &check_program('sort');
	my $cpus = 4;
	my @v=`$sort_exec --version`;

	if ($v[0] && $v[0]=~/(\d+)\.(\d+)\s*$/){
		my $major = $1;
		my $minor = $2;
		if ($major >= 8 && $minor >= 6){
			return "$sort_exec -T $main::tmpdir --parallel $cpus -S $main::sort_buffer";
		}else{
			return "$sort_exec -S $main::sort_buffer -T $main::tmpdir";
		}
	}else{
		die "Sort of coreutils not found!";
	}
}

# ---------------------------------------------------------------------------
# partition_genes: golden/training/test/optimisation partitioning. Replaces
# v1's process_for_gene_prediction (lines 1034-1286) MINUS the fathom (1054-1097),
# geneid (1224-1257), SNAP (1259-1267) and GlimmerHMM (1268-1276) blocks. The
# fathom validator is replaced by the pure-Perl validate_gene_structure call
# wired in here.
#
# Signature:
#   partition_genes($final_gff_nr, $passed_file, $training_set_size,
#                   $augustus_flank_region, $augustus_dir,
#                   $genome_sequence_file, $cdbyank_exec, $gff3_introns_exec,
#                   $outdir, $require_complete)
# Returns: path to the cleaned golden GFF.
# ---------------------------------------------------------------------------
sub partition_genes {
    my ($final_gff_nr, $passed_file, $training_set_size,
        $augustus_flank_region, $augustus_dir,
        $genome_sequence_file, $cdbyank_exec_arg, $gff3_introns_exec_arg,
        $outdir, $require_complete) = @_;

    confess "partition_genes: input GFF '$final_gff_nr' missing or empty"
        unless $final_gff_nr && -s $final_gff_nr;
    confess "partition_genes: pass file '$passed_file' missing or empty"
        unless $passed_file && -s $passed_file;
    confess "partition_genes: outdir '$outdir' not a directory"
        unless $outdir && -d $outdir;

    print "Processing for gene prediction software...\n";
    my $gff_file = $final_gff_nr;
    &order_fasta( $genome_sequence_file, $gff_file );

    my ( %accepted, %training_genes, %aug_optimization_genes );
    open( my $pass_fh, '<', $passed_file )
        || die "Cannot find pass file $passed_file";
    while ( my $ln = <$pass_fh> ) {
        chomp($ln);
        my @data = split( "\t", $ln );
        next unless $data[1];
        $accepted{ $data[0] } = $data[1];
    }
    close $pass_fh;

    print "Creating golden set using main GFF3 file $gff_file\n";

    # -------- Pure-Perl validator (replaces v1 fathom block, lines 1054-1097)
    # Index CDS features per gene/mRNA from the .nr GFF.
    my %cds_by_gene;       # mrna_id => { scaffold => $scaff, cds => \@pieces }
    {
        my $delimiter = &get_gff_delimiter($gff_file);
        my $orig_sep  = $/;
        open( my $gfh, '<', $gff_file ) or confess "Cannot open $gff_file: $!";
        $/ = $delimiter;
        while ( my $record = <$gfh> ) {
            next unless $record =~ /\S/;
            my @lines = split /\n/, $record;
            my ( $mrna_id, $scaffold );
            my @cds;
            for my $line (@lines) {
                next if $line =~ /^\s*#/ || $line !~ /\S/;
                my @f = split /\t/, $line;
                next unless @f >= 9 && defined $f[2];
                if ( $f[2] eq 'mRNA' && $f[8] =~ /ID=([^;]+)/ ) {
                    $mrna_id  = $1;
                    $scaffold = $f[0];
                }
                elsif ( $f[2] eq 'CDS' ) {
                    push @cds, {
                        type   => 'CDS',
                        start  => $f[3] + 0,
                        end    => $f[4] + 0,
                        strand => $f[6],
                    };
                    $scaffold ||= $f[0];
                }
            }
            next unless $mrna_id && @cds;
            $cds_by_gene{$mrna_id} = { scaffold => $scaffold, cds => \@cds };
        }
        close $gfh;
        $/ = $orig_sep;
    }

    # Load genome scaffolds once.
    my ( $scaffold_seq, undef ) = &read_fasta($genome_sequence_file);

    my $validation_log = "$outdir/gene_validation.log";
    open( my $vlog, '>', $validation_log )
        or confess "Cannot write $validation_log: $!";
    my $rejected_count = 0;
    foreach my $gene_id ( keys %accepted ) {
        my $entry = $cds_by_gene{$gene_id};
        unless ($entry) {
            print $vlog "$gene_id\tFAIL\tno_CDS_in_nr_gff\n";
            delete $accepted{$gene_id};
            $rejected_count++;
            next;
        }
        my $scaff_seq = $scaffold_seq->{ $entry->{scaffold} };
        unless ( defined $scaff_seq ) {
            print $vlog "$gene_id\tFAIL\tscaffold_missing:$entry->{scaffold}\n";
            delete $accepted{$gene_id};
            $rejected_count++;
            next;
        }
        my ( $status, $reason ) = validate_gene_structure(
            $entry->{cds}, $scaff_seq,
            { complete => $require_complete,
              liberal  => $main::liberal_cutoffs }
        );
        if ( $status ne 'OK' ) {
            print $vlog "$gene_id\tFAIL\t$reason\n";
            delete $accepted{$gene_id};
            $rejected_count++;
        }
        else {
            print $vlog "$gene_id\tOK\t\n";
        }
    }
    close $vlog;
    print "Pure-Perl validator rejected $rejected_count of "
        . ( scalar(keys %accepted) + $rejected_count )
        . " genes; see $validation_log\n";

    # -------- Training / optimisation / test partitioning (v1 lines 1101-1124)
    my $number_of_passing_genes      = scalar( keys %accepted );
    my $passed_genes_in_training     = int(0);
    my $passed_genes_in_optimization = int(0);

    $training_set_size = int( $number_of_passing_genes * 0.40 )
        if !$training_set_size
        || $training_set_size > int( $number_of_passing_genes * 0.40 );
    foreach my $gene ( shuffle( keys %accepted ) ) {
        $training_genes{$gene} = $accepted{$gene};
        $passed_genes_in_training++;
        last if $passed_genes_in_training >= $training_set_size;
    }
    foreach my $gene ( shuffle( keys %accepted ) ) {
        next if $training_genes{$gene};
        $aug_optimization_genes{$gene} = $accepted{$gene};
        $passed_genes_in_optimization++;
        last if $passed_genes_in_optimization >= $main::aug_optimization_geneset;
    }
    $main::aug_optimization_geneset = scalar( keys %aug_optimization_genes );
    $training_set_size              = scalar( keys %training_genes );
    die "Augustus optimization gene set is zero!\n"
        if $main::aug_optimization_geneset < 1;
    die "Augustus Training gene set is zero!\n" if $training_set_size < 1;

    print "'Golden': $number_of_passing_genes. Filtering GFFs...\n";
    &filter_gff( $gff_file, \%accepted, "$gff_file.golden.gff3" );
    # filter_gff matches on (gene_id|mrna_id) AND source. If gff3_fix_phase or
    # any upstream pass rewrote column 2 (source) after %accepted was built
    # from .passed, every record can silently fall into .golden.gff3.rest and
    # the golden file ends up empty. Surface that loudly instead of letting
    # downstream gff2gbSmallDNA.pl produce mystery-empty .gb files.
    confess "partition_genes: $gff_file.golden.gff3 came out empty even though "
          . "$number_of_passing_genes genes survived validation. Likely a "
          . "source-column mismatch between %accepted (from $passed_file) and "
          . "$gff_file (column 2). Check gff3_fix_phase / sort_gff3 ordering."
        unless -s "$gff_file.golden.gff3";
    my @to_evaluate;
    push( @to_evaluate, "$gff_file.golden.gff3" );

    print
"Processing for gene predictors... - Specific filters may reject some more genes during checking/conversion\n";
    print "\nProcessing for:\n";

    # Augustus
    print "\tAugustus\n";

    # Hints for training (kept from v1; Augustus-shaped emitter).
    Golden::Augustus::gff2hints( "$gff_file.golden.gff3", 1 )
        if -s "$gff_file.golden.gff3";
    Golden::Augustus::gff2hints("$gff_file") if -s "$gff_file";
    Golden::Augustus::gff2hints("$gff_file.golden.gff3.rest")
        if -s "$gff_file.golden.gff3.rest";

    &filter_gff( "$gff_file.golden.gff3", \%training_genes,
                 "$gff_file.golden.train.gff3" );
    &filter_gff( "$gff_file.golden.train.gff3.rest",
                 \%aug_optimization_genes,
                 "$gff_file.golden.optimization.gff3" );
    rename( "$gff_file.golden.optimization.gff3.rest",
            "$gff_file.golden.test.gff3" );

    # Resolve Augustus helper paths from $augustus_dir.
    my $gff2gb_exec               = "$augustus_dir/scripts/gff2gbSmallDNA.pl";
    my $augustus_train_exec_local = "$augustus_dir/bin/etraining";
    my $augustus_filter_exec      = "$augustus_dir/scripts/filterGenes.pl";

    if (    -x $gff2gb_exec
         && -x $augustus_train_exec_local
         && -x $augustus_filter_exec )
    {
        &process_cmd(
"$gff2gb_exec $gff_file.golden.train.gff3 $genome_sequence_file $augustus_flank_region $gff_file.golden.train.gb >/dev/null 2> /dev/null"
        ) if -s "$gff_file.golden.train.gff3";
        &process_cmd(
"$augustus_train_exec_local --species=generic $gff_file.golden.train.gb 2>&1 | grep 'n sequence' | perl -pe 's/.*n sequence (\\S+):.*/\$1/' | sort -u > $gff_file.golden.train.gb.bad.lst 2>/dev/null"
        ) if -s "$gff_file.golden.train.gb";

        if ( -s "$gff_file.golden.train.gb.bad.lst" ) {
            &process_cmd(
"$augustus_filter_exec $gff_file.golden.train.gb.bad.lst $gff_file.golden.train.gb > $gff_file.golden.train.good.gb 2>/dev/null"
            );
        }
        else {
            symlink( "$gff_file.golden.train.gb",
                     "$gff_file.golden.train.good.gb" )
                if -s "$gff_file.golden.train.gb";
        }

        &process_cmd(
"$gff2gb_exec $gff_file.golden.test.gff3 $genome_sequence_file $augustus_flank_region $gff_file.golden.test.gb  >/dev/null 2>/dev/null"
        ) if -s "$gff_file.golden.test.gff3";
        &process_cmd(
"$augustus_train_exec_local --species=generic $gff_file.golden.test.gb 2>&1 | grep 'n sequence' | perl -pe 's/.*n sequence (\\S+):.*/\$1/' | sort -u > $gff_file.golden.test.gb.bad.lst 2>/dev/null"
        ) if -s "$gff_file.golden.test.gb";

        if ( -s "$gff_file.golden.test.gb.bad.lst" ) {
            &process_cmd(
"$augustus_filter_exec $gff_file.golden.test.gb.bad.lst $gff_file.golden.test.gb > $gff_file.golden.test.good.gb 2>/dev/null"
            );
        }
        else {
            symlink( "$gff_file.golden.test.gb",
                     "$gff_file.golden.test.good.gb" )
                if -s "$gff_file.golden.test.gb";
        }

        &process_cmd(
"$gff2gb_exec $gff_file.golden.optimization.gff3 $genome_sequence_file $augustus_flank_region $gff_file.golden.optimization.gb  >/dev/null"
        ) if -s "$gff_file.golden.optimization.gff3";
        &process_cmd(
"$augustus_train_exec_local --species=generic $gff_file.golden.optimization.gb 2>&1 | grep 'n sequence' | perl -pe 's/.*n sequence (\\S+):.*/\$1/' | sort -u > $gff_file.golden.optimization.gb.bad.lst 2>/dev/null"
        ) if -s "$gff_file.golden.optimization.gb";

        if ( -s "$gff_file.golden.optimization.gb.bad.lst" ) {
            &process_cmd(
"$augustus_filter_exec $gff_file.golden.optimization.gb.bad.lst $gff_file.golden.optimization.gb > $gff_file.golden.optimization.good.gb 2>/dev/null"
            );
        }
        else {
            symlink( "$gff_file.golden.optimization.gb",
                     "$gff_file.golden.optimization.good.gb" )
                if -s "$gff_file.golden.optimization.gb";
        }
    }

    # parse GB (kept; sibling Augustus module hosts parse_gb)
    my $train_ref = Golden::Augustus::parse_gb("$gff_file.golden.train.good.gb")
        if -s "$gff_file.golden.train.good.gb";
    my $test_ref  = Golden::Augustus::parse_gb("$gff_file.golden.test.good.gb")
        if -s "$gff_file.golden.test.good.gb";
    my $opt_ref   = Golden::Augustus::parse_gb("$gff_file.golden.optimization.good.gb")
        if -s "$gff_file.golden.optimization.good.gb";

    # Introns -> gmap (v1 line 1279).
    &process_cmd(
        "$gff3_introns_exec_arg < $gff_file.golden.gff3 > $gff_file.golden.splice.gmap"
    ) if -s "$gff_file.golden.gff3" && $gff3_introns_exec_arg;

    print "\tevaluation\n";
    foreach my $file (@to_evaluate) {
        Golden::Augustus::gff_to_gtf($file);
    }

    return "$gff_file.golden.gff3";
}

1;
