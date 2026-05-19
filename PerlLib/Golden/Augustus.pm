package Golden::Augustus;

use strict;
use warnings;
use Carp qw(confess);
use File::Basename;
use Exporter 'import';

# Core Perl utilities used by lifted subs
use List::Util qw(shuffle);
use Pod::Usage;
use Bio::SeqIO;

# JAMg PerlLib modules used by lifted subs
use Gene_obj;
use Gene_obj_indexer;
use GFF3_utils;
use Nuc_translator;
use Fasta_reader;

our @EXPORT_OK = qw(
    wrap_text
    count_seq
    gff_to_gtf
    shuffle_fasta
    split_fasta
    check_augustus
    parse_gb
    gb2gff3
    gff2hints
);

# All globals are owned by the driver (bin/prepare_golden_genes.pl) in package
# `main`. Lifted subs reference them via fully-qualified $main::var syntax.

# ---------------------------------------------------------------------------
# wrap_text: width-wrapping utility used by FASTA writers.
# Lifted verbatim from v1 prepare_golden_genes_for_predictors.pl lines 2783-2790.
# ---------------------------------------------------------------------------
sub wrap_text {
 my $string      = shift;
 my $wrap_length = shift;
 $wrap_length = 120 if !$wrap_length;
 $string =~ s/(.{0,$wrap_length})/$1\n/g;
 $string =~ s/\n{2,}/\n/;
 return $string;
}

# ---------------------------------------------------------------------------
# count_seq: count `>` headers in a FASTA file.
# Lifted verbatim from v1 lines 2792-2802.
# ---------------------------------------------------------------------------
sub count_seq {
 my $fasta = shift;
 my $count = int(0);
 open( IN, $fasta );
 while ( my $ln = <IN> ) {
  next unless $ln =~ /^>/;
  $count++;
 }
 close IN;
 return $count;
}

# ---------------------------------------------------------------------------
# gff_to_gtf: convert GFF3 to GTF using Gene_obj_indexer + GFF3_utils.
# Lifted from v1 lines 2374-2414.
# References $main::scaffold_seq_hashref for genome sequence lookup.
# ---------------------------------------------------------------------------
sub gff_to_gtf {

 my $gff3_file = shift;
 open( OUT, ">$gff3_file.gtf" );

 my $inx_file = "$gff3_file.inx";
 my $gene_obj_indexer = new Gene_obj_indexer( { "create" => $inx_file } );
 my $asmbl_id_to_gene_list_href = &GFF3_utils::index_GFF3_gene_objs( $gff3_file, $gene_obj_indexer );
 die "There was an error indexing the data!" unless $asmbl_id_to_gene_list_href && scalar(keys %{$asmbl_id_to_gene_list_href})>0;

 foreach my $asmbl_id ( sort keys %$asmbl_id_to_gene_list_href ) {

  ## get the genome sequence
  my $genome_seq = $main::scaffold_seq_hashref->{$asmbl_id};
  if ( !$genome_seq ) {
   warn "Cannot find sequence $asmbl_id\n";
   next;
  }
  my @gene_ids = @{ $asmbl_id_to_gene_list_href->{$asmbl_id} };
  foreach my $gene_id (@gene_ids) {
   ## note models of isoforms are bundled into the same gene object.
   my $gene_obj_ref = $gene_obj_indexer->get_gene($gene_id);
   foreach
     my $gene_obj ( $gene_obj_ref, $gene_obj_ref->get_additional_isoforms() )
   {
    $gene_obj->delete_isoforms();    # unbundle the model object!
    my $gtf_text = "";
    eval { $gtf_text = $gene_obj->to_GTF_format( \$genome_seq ); };
    if ($@) {

     # do it in pseudogene mode - if not then UTR is printed as CDS... not good!
     $gene_obj->{is_pseudogene} = 1;
     $gtf_text = $gene_obj->to_GTF_format( \$genome_seq );
    }
    print OUT "$gtf_text\n";
   }
  }
 }
 close OUT;
 unlink $inx_file;
}

# ---------------------------------------------------------------------------
# shuffle_fasta: randomize FASTA order for training-set splitting.
# Lifted verbatim from v1 lines 2076-2115. Uses wrap_text from this module.
# ---------------------------------------------------------------------------
sub shuffle_fasta {
 # shuffling is important to allow aligners process the data
 # in roughly equal time per thread
 my $fasta    = shift;
 my $orig_sep = $/;

 my %sequence_data;
 open( IN, $fasta );
 #stupid files with > in description line
 if ( $orig_sep){
	$/ = '>';
	my $record =<IN>;$record =<IN>;chomp($record);
  	my @data = split( "\n", $record );
	my $id = shift @data;
	if ($id && $id=~/^(\S+)(.+)/){
	  	$sequence_data{$1} = $2."\n".join( '', @data );
	}
 }

 $/ = "\n>";

 while ( my $record = <IN> ) {
  chomp($record);
  my @data = split( "\n", $record );
  next unless $record;
  my $id = shift @data;
	if ($id && $id=~/^(\S+)(.+)/){
	  $sequence_data{$1} = $2."\n".join( '', @data );
	}
 }

 open( OUT, ">$fasta.shuff" );
 foreach my $id ( shuffle( keys %sequence_data ) ) {
  my ($descr,$seq) = split("\n",$sequence_data{$id});
  print OUT $/ . $id . $descr."\n".&wrap_text( $seq );
 }
 $/ = $orig_sep;
 close OUT;
 return "$fasta.shuff";
}

# ---------------------------------------------------------------------------
# split_fasta: split a FASTA into per-N or per-scaffold files.
# Lifted verbatim from v1 lines 2117-2217. Uses shuffle_fasta and wrap_text
# from this module.
# ---------------------------------------------------------------------------
sub split_fasta {
 my @files;
 my $file2split         = shift;
 my $outdir             = shift;
 my $how_many_in_a_file = shift;
 my $pattern            = shift;
 my $shuffle            = shift;
 my $min_seq_size = shift;
 return unless $file2split && -s $file2split && $outdir && $how_many_in_a_file;
 return if -d $outdir;

 print "Splitting $file2split with $how_many_in_a_file sequences that are larger than $min_seq_size bp";
 print " matching $pattern" if $pattern;
 print " (shuffled)" if $shuffle;
 print "\n";

 mkdir($outdir) unless -d $outdir;
 my $filecount;
 my $seqcount = int(0);
 undef ($shuffle) if -s $file2split > 1e8;
 my $shuffled_file = $shuffle ? &shuffle_fasta($file2split) : $file2split;
 my $orig_sep = $/;
 # the issue is that uniref has > characters in the description line which is
 # breaking the parser
 open( FILE, $shuffled_file ) || die;

 # we have two choices. skip the first sequence or assume
 # it doesn't break the parser. we do latter
 if ($orig_sep){
  $/ = ">";
  my $record = <FILE>;$record = <FILE>; # first empty >
  chomp($record);
  my @lines = split( "\n", $record );
  my $id = shift @lines;
  if ((!$pattern || $id =~ /$pattern/) && $id =~ /^(\S+)/){
  	$id = $1;
  }
  my $seq = join( '', @lines );
  $seq =~ s/\s+//g;
  $seq =~ s/\*$//;    # for protein stop codons
  # we will only search those that are big enough.
  if ($id && $seq && length($seq)>=$min_seq_size){
	  $seqcount++;
	  if (   !$filecount
	       || $how_many_in_a_file == 1
	       || $seqcount > $how_many_in_a_file )
	  {
	   $seqcount = int(0);
	   $filecount++;
	   my $outfile = $how_many_in_a_file == 1 ? $id : $file2split . "_" . $filecount;
	   $outfile = $outdir . '/' . basename($outfile);
	   close(OUT);
	   open( OUT, ">$outfile" ) || die("Cannot open $outfile");
	   push( @files, $outfile );
	  }
	  print OUT ">$id\n" . &wrap_text($seq);
  }
 }

 $/ = "\n>";

 while (my $record = <FILE> ) {
  chomp($record);
  next unless $record;
  my @lines = split( "\n", $record );
  my $id = shift @lines;
  next if ( $pattern && $id !~ /$pattern/ );
  if ($id =~ /^(\S+)/){
  	$id = $1;
  }else{
	die "Cannot find ID for a sequence in $file2split ($id)\n$record";
  }
  my $seq = join( '', @lines );
  $seq =~ s/\s+//g;
  $seq =~ s/\*$//;    # for protein stop codons
  # we will only search those that are big enough.
  next unless $seq && length($seq)>=$min_seq_size;
  $seqcount++;

  if (   !$filecount
       || $how_many_in_a_file == 1
       || $seqcount > $how_many_in_a_file )
  {
   $seqcount = int(0);
   $filecount++;
   my $outfile =
     $how_many_in_a_file == 1 ? $id : $file2split . "_" . $filecount;
   $outfile = $outdir . '/' . basename($outfile);
   close(OUT);
   open( OUT, ">$outfile" ) || die("Cannot open $outfile");
   push( @files, $outfile );
  }
  print OUT ">$id\n" . &wrap_text($seq);
 }
 close(FILE);
 close(OUT);
 $/ = $orig_sep;
 unlink($shuffled_file) if $shuffle;
 print "Split $file2split into ".scalar(@files)." segments\n";
 return \@files;
}

# ---------------------------------------------------------------------------
# check_augustus: confirm Augustus install layout and resolve helper scripts.
# Lifted from v1 lines 2058-2074. Accesses $main::augustus_dir,
# $main::gff2gb_exec, $main::augustus_train_exec, $main::augustus_filterGenes_exec.
# ---------------------------------------------------------------------------
sub check_augustus {
 pod2usage
"Can't find the Augustus directory, easiest way is to create a symlink of the augustus executable somewhere in your PATH (e.g. \$HOME/bin) or add the Augustus bin directory in your PATH\n"
   unless $main::augustus_dir && -d $main::augustus_dir;
 $main::gff2gb_exec = $main::augustus_dir . '/scripts/gff2gbSmallDNA.pl'
   if !$main::gff2gb_exec;
 $main::augustus_train_exec = $main::augustus_dir . '/bin/etraining'
   if !$main::augustus_train_exec;
 $main::augustus_filterGenes_exec = $main::augustus_dir . '/scripts/filterGenes.pl'
   if !$main::augustus_filterGenes_exec;
 pod2usage "Can't find etraining from Augustus\n"
   if !-s $main::augustus_train_exec;
 pod2usage "Can't find gff2gbSmallDNA.pl from Augustus\n"
   if !-s $main::gff2gb_exec;
 pod2usage "Can't find filterGenes.pl from Augustus\n"
   if !-s $main::augustus_filterGenes_exec;
}

# ---------------------------------------------------------------------------
# parse_gb: parse a GenBank file produced by Augustus' gff2gbSmallDNA.pl into
# a hashref of per-locus exon coordinates and sequences.
# Lifted from v1 lines 1779-1812. Accesses $main::cdbfasta_exec; the
# `&process_cmd` call is routed to Golden::Filter::process_cmd per spec.
# ---------------------------------------------------------------------------
sub parse_gb {
 my $file      = shift;
 my $fasta_out = $file . ".fasta";
 my %hash;
 my $gb_obj = Bio::SeqIO->new( -file => $file, -format => 'genbank' );
 while ( my $seq_obj = $gb_obj->next_seq() ) {
  my ( $exon_count, $cds_strand );
  my $id  = $seq_obj->id();
  my $seq = $seq_obj->seq();
  foreach my $feat ( $seq_obj->get_SeqFeatures() ) {
   if (    $feat->location->isa('Bio::Location::SplitLocationI')
        && $feat->primary_tag eq 'CDS' )
   {
    $cds_strand = $feat->strand;
    foreach my $loc ( $feat->location->sub_Location ) {
     $exon_count++;
     push( @{ $hash{$id}{'exons'}{$exon_count} }, ( $loc->start, $loc->end ) );
    }
   }
  }
  $hash{$id}{'seq'}    = $seq;
  $hash{$id}{'strand'} = $cds_strand;
 }

 # have to because geneid demands it in order
 open( FSAOUT, ">$fasta_out" );
 foreach my $id ( sort keys %hash ) {
  print FSAOUT ">$id\n" . $hash{$id}{'seq'} . "\n";
 }
 close FSAOUT;
 Golden::Filter::process_cmd("$main::cdbfasta_exec $fasta_out 2>/dev/null");
 return \%hash;
}

# ---------------------------------------------------------------------------
# gb2gff3: convert the parse_gb() hashref to GFF3, sorting via Golden::Filter.
# Lifted from v1 lines 1814-1863. The trailing `&sort_gff3($out)` is routed
# to Golden::Filter::sort_gff3 per spec.
# ---------------------------------------------------------------------------
sub gb2gff3 {
 my $hash_ref = shift;
 my $out      = shift;
 open( GFF, ">$out" );
 foreach my $id ( sort keys %{$hash_ref} ) {
  next unless $hash_ref->{$id}{'strand'};
  my @sorted_exons =
    sort { $a <=> $b } keys %{ $hash_ref->{$id}{'exons'} };
  if ( $hash_ref->{$id}{'strand'} == 1 ) {
   my $strand = '+';
   my ( $gstart, $e ) = @{ $hash_ref->{$id}{'exons'}->{ $sorted_exons[0] } };
   my ( $s, $gend ) = @{ $hash_ref->{$id}{'exons'}->{ $sorted_exons[-1] } };
   print GFF $id
     . "\tGB\texon\t"
     . $gstart . "\t"
     . $gend
     . "\t.\t$strand\t.\t$id.1" . "\n";

   foreach my $exon_counter (@sorted_exons) {
    my ( $start, $stop ) = @{ $hash_ref->{$id}{'exons'}->{$exon_counter} };
    print GFF $id
      . "\tGB\texon\t"
      . $start . "\t"
      . $stop
      . "\t.\t$strand\t.\t$id.1" . "\n";
   }
  }
  elsif ( $hash_ref->{$id}{'strand'} == -1 ) {
   my $strand = '-';
   my ( $gstart, $e ) = @{ $hash_ref->{$id}{'exons'}->{ $sorted_exons[-1] } };
   my ( $s, $gend ) = @{ $hash_ref->{$id}{'exons'}->{ $sorted_exons[0] } };
   print GFF $id
     . "\tGB\texon\t"
     . $gstart . "\t"
     . $gend
     . "\t.\t$strand\t.\t$id.1" . "\n";
   foreach my $exon_counter (@sorted_exons) {
    my ( $start, $stop ) = @{ $hash_ref->{$id}{'exons'}->{$exon_counter} };
    print GFF $id
      . "\tGB\texon\t"
      . $start . "\t"
      . $stop
      . "\t.\t$strand\t.\t$id.1" . "\n";
   }
  }
  print GFF "\n\n";
 }
 close GFF;
 Golden::Filter::sort_gff3($out);
}

# ---------------------------------------------------------------------------
# gff2hints: emit an Augustus extrinsic-hints file from a sorted GFF3.
# Lifted verbatim from v1 lines 2219-2278.
# ---------------------------------------------------------------------------
sub gff2hints {
 my $gff         = shift;
 my $golden      = shift;
 my $src         = $golden ? 'GLD' : 'XNT';
 my $type_suffix = $golden ? '' : 'part';
 my $priority    = $golden ? 7 : 5;
 my $delimiter   = "\n\n";
 my $orig_sep    = $/;
 open( GFF, $gff )          || die;
 open( OUT, ">$gff.hints" ) || die;
 $/ = $delimiter;

 while ( my $record = <GFF> ) {
  my @lines = split( "\n", $record );

  #my $gene_line = $lines[0];
  my $mRNA_line = $lines[1];
  $mRNA_line =~ /ID=([^;]+)/;
  my $mRNA_id = $1;
  $mRNA_id =~ s/[^\w\.]+/_/g;

  # CDS intron
  for ( my $i = 2 ; $i < scalar(@lines) ; $i++ ) {
   my @data = split( "\t", $lines[$i] );
   next unless $data[8];
	#seems that some programs do not set start stop for inverse strands correctly!
	if ($data[4] < $data[3]){
		my $te = $data[3];
		$data[3] = $data[4];
		$data[4] = $te;
	}
   if ( $data[2] eq 'mRNA' ) {
    $data[2] = 'genicpart';
    print OUT join( "\t", @data[ 0 .. 7 ] )
      . "\tsrc=$src;pri=$priority;grp=$mRNA_id\n";
   }
   elsif (    $data[2] eq 'exon'
           || $data[2] eq 'CDS'
           || $data[2] eq 'intron' )
   {
    $data[2] .= $type_suffix;
    print OUT join( "\t", @data[ 0 .. 7 ] )
      . "\tsrc=$src;pri=$priority;grp=$mRNA_id\n";
   }
   elsif ( $data[2] =~ /splice/ ) {
    $data[2] = $data[2] eq 'three_prime_cis_splice_site' ? 'ass' : 'dss';
    print OUT join( "\t", @data[ 0 .. 7 ] )
      . "\tsrc=$src;pri=$priority;grp=$mRNA_id\n";
   }
   elsif ( $data[2] =~ /UTR/i ) {
    $data[2] = 'UTRpart';
    print OUT join( "\t", @data[ 0 .. 7 ] )
      . "\tsrc=$src;pri=$priority;grp=$mRNA_id\n";
   }
  }
 }
 close GFF;
 close OUT;
 $/ = $orig_sep;
}

1;
