#!/usr/bin/env perl

=pod

=head1 NAME

 extract_cds_on_frameshifts.pl
 
=head1 USAGE

Mandatory:

  input|fasta :s => Inputfile
  mtdna :s       => BLAST database with mitochondrial proteins
  rrna :s        => rDNA gene BLAST database  
  nuclear :s     => BLAST database with proteins (can also including mtDNA if they have been also included in mtDNA above)
  debug          => debug verbosity
  threads :i     => Number of threads
  gencon_mt :i   => Genetic code for mitochondrial (def 5; Insect)
  gencon_nt :i   => Genetic code for nuclear genes (def 1; Standard)

Optional:

 codon_cusp_nuc :s    => Nuclear codon usage file
 codon_cusp_mt  :s    => Mitochondrial codon usage file
 blast_format   :s    => BLAST format 
 blast_nuclear  :s    => A BLASTX for nuclear genes
 blast_rrna     :s    => A BLASTN for rRNA genes
 blast_mito     :s    => A BLASTX for mtDNA genes

=head1 AUTHORS

 Alexie Papanicolaou

        1 CSIRO Ecosystem Sciences, GPO 1700, Canberra 2601, Australia
        alexie@butterflybase.org
        
 Many thanks/credit to James Wasmuth (Toronto) for code bundled in prot4EST

=head1 DISCLAIMER & LICENSE

This software is released under the GNU General Public License version 3 (GPLv3).
It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html.
Please note that incorporating the whole software or parts of its code in proprietary software
is prohibited under the current license.

=head1 BUGS & LIMITATIONS

 None known. Probably lots

 BLAST with normal output can be erroneous in NFS systems. Use XML!
 BLAST engine error: Failed to parse sequence range (range cannot be empty)
 
 
=cut

use strict;
use warnings;
use Data::Dumper;
use Time::localtime;
use Bio::SearchIO;
use Getopt::Long;
use Pod::Usage;
use Tie::IxHash;
use Algorithm::Loops qw( NestedLoops );
use Statistics::Descriptive;
use FindBin qw/$RealBin/;
$|=1;

my ( $debug, $infile, $mitochondrial_db, $rrna_db, $nuclear_db, $nuc_cusp_file,
     $mt_cusp_file );
my ( $mito_blast, $rrna_blast, $nuclear_blast );
my $threads        = 4;
my $mtDNA_gencon   = 5;         # insects
my $nuclear_gencon = 1;         #std
my $blast_user_format   = 'blast';
my $nuc_cutoff = '1e-10';
my $mt_cutoff = '1e-20'; 
GetOptions(
 'input|fasta=s'    => \$infile,
 'mtdna:s'          => \$mitochondrial_db,
 'rrna:s'           => \$rrna_db,
 'nuclear=s'        => \$nuclear_db,
 'debug'            => \$debug,
 'threads:i'        => \$threads,
 'gencon_mt:i'      => \$mtDNA_gencon,
 'gencon_nuc:i'     => \$nuclear_gencon,
 'codon_cusp_nuc:s' => \$nuc_cusp_file,
 'codon_cusp_mt:s'  => \$mt_cusp_file,
 'blast_format:s'   => \$blast_user_format,
 'blast_nuclear:s'  => \$nuclear_blast,
 'blast_rrna:s'     => \$rrna_blast,
 'blast_mito:s'     => \$mito_blast,
 'cutoff_nuc:s'     =>\$nuc_cutoff,
 'cutoff_mt:s'     =>\$mt_cutoff,
);

pod2usage "No input!\n" unless $infile && -s $infile && $nuclear_db;

pod2usage "Please provide genetic codes!\n"
  unless $mtDNA_gencon && $nuclear_gencon;

my $cwd = `pwd`;
chomp($cwd);
$blast_user_format = lc($blast_user_format);
die "BLAST format can only be blastxml or blast\n" unless ($blast_user_format eq 'blast' || $blast_user_format eq 'blastxml');
my $blastx_exec      = &get_prog('blastx');
my $blastn_exec      = &get_prog('blastn');
my $transeq_exec     = &get_prog('transeq');
my $blastdbcmd_exec  = &get_prog('blastdbcmd');
my $makeblastdb_exec = &get_prog('makeblastdb');
my $cd_hit_est_exec  = &get_prog('cd-hit-est');
my $gencode_ncbi     = $RealBin . '/lib/gencode.dmp';

my %one2Three = (
                  'A' => 'Ala',
                  'C' => 'Cys',
                  'D' => 'Asp',
                  'E' => 'Glu',
                  'F' => 'Phe',
                  'G' => 'Gly',
                  'H' => 'His',
                  'I' => 'Ile',
                  'K' => 'Lys',
                  'L' => 'Leu',
                  'M' => 'Met',
                  'N' => 'Asn',
                  'P' => 'Pro',
                  'Q' => 'Gln',
                  'R' => 'Arg',
                  'S' => 'Ser',
                  'T' => 'Thr',
                  'V' => 'Val',
                  'W' => 'Trp',
                  'Y' => 'Tyr',
                  '*' => 'End'
);

$mito_blast = $mito_blast ? $mito_blast : $infile . "_vs_" . 'mitochondrial';
$rrna_blast = $rrna_blast ? $rrna_blast : $infile . "_vs_" . 'rrna';
$nuclear_blast = $nuclear_blast ? $nuclear_blast : $infile . "_vs_" . 'nuclear';

my ($nucSeqIDs,$seq_hash_ref) = &prepare_ids($infile);
foreach my $v ( \( values %{$nucSeqIDs} ) ) {
 $$v = 'none';
}
my $number_seqs = scalar( keys %{$nucSeqIDs} );
print &mytime . $number_seqs . " sequences to be translated\n";
die "No sequences!\n" unless $number_seqs;

# do blasts
print "If not all BLASTs provided, doing BLASTs...\n";
&do_blast( $infile, $blast_user_format, $mitochondrial_db, 'mitochondrial', $mito_blast )
  if $mitochondrial_db && !-s $mito_blast;
&do_blast( $infile, $blast_user_format, $rrna_db, 'rrna', $rrna_blast )
  if $rrna_db && !-s $rrna_blast;
&do_blast( $infile, $blast_user_format, $nuclear_db, 'nuclear', $nuclear_blast )
  if !-s $nuclear_blast;

print "Continuing...\n";
unless ( $nuc_cusp_file && -s $nuc_cusp_file ) {

 # reduce redundancy; long ORFs only
 &process_cmd(
"$cd_hit_est_exec -r 1 -i $infile -o $infile.1000bp.nr90 -M 0 -T $threads -d 0 -l 1000 >/dev/null"
 ) unless -s "$infile.1000bp.nr90";
 die "Clustering failed.\n" unless -s "$infile.1000bp.nr90";

 #NUCLEAR cusp file; biased for conserved proteins:
 $nuc_cusp_file =
   &blast2codon( "$infile.1000bp.nr90", $nuclear_blast, '1e-30', 'nuclear' );
 die
"Nuclear gene BLASTX file $nuclear_blast does not have enough significant hits to produce a nuclear gene codon usage table...\n"
   unless $nuc_cusp_file && -s $nuc_cusp_file;
}

unless ( $mt_cusp_file && -s $mt_cusp_file ) {

 # reduce redundancy; long ORFs only
 &process_cmd(
"$cd_hit_est_exec -r 1 -i $infile -o $infile.300bp.nr90 -M 0 -T $threads -d 0 -l 300 >/dev/null"
 ) unless -s "$infile.300bp.nr90";

 die "Clustering failed.\n" unless -s "$infile.300bp.nr90";

 #mtDNA cusp file; biased for conserved proteins:
 $mt_cusp_file =
   &blast2codon( "$infile.300bp.nr90", $mito_blast, '1e-30', 'mitochondrial' );
 unless ( $mt_cusp_file && -s $mt_cusp_file ) {
  warn
"mtDNA gene BLASTX file $mito_blast does not have enough significant hits to produce a mtDNA gene codon usage table. Will use the nuclear one (probably will not matter as you don't have many (if any) mtDNA genes)\n";
  $mt_cusp_file = $nuc_cusp_file;
 }
}

print "Translating...\n";
my $nt_hits   = int(0);
my $mito_hits   = int(0);
my $rrna_hits = int(0);

# rRNA
my $rna_outFile = $infile . '_protpred_rRNA_nuc.fsa';
if ( !-f $rna_outFile.'.completed' ) {
 my $rrna_parse = &parse_blast( $rrna_blast, $blast_user_format,1, '1e-65' );
 $rrna_hits = scalar keys %{$rrna_parse};
 if ( $rrna_hits && $rrna_hits > 0 ) {
  open( OUT, ">$rna_outFile.dummy" );
  foreach my $seqID ( keys %{$rrna_parse} ) {
   $nucSeqIDs->{$seqID} = 'rRNA';
   print OUT $seqID . "\n";
  }
  close OUT;
  &process_cmd("$blastdbcmd_exec -db $infile  -outfmt \%f -entry_batch $rna_outFile.dummy > $rna_outFile") if -s "$rna_outFile.dummy";
  unlink( $rna_outFile . '.dummy' );
  print &mytime . "Kept as rRNA genes: $rrna_hits\n";
  &touch_empty($rna_outFile.'.completed');
 }
 else {
  print &mytime . "Found no rRNA genes\n";
  &touch_empty($rna_outFile.'.completed');
 }
}else {
 my @hits = `grep '^>' $rna_outFile`;
 foreach my $hit (@hits){
  $hit=~/^>(\S+)/;
  if ($1){
   $nucSeqIDs->{$1} = 'rRNA';
   $rrna_hits++;
  }
 }
 print "Existing output with $rrna_hits rRNA genes found.\n";
}

#MTDNA
my $mito_outFile = $infile . '_protpred_mtDNA_nuc.fsa';
my $mito_transFile = $infile . '_protpred_mtDNA_pro.fsa';
if ( !-f $mito_transFile.'.completed' ) {
 $mito_hits =
   &parse_tiled_blastx_ala_prot4est(
                                     $mito_blast,    $mt_cutoff,
                                     'mitochondrial', $mt_cusp_file,
                                     $mito_outFile, $mito_transFile
   );
 if ( $mito_hits ) {
  &touch_empty($mito_transFile.'.completed');
  print &mytime . "Translated as mtDNA genes: $mito_hits\n";
 }
}else {
 my @hits = `grep '^>' $mito_outFile`;
 foreach my $hit (@hits){
  $hit=~/^>(\S+)/;
  if ($1){
   $nucSeqIDs->{$1} = 'mitochondrial';
   $mito_hits++;
  }
 }
 print "Existing output with $mito_hits mtDNA genes found.\n";
}

#NUC
my $nuclear_outFile = $infile . '_protpred_nuclear_nuc.fsa';
my $nuclear_transFile = $infile . '_protpred_nuclear_pro.fsa';
if ( !-f $nuclear_transFile.'.completed' ) {
 $nt_hits =
   &parse_tiled_blastx_ala_prot4est( $nuclear_blast, $nuc_cutoff, 'nuclear',
                                     $nuc_cusp_file, $nuclear_outFile,$nuclear_transFile );
 if ( $nt_hits ) {
  &touch_empty($nuclear_transFile.'.completed');
  print &mytime . "Translated as nuclear genes: $nt_hits\n";
 }
}
else {
 my @hits = `grep '^>' $nuclear_outFile`;
 foreach my $hit (@hits){
  $hit=~/^>(\S+)/;
  if ($1){
   $nucSeqIDs->{$1} = 'nuclear';
   $nt_hits++;
  }
 }
 print "Existing output with $nt_hits nuclear genes found.\n";
}

my $total = $nt_hits + $mito_hits + $rrna_hits;
print "A total of $total out of $number_seqs has been translated ("
  . sprintf( "%.2f %%", ( $total / $number_seqs ) * 100 )
  . ")\n";

if (!-s "$infile.untranslated.fsa" && $total < scalar( keys %$nucSeqIDs)){  
 print "Printing out untranslated: $infile.untranslated.fsa\n" ;
  open( FILTOUT, ">$infile.untranslated" ) || die("Cannot write to $infile.untranslated\n".$!);
  foreach my $id ( keys %$nucSeqIDs ) {
   print FILTOUT "$id\n" if $nucSeqIDs->{$id} eq 'none';
  }
  close FILTOUT;
  &process_cmd( "$blastdbcmd_exec -db $infile  -outfmt \%f -entry_batch $infile.untranslated > $infile.untranslated.fsa"  ) if -s "$infile.untranslated";
  unlink("$infile.untranslated");
}
#################################
sub get_prog() {
 my $prog = shift;
 my $path = `which $prog`;
 die "Program $prog not found\n" if ( !$path );
 chomp($path);
 print "Using $path\n" if $debug;
 return $path;
}

sub touch_empty(){
  my $file = shift ||return;
  system("date > $file")
}

sub process_cmd {
 my ( $cmd, $dir ) = @_;
 print "CMD: $cmd\n" if $debug;
 chdir($dir) if $dir;
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  die "Error, cmd died with ret $ret\n";
 }
 chdir($cwd) if $dir;
 return;
}

sub mytime() {
 my $hour =
   localtime->hour() < 10 ? '0' . localtime->hour() : localtime->hour();
 my $min = localtime->min() < 10 ? '0' . localtime->min() : localtime->min();
 my $sec = localtime->sec() < 10 ? '0' . localtime->sec() : localtime->sec();
 return "$hour:$min:$sec\t";
}

sub do_blast() {
 my ( $input, $blast_format, $db, $type, $outfile ) = @_;
 my $exec =
     $type eq 'rrna'
   ? $blastn_exec
   : $blastx_exec;
 return $outfile if -s $outfile;
 my $cmd = $exec
   . " -num_descriptions 10 -num_alignments 10 -db $db -query $input -num_threads $threads -out $outfile ";
 $cmd .= " -evalue 1e-20 -query_gencode $mtDNA_gencon -soft_masking true "
   if $type eq 'mitochondrial';
 $cmd .= $blast_format eq 'blast' ? ' -outfmt 0 ' : ' -outfmt 5 ' ;
 $cmd .= " -task megablast -evalue 1e-40 -dust no " if $type eq 'rrna';
 $cmd .= " -evalue 1e-10 -query_gencode $nuclear_gencon -soft_masking true " if $type eq 'nuclear';
 &process_cmd($cmd);
 return $outfile;
}

sub parse_blast() {

 # adapted from James Wasmuth's prot4est
 my ( $blsFile,$blast_format, $hspLimit, $eCut, $allowed_ids ) = @_;
 $eCut = 1e-1 unless $eCut;
 my $blast_check = `head -n 1 $blsFile`;
 die "I don't think $blsFile is of blast format $blast_format" unless (($blast_check=~/^T?BLAST/ && $blast_format eq 'blast') || ($blast_check=~/xml/ && $blast_format eq 'blastxml'));
 my $blsIO = Bio::SearchIO->new( -file => $blsFile, -format => $blast_format, -signif => $eCut );
 my %hspData;
 $allowed_ids = $nucSeqIDs if !$allowed_ids;
RES: while ( my $res = $blsIO->next_result ) {
  next unless $res->num_hits() > 0;
  my $qName = $res->query_name;
  # blastxml issue:
  if ($qName=~/^Query_\d+/){
   $qName = $res->query_description();
   $qName =~/^(\S+)/;
   $qName=$1;
  } 
  next if !$allowed_ids->{$qName};
  my $hspCount = 0;
HIT: while ( my $hit = $res->next_hit ) {
   next HIT unless ( $hit->significance <= $eCut );
   my $hName = $hit->name;
   my @hsps;
 HSP: while ( my $hsp = $hit->next_hsp ) {
    next HSP unless ( $hsp->evalue <= $eCut );

    my $qStart = $hsp->start('query');
    my $qEnd   = $hsp->end('query');
    push @{ $hspData{$qName} },
      {
        Evalue  => $hsp->evalue,
        Score   => $hsp->score,
        QStart  => $qStart,
        QEnd    => $qEnd,
        QLength => $hsp->length('query'),
        QStrand => $hsp->strand,
        HStart  => $hsp->start('hit'),
        HEnd    => $hsp->end('hit'),
        HName   => $hName
      };

    if ( ++$hspCount >= $hspLimit ) {
     next RES;
    }
   }
  }
 }
 return \%hspData;
}

sub blast2codon() {

 my ( $allowed_fasta, $blastfile, $evalue, $type ) = @_;
 my (%allowed_ids);
 my $codonfile = $infile . '.' . $type . '_codontable';
 print "Using BLAST (eval < $evalue) to create a codon usage table for $type genes\n";
 my ($gcTT,$stop_codon_table) = $type eq 'nuclear' ? &getTransTable( $gencode_ncbi, $nuclear_gencon ) :  &getTransTable( $gencode_ncbi, $mtDNA_gencon );
 my @ids = `grep '^>' $allowed_fasta`;
 foreach my $id (@ids) {
  $id =~ /^>(\S+)/;
  $allowed_ids{$1} = 1 if $1;
 }
 my $seqnumber      = scalar( keys %allowed_ids );
 my $blsData        = &parse_blast( $blastfile,$blast_user_format, 10, $evalue, \%allowed_ids );
 my $seqnumber_hits = scalar( keys %{$blsData} );
 
 if (!$seqnumber_hits || $seqnumber_hits == 0){
  warn "No hits found for $type. Skipping this type of sequences\n";
  return;
 }
 
 warn "You have few sequences with significant BLASTX hits ($seqnumber_hits out of $seqnumber) for creating the codon usage table. Careful...\n"
   if $seqnumber_hits < 5 && $type eq 'mitochondrial';
 warn "You have very few sequences with significant BLASTX hits ($seqnumber_hits out of $seqnumber) for creating the codon usage table. Careful...\n"
   if $seqnumber_hits < 100 && $type eq 'nuclear';
 my $codonUsage = &codonUsageFromBLAST( $blsData, $infile, $gcTT );
 my ( $codonTable, $totalCodons ) = &calculateFreqTable( $codonUsage, $gcTT );
 &writeCodFile( $codonTable, $codonfile . '.cod' );
 &gcg2cusp( $codonfile . '.cod', $codonfile . '.cusp' );
 return $codonfile . '.cusp';
}

sub prepare_seq(){
 my $file = shift;
 my %hash;
 open (my $fh,$file);
 my ($seq,$id);
 while (my $ln=<$fh>){
  if ($ln=~/^>(\S+)/){
   $hash{$id} = &iupac_replace(uc($seq)) if $id && $seq;
   $id = $1;
   $seq = '';
  }else{
   chomp($ln);
   $seq .=$ln;
  }
 }
 $hash{$id} = &iupac_replace(uc($seq)) if $id && $seq;
 close $fh;
 return \%hash;
}

sub prepare_ids() {
 print "Preparing and counting seqs...\n";
 my %ids;
 &process_cmd(
          "$makeblastdb_exec -in $infile -dbtype nucl -parse_seqids >/dev/null")
   unless -s "$infile.nhr";
   my $cmd = "$blastdbcmd_exec -db $infile -outfmt %a -entry all";
 print $cmd."\n" if $debug;
 my @ids = `$cmd`;
 foreach my $id (@ids) {
  chomp($id);
  $ids{$id}++;
  die "Sequence $id exists twice!\n" if $ids{$id} > 1;
 }
 my $seq_hash_ref = &prepare_seq($infile);
 return (\%ids,$seq_hash_ref);
}

sub parse_tiled_blastx_ala_prot4est() {

 # adapted from James Wasmuth's prot4est
 my ( $blastfile, $evalue, $type, $cusp_file_used, $nuc_outfile, $protein_outfile ) = @_;
 print "Parsing $blastfile...\n";
 my $hspLimit         = 10;       #for tiling
 my $significant_hits = int(0);
 my $blsData = &parse_blast( $blastfile,$blast_user_format, $hspLimit, $evalue );
 my $codon_table = $type eq 'nuclear' ? $nuclear_gencon : $mtDNA_gencon;
 my $gcTT = &getTransTable( $gencode_ncbi, $codon_table);
 # Parse through each entry and ignore sequences already dealt with
 # then create a tile path through them
 my $gapLimit = 10;
 open (NUCOUT, ">$nuc_outfile") || die;
 open (PROTOUT, ">$protein_outfile") ||die;
 my $counter = int(0);
QUERY: foreach my $qID ( keys %$blsData ) {
  next QUERY unless ( exists $nucSeqIDs->{$qID} );
  $counter++;
  print "$counter/$number_seqs sequences processed\r" if ($counter =~/000$/);
  next QUERY unless $nucSeqIDs->{$qID} eq 'none';
  my $data = $blsData->{$qID};
  my $tileData =
    &buildTilePath( $qID, $data, $gapLimit, $infile, $cusp_file_used );
  $nucSeqIDs->{$qID} = $type;
  my @seq;
  $nucSeqIDs->{$qID} = $type;
  foreach my $tile (@$tileData) {
   if ( my $insert = $tile->{QInsert} ) {
    push (@seq, $insert);
   }
   else {
    my $nucSeq = &get_indexed_seq( $qID, $tile->{QStart}, $tile->{QEnd}, $tile->{QStrand} );
    push (@seq, $nucSeq);
   }
  }
  my $desc = ' len:'.length(join( "", @seq ));
  # create sequence in the correct frame
  if ( $tileData->[-1]->{QStrand} == -1 ) {
   @seq = reverse @seq;
   $desc .= " (-) $qID:".$tileData->[-1]->{QStart}.'-'.$tileData->[0]->{QEnd}." (-)";
  }else {
   $desc .= " (+) $qID:".$tileData->[0]->{QStart}.'-'.$tileData->[-1]->{QEnd}." (+)";
  }
  
  my ($has_stop_codon);
  my $putative_stop_codon = $tileData->[-1]->{QStrand} == 1 ? &get_indexed_seq( $qID, ($tileData->[-1]->{QEnd})+1, ($tileData->[-1]->{QEnd})+3, $tileData->[-1]->{QStrand} ) 
    : &get_indexed_seq( $qID, ($tileData->[-1]->{QStart})-3, ($tileData->[-1]->{QStart})-1, $tileData->[-1]->{QStrand} );
    if ($putative_stop_codon){
     if ($putative_stop_codon && length($putative_stop_codon)==3 && &is_stop_codon($putative_stop_codon,$type)){ 
      $has_stop_codon = 1;
      push (@seq, $putative_stop_codon) ;
    }
  }
  my $seq = join( "", @seq );
  my $orf_part='';
  # check for pseudogene (slows it down )
  my $protein_seq = &translate($seq,$codon_table);
  my $stop_codon_number = ( $protein_seq =~ tr/\*// ); 
  #print "$qID Has $stop_codon_number stopcodons\n" if $debug ;
  if ($stop_codon_number > 1 || ($stop_codon_number == 1 && !$has_stop_codon)){
    $orf_part = ' type:pseudogene';
  }elsif ($seq=~/^ATG/){
   if ($has_stop_codon){
    $orf_part = ' type:complete';
   }else{
    $orf_part = ' type:3prime_partial';
   }
  }elsif($has_stop_codon){
   $orf_part = ' type:5prime_partial';
  }else{
   $orf_part = ' type:internal';
  }
  
  $desc .= " localization:$type method:similarity";
  print NUCOUT ">$qID CDS $orf_part $desc\n";
  print PROTOUT ">$qID polypeptide $orf_part $desc\n";
  
  while ( $seq =~ m/(\S{1,60})/g ) {
   print NUCOUT $1, "\n";
  }
  
  while ( $protein_seq =~ m/(\S{1,60})/g ) {
   print PROTOUT $1, "\n";
  }
  $nucSeqIDs->{$qID} = $type;
  $significant_hits++;
 }
 close NUCOUT;
 close PROTOUT;
 return int($significant_hits);
}

sub get_indexed_seq {
 my ( $seqID, $start, $end, $strand ) = @_;
 $strand = 1 if !$strand; 
 my $size = abs($end - $start)+1;
 my $seq; 
 my $whole_seq = $seq_hash_ref->{$seqID};
 return if !$whole_seq;
 if ($strand == -1){
  $whole_seq = &revcompl($whole_seq);
  my $start = length($whole_seq) - $end+1;
  $seq =  substr($whole_seq,($start-1),$size);
 }else{
  $seq =  substr($whole_seq,($start-1),$size);
 } 
 return $seq;
}

sub iupac_replace($) {
    my $string = shift;
    $string =~ s/\s+//g;
    $string =~ s/S/G/g;
    $string =~ s/R/G/g;
    $string =~ s/Y/T/g;
    $string =~ s/M/C/g;
    $string =~ s/K/T/g;
    $string =~ s/W/T/g;
    $string =~ s/B/T/g;
    $string =~ s/D/T/g;
    $string =~ s/H/T/g;
    $string =~ s/V/G/g;
    $string =~ s/X/N/g;
    return $string;
}

sub revcompl($) {
  my $seq=shift;
  my $rev = reverse($seq); #reverse
  $rev =~ tr/ACGTacgt/TGCAtgca/; #complement
  return $rev;
}

sub is_stop_codon(){
 my ($codon,$type) = @_;
 my ($gcTT,$stop_codon_table) = $type eq 'nuclear' ? &getTransTable( $gencode_ncbi, $nuclear_gencon ) :  &getTransTable( $gencode_ncbi, $mtDNA_gencon );
 return 1 if $stop_codon_table->{$codon};
}

sub translate(){
   my ($seq,$codon_table,$frame)=@_;
   $codon_table = 1 unless $codon_table;
   $frame = 1 unless $frame;
   #embos 6.4 requires sequence ID
   my $cmd = "echo -e '>id\n$seq'| $transeq_exec -sequence stdin -outseq stdout -supper1 -table $codon_table -frame $frame 2>/dev/null|tail -n +2";
   #print $cmd."\n" if $debug;
   my @protein = `$cmd`;
   my $protein_seq = join('',@protein);
   $protein_seq=~s/\s+//g;
   return $protein_seq;
}

# everything below here is from prot4EST; mostly unadulterated


sub buildTilePath {

 my ( $qID, $hspData, $gapLimit, $nucFile, $cuspFile ) = @_;

 # $hspData is a reference to an array of hashes
 # keys that are present in each hash are:
 #					Evalue
 #					Score
 #					QStart
 #					QEnd
 #					QLength
 #					QStrand
 #					HStart
 #					HEnd
 #					HName

 # All hits should be on the same strand of the nucleotide.
 # The 'correct' frame will be that of the highest
 # scoring (evalue) hit sequence. This will also give the
 # first protein to scan.

 my @sortedHSPs = sort { $a->{Evalue} <=> $b->{Evalue} } @$hspData;

 # Get the order of the protein hits wrt their BLAST score

 my %protOrder;
 tie %protOrder, "Tie::IxHash";
 my $bestStrand;
 foreach my $hsp (@sortedHSPs) {
  $bestStrand = $hsp->{QStrand} unless $bestStrand;
  push @{ $protOrder{ $hsp->{HName} } }, $hsp;
 }

 # Cycle through the matches. Check to see if they overlap a current
 # match. If they do then accept the segment that is already present.
 # Once these are done, then need to see if any gaps between sections
 # are longer than the gaplimit (nucleotides)

 # Store the terminus nucleotide coords of the growing tile

 # Piece together the tile path for each protein hit
 # then join with previous ones

 my @tile;

PROT: while ( my ( $protID, $hsps ) = each %protOrder ) {

  # Want to consider the HSPs in order from 5'-->3'
  # Does this work for the minus strand sequences?
  # Currently they are ordered by bit score.
  # The match with the highest bit score will give
  # the frame to be translated in.
  #print $protID,"\n";
  my @seq;

  # Put the HSPs for this query-protein match in
  # a linear order
  my %orderedHsps;
  my $i = 0;

  my @orderedHsps = sort { $b->{Score} <=> $a->{Score} } @$hsps;

HSP: for ( my $i = 0 ; $i <= $#orderedHsps ; $i++ ) {
   my $thisHsp = $orderedHsps[$i];

   #print $thisHsp->{QStart},"\n";
   next unless ( $thisHsp->{QStrand} == $bestStrand );
   if ( scalar @tile > 0 ) {

    my ( $thisQStart, $thisQEnd, $thisHStart, $thisHEnd, $thisScore,
         $thisHName ) = @{$thisHsp}{qw/QStart QEnd HStart HEnd Score HName/};

    # Is the new HSP 5' or 3' of the tile?
    my $tileStart = $tile[0]->{QStart};
    my $tileEnd   = $tile[-1]->{QEnd};

    if ( $tileStart > $thisQStart ) {

     # The new tile is 5' of the current tile
     my $prevHsp = $tile[0];
     my (
          $prevQStart, $prevQEnd,  $prevHStart,
          $prevHEnd,   $prevScore, $prevHName
     ) = @{$prevHsp}{qw/QStart QEnd HStart HEnd Score HName/};

     # Next question is does it overlap or is there a gap?
   OVERLAP: if ( $prevQStart <= $thisQEnd ) {    # Yes they do overlap

      # One of the HSPs will need to be hacked back.
      # Assumptions:
      # i) The higher scoring (bits) is the correct one
      # ii) In hacking back the lower sequence it is likely
      # that nucleotides will be lost.

      if ( $prevScore > $thisScore ) {

       # Remove from $thisHsp

       while ( $prevQStart <= ( $thisQEnd -= 3 ) ) { }

       if ( $thisQEnd > $thisQStart ) {
        $thisHsp->{QEnd} = $thisQEnd;
       }
      }
      else {

       # Edit the previous HSP
       while ( $thisQEnd >= ( $prevQStart += 3 ) ) { }

       # This may role back the entire hsp
       if ( $prevQStart >= $prevQEnd ) {
        shift @tile;
       }
       else {
        $tile[0]->{QStart} = $prevQStart;
       }
      }
      unshift @tile, $thisHsp;
      $i = 0;
     }

     else {

      # There's a gap
      my $Ngap = $prevQStart - $thisQEnd - 1;    # Size of the gap
      next HSP if $Ngap > $gapLimit;
      my $Pgap = 0;
      if ( $thisHName eq $prevHName ) {
       if ( $bestStrand > 0 ) {
        $Pgap = $prevHStart - $thisHEnd - 1;
       }
       else {
        $Pgap = $thisHStart - $prevHEnd - 1;
       }
      }

      # If Pgap==0 and Ngap<=2 then most probably an incorrect insertion.
      # So ignore the nucleotides
      if ( $Pgap == 0 && $Ngap <= 2 ) {
       unshift @tile, $thisHsp;
      }
      elsif ( $Pgap == 0 && $Ngap >= 3 ) {

       # If $Ngap%3==0 simply add in co-ordinates into the tile
       # Else need to use the cusp table to find the most likely codon(s)
       my $mod = $Ngap % 3;
       my %newHSP;
       if ( $mod == 0 ) {
        %newHSP = (
                    QStart  => $thisQEnd + 1,
                    QEnd    => $prevQStart - 1,
                    QStrand => $bestStrand
        );
       }
       else {
        my $corrected = &likelyCodonCusp( $qID,
                                       $thisQEnd + 1,
                                       $prevQStart - 1,
                                       $bestStrand, $mod, $nucFile, $cuspFile );
        next PROT if !$corrected;
        %newHSP = ( QInsert => $corrected );
       }
       unshift @tile, $thisHsp, \%newHSP;
      }

      # If Pgap>1 see what Ngap is.
      # If Ngap is divisible by three then accept wholesale
      # Else will need to pick the codons either:
      # i) based on their predicted frequency OR
      # ii) by comparison to the hit sequence
      elsif ( $Pgap && $Ngap == 0 ) {
       unshift @tile, $thisHsp;
      }
      elsif ( $Pgap && $Ngap % 3 == 0 ) {
       my %newHSP = (
                      QStart  => $thisQEnd + 1,
                      QEnd    => $prevQStart - 1,
                      QStrand => $bestStrand
       );
       unshift @tile, $thisHsp, \%newHSP;
      }

      elsif ( $Pgap && $Ngap % 3 > 0 ) {
       my $Pngap = $Pgap * 3;    # Number of nucleotides

       # This is where either the best scoring codons are accepted,
       # or the BLAST match itself can be used.
       # The former will be used as default and the latter coded,
       # at a later point.

       my $mod       = $Ngap % 3;
       my $corrected = &likelyCodonCusp( $qID,
                                       $thisQEnd + 1,
                                       $prevQStart - 1,
                                       $bestStrand, $mod, $nucFile, $cuspFile );
       next PROT if !$corrected;
       my %newHSP = ( QInsert => $corrected );
       unshift @tile, $thisHsp, \%newHSP;
      }
      else {
       print "\n#Internal Error: - Cannot figure out where the tile fits.\n",
         "Danger of continuous looping - $qID => $protID\n",
         "\n";
       exit;
      }
      $i = 0;
     }
     next HSP;
    }

    elsif ( $thisQEnd > $tileEnd ) {

     # The new sequence is 3' of the current tile
     my $prevHsp = $tile[-1];
     my (
          $prevQStart, $prevQEnd,  $prevHStart,
          $prevHEnd,   $prevScore, $prevHName
     ) = @{$prevHsp}{qw/QStart QEnd HStart HEnd Score HName/};

     # Next question is does it overlap or is there a gap?
   OVERLAP: if ( $prevQEnd >= $thisQStart ) {    # Yes they do overlap
          # One of the HSPs will need to be hacked back.
          # Assumptions:
          # i) The higher scoring (bits) is the correct one
          # ii) In hacking back the lower sequence it is likely
          # that nucleotides will be lost.
      if ( $prevScore >= $thisScore ) {

       # Remove from $thisHsp
       while ( $prevQEnd >= ( $thisQStart += 3 ) ) { }
       if ( $thisQStart < $thisQEnd ) {
        $thisHsp->{QStart} = $thisQStart;
       }
      }
      else {

       # Edit the previous HSP
       while ( $thisQStart <= ( $prevQEnd -= 3 ) ) { }
       if ( $prevQStart >= $prevQEnd ) {
        pop @tile;
       }
       else {
        $tile[-1]->{QEnd} = $prevQEnd;
       }
      }
      push @tile, $thisHsp;
      $i = 0;
     }
     else {    # There is a nice gap
      my $Ngap = $thisQStart - $prevQEnd - 1;    # Size of the gap
      next HSP if $Ngap > $gapLimit;

      my $Pgap = 0;
      if ( $thisHName eq $prevHName ) {
       if ( $bestStrand > 0 ) {
        $Pgap = $thisHStart - $prevHEnd - 1;    # Size of the gap on the protein
       }
       else {
        $Pgap = $prevHStart - $thisHEnd - 1;
       }
      }

      # If Pgap==0 and Ngap<=2 then most probably an incorrect insertion.
      # So ignore the nucleotides
      if ( $Pgap == 0 && $Ngap <= 2 ) {
       push @tile, $thisHsp;
      }

      # Elsif Pgap==0 and Ngap>=3 then we're looking at a true insertion
      elsif ( $Pgap == 0 && $Ngap >= 3 ) {

       # If $Ngap%3==0 simply add in co-ordinates into the tile
       # Else need to use the cusp table to find the most likely codon(s)
       my $mod = $Ngap % 3;
       my %newHSP;
       if ( $mod == 0 ) {
        %newHSP = (
                    QStart  => $prevQEnd + 1,
                    QEnd    => $thisQStart - 1,
                    QStrand => $bestStrand
        );
       }
       else {
        my $corrected = &likelyCodonCusp( $qID,
                                       $prevQEnd + 1,
                                       $thisQStart - 1,
                                       $bestStrand, $mod, $nucFile, $cuspFile );
        next PROT if !$corrected;
        %newHSP = ( QInsert => $corrected );
       }
       push @tile, \%newHSP, $thisHsp;
      }

      # If Pgap>1 see what Ngap is.
      # If Ngap is divisible by three then accept wholesale
      # Else will need to pick the codons either:
      # i) based on their predicted frequency OR
      # ii) by comparison to the hit sequence
      elsif ( $Pgap && $Ngap == 0 ) {
       push @tile, $thisHsp;
      }
      elsif ( $Pgap && $Ngap % 3 == 0 ) {
       my %newHSP = (
                      QStart  => $prevQEnd + 1,
                      QEnd    => $thisQStart - 1,
                      QStrand => $bestStrand
       );
       push @tile, \%newHSP, $thisHsp;
      }

      elsif ( $Pgap && $Ngap % 3 > 0 ) {
       my $Pngap = $Pgap * 3;    # Number of nucleotides
            # This is where either the best scoring codons are accepted,
            # or the BLAST match itself can be used.
            # The former will be used as default and the latter coded,
            # at a later point.

       my $mod       = $Ngap % 3;
       my $corrected = &likelyCodonCusp( $qID,
                                       $prevQEnd + 1,
                                       $thisQStart - 1,
                                       $bestStrand, $mod, $nucFile, $cuspFile );
       next PROT if !$corrected;
       my %newHSP = ( QInsert => $corrected );
       push @tile, \%newHSP, $thisHsp;
      }
      else {
       print "\n#Internal Error: - Cannot figure out where the tile fits.\n",
         "Danger of continuous looping - $qID => $protID\n",
         "\n";
       exit;
      }
      $i = 0;
     }
     next HSP;
    }
   }
   else {    # It's the first one
    push @tile, $thisHsp;
   }
  }

  #print scalar @tile,"\n";
  #foreach my $tile (@tile)	{
  #	while (my ($k,$v)=each %$tile)	{
  #		print "$k=>$v\n";
  #	}
  #	print "//\n";
  #}
  #exit;
 }
 return \@tile;
}

sub likelyCodonCusp {
 my ( $qID, $start, $end, $strand, $mod, $nucFile, $cuspFile ) = @_;

 my $codonFreqs  = &getThousand($cuspFile);
 my $wobbleFreqs = &getWobbleThousand($cuspFile);

 #merge the hashes

 my %freqs = ( %$codonFreqs, %$wobbleFreqs );
 my %best;
 $best{score} = 0.1;

 #need to create each posible codon string and score it
 my $nucSeq = &get_indexed_seq( $qID, $start, $end, $strand);
 my $addNs = 3 - $mod;    #need to add this many Ns
 my $minus = $mod;        #or remove this many nucleotides

 my $possibleSeqs = &modifiedSeqs( $nucSeq, $addNs, $minus );
 return if !$possibleSeqs;
 foreach my $seq2score (@$possibleSeqs) {
  my $score = &scoreSeq( $seq2score, \%freqs );

  #print "$seq2score\t$score\n";
  if ( $score > $best{score} ) {
   $best{score} = $score;
   $best{seq}   = $seq2score;
  }
 }
 return $best{seq};
}

sub modifiedSeqs {
 my ( $seq, $addNs, $minusNs ) = @_;
 my @seq = split //, $seq;

 #first add
 my $adds = &addNucs( \@seq, $addNs );
 return if !$adds;
 my $minus = &minusNucs( \@seq, $minusNs );

 my @possible = ( @$adds, @$minus );
 return \@possible;
}

sub addNucs {
 my ( $seq, $N ) = @_;
 my @adds;
 my $i =
   NestedLoops( [ [ 0 .. @$seq ], ( sub { [ $_ .. @$seq ] } ) x ( $N - 1 ) ] );

 while ( my (@positions) = $i->() ) {
  my @tmp = @$seq;
  splice( @tmp, $_, 0, 'N' ) for reverse @positions;
  if ( ( scalar @tmp ) % 3 != 0 ) {
   #warn "\n#Internal Error: - introducing N nucleotides gives incorrect length codons!\n";
   return;

  }
  push @adds, join( '', @tmp );
 }
 return \@adds;
}

sub minusNucs {
 my ( $seq, $N ) = @_;
 my %minus;

 my $i =
   NestedLoops( [ [ 0 .. @$seq ], ( sub { [ $_ .. @$seq ] } ) x ( $N - 1 ) ] );

I: while ( my (@positions) = $i->() ) {
  my @tmp = @$seq;
  for ( reverse @positions ) {
   next I if $_ > $#tmp;
   splice( @tmp, $_, 1 );
  }
  if ( ( scalar @tmp ) % 3 != 0 ) {
   print
"\n#Internal Error: - removing N nucleotides gives incorrect length codons!\n",
     "\n";
   exit;
  }
  $minus{ join( '', @tmp ) }++;
 }
 my @minus = keys %minus;
 return \@minus;
}

sub scoreSeq {
 my ( $seq, $freqs ) = @_;

 my $score = 0;
 if ( length($seq) == 0 ) {
  return $score;
 }
 while ( $seq =~ m/(\w{3})/g ) {
  my $codon = $1;

  if ( $codon eq 'NNN' ) {
   $score += 1;    # nominal score
  }
  elsif ( $codon =~ m/^T(AA|AG|GA)$/ ) {
   $score += -10;    # penalty for stop codon
  }
  else {
   $score += $freqs->{$1} || '0';
  }
 }

 #want the score normalised by the number of codons
 my $numOfCodons = length($seq) / 3;
 unless ( $numOfCodons =~ m/^-?\d+$/ ) {
  print "\n#Internal Error: - the returned string is not divisible by three.\n",
    "Expecting an integer.\n",
    "\n";
  exit;
 }
 $score /= $numOfCodons;
 return $score;
}

sub getTransTable() {
 my ( $gcFile, $genCode ) = @_;
 my $aaString;
 my (%tt,%stop_codons);
 open GC, $gcFile or die;
 while (<GC>) {
  chomp;
  my @entry = split /\t\|\t/;
  if ( $entry[0] == $genCode ) {
   $aaString = $entry[3];
  }
 }
 close GC;

 #The order is qw/T C A G/;
 my @ntides = qw/T C A G/;
 my @aacids = split //, $aaString;
 foreach my $first (@ntides) {
  foreach my $second (@ntides) {
   foreach my $third (@ntides) {
    my $codon = join( '', ( $first, $second, $third ) );
    my $aa = shift @aacids;
    $tt{$codon} = $aa;
    $stop_codons{$codon} = $aa if $aa eq '*';
   }
  }
 }
 return (\%tt,\%stop_codons);
}

sub codonUsageFromBLAST {
 my ( $hspData, $nFile, $gcTT ) = @_;

 #go through each stored HSP and pick the longest one
 #then read through the ntide sequence and count the codons
 #for each amino acid given the translation table

 my %usage;
 while ( my ( $seqID, $hsps ) = each %$hspData ) {

  #find the longest
  my @hsps = sort { $b->{QLength} <=> $a->{QLength} } @$hsps;
  my $hsp = $hsps[0];

  my $start  = $hsp->{QStart};
  my $end    = $hsp->{QEnd};
  my $strand = $hsp->{QStrand};

  my $S = 'plus';    #Strand
  if ( $strand == '-1' ) {
   $S = 'minus';     #reverse complement
  }
  my $L = $start . '-' . $end;
  my $cmd = 
    "$blastdbcmd_exec -entry '$seqID' -db $nFile -range $L -strand $S";
    #print $cmd."\n" if $debug;
  my $seq = `$cmd`;
  
  $seq =~ s/>lcl.*//;
  $seq =~ s/\s+//g;
  foreach my $codon ( unpack "a3" x ( length($seq) / 3 ), $seq ) {
   next unless $codon =~ m/^[AGCT]{3}$/i;
   my $aa = $gcTT->{$codon};
   unless ($aa) {
    print
"\n#Internal Error: %package - Cannot find the amino acid for codon: $codon\n",
      "\n";
   }
   $usage{$aa}->{$codon}++;
  }
 }
 return \%usage;
}

sub calculateFreqTable {
 my ( $codonUsage, $gcTT ) = @_;

 #$cuTable{threeLetterCode}->{codon}->[raw_count,normalised,fraction]
 my %cuTable;

 my $totalCodons = 0;
 my %totalForAA;
 while ( my ( $singleAA, $codons ) = each %$codonUsage ) {
  while ( my ( $codon, $rawCount ) = each %$codons ) {

   #get three letter code
   my $threeAA = $one2Three{$singleAA};
   $cuTable{$threeAA}->{$codon}->[0] = $rawCount;
   $totalCodons += $rawCount;
   $totalForAA{$threeAA} += $rawCount;
  }
 }
 my $coeff = 1000 / $totalCodons;

 #now cycle through and add the normalised number, the fraction
 #and add missing codons

 my %aa2Codon;
 while ( my ( $codon, $aa ) = each %$gcTT ) {
  $aa = $one2Three{$aa};
  $aa2Codon{$aa}->{$codon}++;
 }

 while ( my ( $aa, $codons ) = each %cuTable ) {
  while ( my ( $codon, $nums ) = each %$codons ) {
   $nums->[1] = $nums->[0] * $coeff;
   $nums->[2] = $nums->[0] / $totalForAA{$aa};
   delete $aa2Codon{$aa}->{$codon};
  }

  #are there any codons no represented?
  foreach my $codonLeft ( keys %{ $aa2Codon{$aa} } ) {
   $cuTable{$aa}->{$codonLeft} = [ '0', '0', '0' ];
   delete $aa2Codon{$aa}->{$codonLeft};
  }
 }

 #What if an amino acid has not been included?
 #Unlikely but possible in small datasets

 while ( my ( $aa, $codons ) = each %aa2Codon ) {
  my $numOfCodons = scalar keys %$codons;
  next if $numOfCodons == 0;
  my $frac = ( 1 / $numOfCodons );

  foreach my $codon ( keys %$codons ) {
   $cuTable{$aa}->{$codon} = [ '0', '0', $frac ];
  }
 }

 return ( \%cuTable, $totalCodons );
}

sub writeCodFile {
 my ( $codonTable, $codFile ) = @_;
 open COD, ">", $codFile or die;
 while ( my ( $aa, $codons ) = each %$codonTable ) {
  while ( my ( $codon, $nums ) = each %$codons ) {
   my $outString = sprintf( "%-9s%-9s%15.2f%12.2f%12.2f\n",
                            $aa, $codon,
                            $nums->[0] || 0.00,
                            $nums->[1] || 0.00,
                            $nums->[2] || 0.00 );
   print COD $outString;
  }
 }
 close COD;
 return 1;
}

sub distributionTable {
 my $codonTable = shift @_;
 my %distTable;

 while ( my ( $aa, $codons ) = each %$codonTable ) {
  my $count = '0';
  while ( my ( $codon, $nums ) = each %$codons ) {
   my $freq = $nums->[2];
   $distTable{$aa}->{$codon}->{start} = $count;
   $count += $freq;
   $distTable{$aa}->{$codon}->{end} = $count;

   #if ($aa eq 'Cys')	{
   #	print $distTable{$aa}->{$codon}->{start}," is the start\n";
   #	print $distTable{$aa}->{$codon}->{end}," is the end\n";
   #}
  }
 }
 return \%distTable;
}

sub pickFromDist {
 my ( $distTable, $aa ) = @_;

 #Given a distribution table and an amino acid randomly pick
 #a codon from this distribution

 my $randPick = rand(1);

 my $chosenCodon;

 if ( length($aa) == 1 ) {
  $aa = $one2Three{$aa};
 }

 my $aaDist = $distTable->{$aa};

 while ( my ( $codon, $dist ) = each %$aaDist ) {

  if ( $dist->{start} < $randPick && $dist->{end} > $randPick ) {
   $chosenCodon = $codon;
  }
 }
 unless ($chosenCodon) {
  print "#Internal Error: - No codon picked.\n", "\n";
  exit;
 }
 return $chosenCodon;

}

sub getThousand {
 my $cuspFile = shift;

 my %freq;
 open CUSP, $cuspFile or die;
 while (<CUSP>) {
  chomp;
  next if m/^#/;
  my @entry = split /\s+/;
  unless ( scalar @entry == 5 ) {
   print "\n#Error: $cuspFile does not appear to be in cusp format\n", "\n";
   exit;
  }
  $freq{ $entry[0] } = $entry[3];
 }
 close CUSP;
 return \%freq;
}

sub getWobbleThousand {
 my $cuspFile = shift;

 #read through the frequency table and
 #consider only two bases then take the mean

 my %pairs;
 open CUSP, $cuspFile or die;
 while (<CUSP>) {
  chomp;
  next if m/^#/;
  my @entry = split /\s+/;
  unless ( scalar @entry == 5 ) {
   print "\n#Error: $cuspFile does not appear to be in cusp format\n", "\n";
   exit;
  }
  my $codon = $entry[0];
  my $freq  = $entry[3];

  #build up all possible wobble positions
  my @codon = split //, $codon;

  foreach my $N ( 1, 2 ) {
   my $wobble;
   my $i = NestedLoops(
               [ [ 0 .. @codon ], ( sub { [ $_ .. @codon ] } ) x ( $N - 1 ) ] );

   while ( my (@positions) = $i->() ) {
    my @wob = @codon;
    for my $pos (@positions) {
     $wob[$pos] = 'N';
    }
    next if scalar @wob > scalar @codon;
    $wobble = join( "", @wob );

    #print "$wobble ",$entry[3],"\n";
    unless ( exists $pairs{$wobble} ) {
     $pairs{$wobble} = Statistics::Descriptive::Full->new();
    }
    $pairs{$wobble}->add_data( $entry[3] );
   }
  }
 }
 close CUSP;
 my %freqs;
 while ( my ( $k, $v ) = each %pairs ) {

  #print "$k => ",$v->mean,"\n";
  $freqs{$k} = $v->mean;
 }
 return \%freqs;
}

sub cusp2gcg {
 my ( $cuspFile, $codFile ) = @_;

 open CUSP, $cuspFile or die;
 open COD, ">", $codFile or die;

 while (<CUSP>) {
  next if ( m/^#/ || m/^\s$/ );    #ignore comments
  unless (m/(\w{3})\s+([\w*])\s+([.\d]+)\s+([.\d]+)\s+([.\d]+)\n/) {
   warn "#Internal Error: - The file is not 'cusp' format.\n$_:\n", "\n";
   exit;
  }
  my $codon     = $1;
  my $aa_single = $2;
  my $frac      = $3;
  my $normal    = $4;
  my $num       = $5;

  my $aa_triple = $one2Three{$aa_single};

  my $outstring = sprintf( "%-9s%-9s%15.2f%12.2f%12.2f\n",
                           $aa_triple, $codon,
                           $num    || 0.00,
                           $normal || 0.00,
                           $frac   || 0.00 );

  print COD $outstring;
 }
 close CUSP;
 close COD;
}

sub gcg2cusp {
 my ( $codFile, $cuspFile ) = @_;

 my %three2One = reverse %one2Three;    #reverse the hash

 open COD, $codFile or die;
 open CUSP, ">", $cuspFile or die;

 print CUSP
   "# CUSP codon usage file\n# Codon	Amino acid	Fract   /1000	Number\n";

 while (<COD>) {
  next if ( m/^#/ || m/^\s$/ );         #ignore comments
  unless (m/(\w{3})\s+(\w{3})\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)/) {
   warn "#Internal Error: - The file is not 'gcg' format.\n$_:\n", "\n";
   exit;
  }
  my $aa_triple = $1;
  my $codon     = $2;
  my $num       = $3;
  my $norm      = $4;
  my $frac      = $5;

  my $aa_single = $three2One{$aa_triple};

  my $outstring = sprintf( "%-8s%-9s%12.3f%9.3f%-2s%1.0f\n",
                           $codon, $aa_single,
                           $frac     || 0.00,
                           $norm     || 0.00,
                           ' ', $num || 0 );

  print CUSP $outstring;
 }
 close COD;
 close CUSP;
}

