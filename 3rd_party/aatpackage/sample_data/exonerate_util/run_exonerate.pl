#!/usr/bin/env perl

=pod

=head1 NAME

run_exonerate

=head1 USAGE

  'in|fasta:s'=> Fasta file (cDNA or protein)
  'protein' => File is a protein (otherwise assumed to be cDNA)
  'minorf:i'=> minimum ORF size for cDNA
  'reference:s'=> Fasta with reference sequences (genome)
  'blast_hash:s'=> Blast HASH file from analyze_blast
  'only_coding'=> For cDNA: Use only the coding sequence (coding2genome model). Otherwise (default) use cdna2genome (includes UTR)
  'exclude:s' => Exclude file (one ID per line) or comma separated string of IDs 
  'threads' => Threads to use
  'eval_cutoff:s' => Evalue cutoff for BLAST (def 1e-80)
  'bit_cutoff:i' => Bit score cutoff for BLAST (def 100 for protein tblastn, 500 for nucleotide megablast)
  'prep_only:i' => Set to 1 or more. Do BLAST and prepare commands only. If number >1 is given, then split cmds to that many subfiles (e.g. to run with Parafly)
  'norefine' => Don't refine region in exonerate. Prevents segmentation faults in some searches
  'aat:s' => AAT FILTER output files (one or more). Do not use blast hash
  'score_dps:i' => Minimum DPS score (def 100)
  intron_max => Def 70000
  same_species => FASTA file is from same species as referene genome
 'softmask' => genome has been softmasked with repeatmasker
 'padding' => genome padding for exonerate alignment (2000)
 
 Example cmd with existing AAT filter output:
  run_exonerate.pl -in annotations_merged.pasa_assemblies.denovo_transcript_isoforms.fasta.transdecoder.pep.nr70 -ref csiro4b.10000.assembly_repeatmasked  -soft -thread 4 -prep_only 1 -same -protein -aat *filter -separ
 
 The result will not be run yet. Also the co-ordinates have not been corrected (post-process in another script. See reference ID in GFF output)
 
=cut

# remember to use repeatmasker ~/software/RepeatMasker/RepeatMasker -e abblast -pa 35 -s -species insecta -no_is -gccalc -gff  -xsmall

use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use Bio::Index::Fasta;

my $shuf_exec = `which shuf 2>/dev/null`;    #coreutils
my (
     $fasta_file,       $reference_file, $blast_file, $coding_only,
     $exclude_file,     $is_protein,     $prep_only,  $no_refine,
     @aat_filter_files, $same_species,   $softmasked, $separate
);
my $minorf        = 501;         #minimum orf size in bp
my $threads       = 4;
my $eval_cut      = '1e-80';
my $bit_cut       = 500;
my $dps_min_score = 100;
my $max_intron    = 70000;
my $padding       = 2000;
GetOptions(
            'minorf:i'      => \$minorf,
            'in|fasta:s'    => \$fasta_file,
            'reference:s'   => \$reference_file,
            'blast_hash:s'  => \$blast_file,
            'only_coding'   => \$coding_only,
            'exclude:s'     => \$exclude_file,
            'threads:i'     => \$threads,
            'protein'       => \$is_protein,
            'eval_cutoff:s' => \$eval_cut,
            'bit_cutoff:i'  => \$bit_cut,
            'prep_only:i'   => \$prep_only,
            'norefine'      => \$no_refine,
            'aat:s{,}'      => \@aat_filter_files,
            'score_dps:i'   => \$dps_min_score,
            'intron_max:i'  => \$max_intron,
            'same_species'  => \$same_species,
            'softmask'      => \$softmasked,
            'padding:i'     => \$padding,
            'separate'      => \$separate,
);

pod2usage
  unless $fasta_file && $reference_file && -s $fasta_file && -s $reference_file;
$bit_cut = 100 if $is_protein && $bit_cut == 500;
my $minaa = int( $minorf / 3 );
my %queries;

#my $cmd_processor = `which cmd_process_forker.pl`;chomp ($cmd_processor);
my $cmd_processor = `which ParaFly`;
chomp($cmd_processor);
die "Cannot find ParaFly\n" unless $cmd_processor;
my %exclude_ids;
if ( $exclude_file && -s $exclude_file ) {
 open( F, $exclude_file );
 while (<F>) {
  $_ =~ s/^\s+//;
  chomp($_);
  $exclude_ids{$_} = 1;
 }
 close(F);
}
elsif ( $exclude_file && $exclude_file =~ /,/ ) {
 my @a = split( ',', $exclude_file );
 foreach (@a) {
  $exclude_ids{$_} = 1;
 }
}
elsif ($exclude_file) {
 $exclude_ids{$exclude_file} = 1;
}
print "Will exclude " . scalar( keys %exclude_ids ) . " queries"
  if %exclude_ids;

my $exonerate = `which exonerate`;
chomp($exonerate);
die "Exonerate not found in path\n" unless $exonerate && -e $exonerate;

# make indexes
my ( $inx, $ref_inx );
unless ( -f "$fasta_file.index" ) {
 print "Indexing $fasta_file...\n";
 my $i =
   Bio::Index::Fasta->new( -filename => "$fasta_file.index", -write_flag => 1 );
 $i->make_index($fasta_file);
}
$inx = Bio::Index::Fasta->new( -filename => "$fasta_file.index" )
  || die("Could not get index for $fasta_file\n");

unless ( -f "$reference_file.index" ) {
 print "Indexing $reference_file...\n";
 my $i = Bio::Index::Fasta->new( -filename   => "$reference_file.index",
                                 -write_flag => 1 );
 $i->make_index($reference_file);
}
$ref_inx = Bio::Index::Fasta->new( -filename => "$reference_file.index" )
  || die("Could not get index for $reference_file\n");

print "Searching $fasta_file vs $reference_file\n";
if (@aat_filter_files) {
 print "Temporary output from dps. Using it...\n";
 foreach my $file (@aat_filter_files) {
  &parse_aat_filter($file);
 }
}
elsif ( !$blast_file ) {
 my $reference_base = basename($reference_file);
 my $base           = basename($fasta_file);
 my $blast_file     = $base . '_vs_' . $reference_base . '.blast';
 if ( !-s $blast_file ) {
  print "Blasting input FASTA vs reference\n";
  system(
"makeblastdb -in $reference_file -dbtype nucl -title $reference_file -parse_seqids -out $reference_file"
  ) unless -s $reference_file . '.nin';
  system(
"blastn -outfmt 6 -max_target_seqs 1  -evalue $eval_cut -task megablast -num_threads $threads -dust no -query $fasta_file -db $reference_file -out $blast_file"
  ) if !$is_protein;
  system(
"tblastn -outfmt 6 -max_target_seqs 1 -evalue $eval_cut -num_threads $threads -seg no -query $fasta_file -db $reference_file -out $blast_file"
  ) if $is_protein;
 }
 &parse_analyze_blast($blast_file);
}
elsif ( !-s $blast_file ) {
 die "Could not find blast file $blast_file";
}
else {
 &parse_analyze_blast($blast_file);
}

############################
if ( !$is_protein ) {
 unless ( -s "$fasta_file.orf" ) {
  print "Preparing cDNA file ORFs\n";
  system("getorf $fasta_file $fasta_file.orf.u  -minsize $minorf -find 3");
  system("fasta_formatter -i $fasta_file.orf.u  -o $fasta_file.orf");
  unlink("$fasta_file.orf.u");
  die "Cannot translate $fasta_file\n" unless ( -s "$fasta_file.orf" );
 }
 print "Parsing ORFs\n";
 open( ORF, "$fasta_file.orf" );
 while ( my $ln = <ORF> ) {
  if ( $ln =~ /^>(\S+)_(\d+)\s\[(\d+)\s\-\s(\d+)\]/ ) {
   my $id     = $1;
   my $length = $queries{$id}{'length'};

   # some sequences may not have a blast hit:
   next if !$length;
   my $orf    = $2;
   my $start  = $3;
   my $end    = $4;
   my $strand = '+';
   my $size   = abs( $end - $start ) + 1;
   if ( $end < $start ) {
    $strand = '-';
    $start  = $length - $start + 1;
    $end    = $start + $size - 1;
   }

   # does not add stop codon
   $end  += 3;
   $size += 3;
   my $seq = <ORF>;
   if ( exists $queries{$id} && exists $queries{$id}{'orf'}{'size'} ) {

    # longest only
    if ( $queries{$id}{'orf'}{'size'} < $size ) {
     delete $queries{$id}{'orf'};
    }
    else {
     next;
    }
   }
   chomp($seq);
   $queries{$id}{'orf'}{'size'}   = $size;
   $queries{$id}{'orf'}{'start'}  = $start;
   $queries{$id}{'orf'}{'end'}    = $end;
   $queries{$id}{'orf'}{'seq'}    = $seq;
   $queries{$id}{'orf'}{'strand'} = $strand;
  }
 }
 close ORF;
 open( ANNOT, ">$fasta_file.annotations.txt" );
 foreach my $id ( keys %queries ) {
  print ANNOT "$id "
    . $queries{$id}{'orf'}{'strand'} . ' '
    . $queries{$id}{'orf'}{'start'} . ' '
    . $queries{$id}{'orf'}{'size'} . "\n"
    if $queries{$id}{'orf'}{'size'};
 }
 close ANNOT;
}
else {

 #is protein, add stop codon at the end
 #min aa = $minaa
 foreach my $id ( keys %queries ) {
  delete $queries{$id} if $queries{$id}{'length'} < $minaa;
  $queries{$id}{'seq'} .= '*' unless $queries{$id}{'seq'} =~ /\*$/;

 }
}

#exonerate efficiency
my ( %reference_hash, $no_orf );
foreach my $id ( keys %queries ) {
 unless ( $queries{$id}{'seq'}
          && ( $is_protein || exists $queries{$id}{'orf'} ) )
 {

  #    warn "Query $id: no ORF found\n";
  $no_orf++;
  delete $queries{$id};
  next;
 }
 next unless $queries{$id}{'hit'};
 push( @{ $reference_hash{ $queries{$id}{'hit'} } }, $id );
}
die "No references have been matched.\n"
  unless scalar( keys %reference_hash ) > 0;
warn "$no_orf queries had no ORF\n" if $no_orf;
my $reference_dir = $reference_file . "_dir1/";
unless ( -d $reference_dir ) {
 print "Splitting reference file\n";
 my $files_ref = &splitfasta($reference_file,$reference_dir,1);
 die "Cannot split $reference_file\n" unless $files_ref && scalar(@$files_ref>0);
}

my $query_dir = basename($fasta_file) . '_queries';

if ( -d $query_dir ) {
 print "Cleaning uncompleted output files from previous runs\n";
 my @todelete1 = glob( $fasta_file . "_queries/*results" );
 my @todelete2 = glob( $fasta_file . "_queries/*queries" );
 foreach my $fasta_file (@todelete1) {
  my $check = `tail -n 1 $fasta_file`;
  unlink($fasta_file) unless $check =~ /completed exonerate analysis\s*$/;
 }
 foreach my $fasta_file (@todelete2) {
  unlink($fasta_file) unless -s $fasta_file . '.exonerate_results';
 }
}
else {
 mkdir($query_dir);
}
unlink("exonerate_results.txt");

&separate_exonerate() if $separate;
&scaffold_exonerate() if !$separate;
unlink("multithread.cmd")  unless -s "multithread.cmd";
die "No commands to run\n" unless -s "multithread.cmd";
if ($shuf_exec) {
 chomp($shuf_exec);
 system("$shuf_exec multithread.cmd > multithread.cmd.");
 rename( "multithread.cmd.", "multithread.cmd" );
}

if ($prep_only) {
 if ( $prep_only > 1 ) {
  my $lns = `wc -l < multithread.cmd`;
  chomp($lns);
  my $file_lns = int( ( $lns + 10 ) / $prep_only );
  system("rm -f multithread.cmd.? multithread.cmd.??");
  system( "split -d -l $file_lns -a 3 multithread.cmd multithread.cmd." );
  warn "Data prepared, see multithread.cmd.*\n";
 }
 else {
  warn "Data prepared, see multithread.cmd\n";
 }
 exit();
}
print "Starting exonerate with $threads runs versus "
  . scalar( keys %reference_hash )
  . " reference targets in $fasta_file"
  . "_queries/\n";
system(
"$cmd_processor -shuffle -CPU $threads -c multithread.cmd -failed_cmds multithread.cmd.failed -v"
) if -s "multithread.cmd";
print "Finished!\n" if -s "multithread.cmd";
print "Something went wrong or no exonerates jobs to run run...\n"
  if !-s "multithread.cmd";

######################################3
sub parse_analyze_blast() {
 my $blast_file = shift;
 print "Parsing BLAST tabular file\n";
 open( BLAST, $blast_file );
 while ( my @data = split( "\t", <BLAST> ) ) {
  chomp( $data[-1] );
  my $query_id = $data[0];
  next if exists $exclude_ids{$query_id};
  next if $data[-1] < $dps_min_score;
  next
    if $queries{$query_id}{'score'} && $data[-1] < $queries{$query_id}{'score'};

  my $gene_obj = $inx->fetch($query_id);
  if ( !$gene_obj ) {
   warn "Sequence not found for $query_id. Skipping...\n";
   next;
  }

  $queries{$query_id}{'score'}  = $data[-1];
  $queries{$query_id}{'seq'}    = $gene_obj->seq();
  $queries{$query_id}{'length'} = $gene_obj->length();
  $queries{$query_id}{'hit'}    = $data[1];
 }

 close BLAST;
}

#################3
sub parse_aat_filter() {
 my $file = shift || die;
 open( IN, $file ) || die "Cannot open $file";
 $file =~ /^([^\.]+)\./;    #TODO GROSS ASSUMPTION no other way!!!
 my $ref = $1 || die;
 print "Assuming file's reference is $ref\n";

 #   $dstart, $dend, $score, $astart, $aend, $orient, $zero1, $zero2, $acc
 # want $acc, and $dstart, $dend
 my $header = <IN>;
 while ( my $ln = <IN> ) {
  $ln =~
    /^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)/;
  my ( $dstart, $dend, $score, $astart, $aend, $orient, $zero1, $zero2, $query )
    = ( $1, $2, $3, $4, $5, $6, $7, $8, $9 );
  die unless $query;

  #print "Processing $query\n";
  next if $score < $dps_min_score;
  next if $queries{$query}{'score'} && $score < $queries{$query}{'score'};
  my $gene_obj = $inx->fetch($query);
  if ( !$gene_obj ) {
   warn "Sequence not found for $query. Skipping...\n";
   $exclude_ids{$query} = 1;
   next;
  }

  $queries{$query}{'seq'}         = $gene_obj->seq();
  $queries{$query}{'length'}      = $gene_obj->length();
  $queries{$query}{'hit'}         = $ref;
  $queries{$query}{'score'}       = $score;
  $queries{$query}{'lower_bound'} = $dstart;
  $queries{$query}{'upper_bound'} = $dend;

 }
 close IN;
}

sub parse_aat_core() {

=cut
filter is this 
                       327013     3162
    2194   301443    103     121   333 1     0     0 gi|157103285|ref|XP_001647909.1|
    3862   273047    106     146   293 1     0     0 gi|157103227|ref|XP_001647880.1|
    7015   189003    122     151   251 1     0     0 gi|157103225|ref|XP_001647879.1|
    7945   299406    116     102   263 0     0     0 gi|157103213|ref|XP_001647874.1|
   11230   232423    118    1327  1749 1     0     0 gi|157103161|ref|XP_001647848.1|
   14497   324075    111     181   355 1     0     0 gi|157103219|ref|XP_001647877.1|
   16461   251710    103     155   332 0     0     0 gi|157103295|ref|XP_001647914.1|
   20448   282513    118      68   974 1     0     0 gi|157103179|ref|XP_001647857.1|
   26641   292605    134     188   348 0     0     0 gi|157103229|ref|XP_001647881.1|
   29869   283760    103     373   510 0     0     0 gi|157103446|ref|XP_001647986.1|
   88778   280235    167    1087  1616 1     0     0 gi|157103363|ref|XP_001647947.1|
   90202   260775    118      56   180 0     0     0 gi|157103307|ref|XP_001647920.1|
  114549   274786    105     283   341 1     0     0 gi|157103261|ref|XP_001647897.1|
  140459   325736    109      29   144 1     0     0 gi|157103283|ref|XP_001647908.1|
  189220   301479    107     330   413 1     0     0 gi|157103479|ref|XP_001647999.1|

but lets use the DPS because we could cocncat the files (filter does not report scaffold name
Query sequence (top one on the alignment): >scaffold_266
Chain   40658   221133    105    8      68   235 gi|157103179|ref|XP_001647857.1|
need: 
  $queries{$query_id}{'seq'}= sequence
  $queries{$query_id}{'length'}= gene length
  $queries{$query_id}{'hit'}= sacffold

=cut

 # query is gene and ref is scaffold (but dps reports opposite)
 my $file = shift || die;
 open( IN, $file ) || die "Cannot open $file";
 my $ref;
 while ( my $ln = <IN> ) {
  next unless $ln =~ /^Chain|^Query sequence \(top one/;
  chomp($ln);
  if ( $ln =~ /^Query sequence/ ) {
   $ln =~ />(\S+)$/;
   $ref = $1 || die;
  }
  else {
   $ln =~ /(\S+)$/;
   my $query = $1 || die;
   next if exists $exclude_ids{$query};
   $ln =~ /^Chain\s+(\d+)\s+(\d+)\s+(\d+)/;
   my $s     = $1;
   my $e     = $2;
   my $score = $3;
   next if $score < $dps_min_score;
   next if $queries{$query}{'score'} && $score < $queries{$query}{'score'};
   my $strand = 1;

   if ( $s > $e ) {
    $strand = -1;
    my $t = $e;
    $e = $s;
    $s = $t;
   }

   my $gene_obj = $inx->fetch($query);
   if ( !$gene_obj ) {
    warn "Sequence not found for $query. Skipping...\n";
    $exclude_ids{$query} = 1;
    next;
   }
   $queries{$query}{'seq'}         = $gene_obj->seq();
   $queries{$query}{'length'}      = $gene_obj->length();
   $queries{$query}{'hit'}         = $ref;
   $queries{$query}{'score'}       = $score;
   $queries{$query}{'lower_bound'} = $s;
   $queries{$query}{'upper_bound'} = $e;
  }
 }
 close IN;
}

sub scaffold_exonerate() {
 open( CMD, ">multithread.cmd" );

 #only process those with ORF
 foreach my $target ( sort keys %reference_hash ) {
  my @queries    = @{ $reference_hash{$target} };
  my $query_file = $fasta_file . '_queries/' . $target . '.queries';
  open( QUERY, ">$query_file" );
  foreach my $id (@queries) {
   if ($coding_only) {
    next
      if !$queries{$id}{'orf'}{'seq'} || $queries{$id}{'orf'}{'seq'} =~ /^\s*$/;
    print QUERY ">$id\n" . $queries{$id}{'orf'}{'seq'} . "\n" if $coding_only;
   }
   else {
    next if !$queries{$id}{'seq'} || $queries{$id}{'seq'} =~ /^\s*$/;
    print QUERY ">$id\n" . $queries{$id}{'seq'} . "\n";
   }
  }
  close QUERY;
  my $target_file = $reference_dir . "/$target";
### -D -M control memory usage!
  #"echo '$header' > $query_file.exonerate_results && $exonerate"
  my $common =
      "$exonerate "
    . ' --ryo "RYOAP_START\nRYOAP_STATS_D\tAlignment length\tidentical\tsimilar\tmismatch\tidentical%%\tsimilar%%\nRYOAP_STATS\t%s\t%et\t%ei\t%es\t%em\t%pi\t%ps\nRYOAP_CODING_QUERY\t%qi\t%qab\t%qae\n>%qi\n%qas\nRYOAP_CODING_GENOME\t%ti\t%tcb\t%tce\n>%ti\n%tcs\nRYOAP_END\n" '
    . " --targettype dna -D 8192 -M 8192 --showquerygff y --showtargetgff y -n 1 -S n --maxintron $max_intron  ";
  my $cmd;
  if ($same_species) {
   $cmd =
     $coding_only
     ? " -s 2000 -Q dna -m coding2genome --percent 90"
     : " -s 2000 -Q dna -m cdna2genome --annotation $fasta_file.annotations.txt";
   $cmd =
     " -s 1000 -Q protein --exhaustive 1 -m protein2genome:bestfit --percent 70"
     if $is_protein;
   $cmd .= " --hspfilter 100 --geneseed 250 ";
  }
  else {

   # I DON'T KNOW WHAT GOOD SCORE/PERCENT VALUES ARE HERE
   $cmd =
     $coding_only
     ? " -s 200 -Q dna -m coding2genome --percent 30"
     : " -s 200 -Q dna -m cdna2genome --annotation $fasta_file.annotations.txt";
   $cmd =
     " -s 100 -Q protein --exhaustive 1 -m protein2genome:bestfit --percent 20 "
     if $is_protein;
  }
  $cmd .= ' --refine region ' if !$no_refine && !$is_protein;
  $cmd .= " --softmasktarget y " if $softmasked;
  $cmd .= " $query_file $target_file ";
  print CMD (

   #"test ! -e $query_file.exonerate_results  && " .
   $common . $cmd . " > $query_file.exonerate_results\n"
  );
 }
 close(CMD);

}

sub separate_exonerate() {
 open( CMD, ">multithread.cmd" );

 #only process those with ORF
 foreach my $target ( sort keys %reference_hash ) {
  my @queries         = @{ $reference_hash{$target} };
  my $base_query_file = $fasta_file . '_queries/' . $target . '.queries';
  my $scaffold_obj    = $ref_inx->fetch($target);
  if ( !$scaffold_obj ) {
   warn "Cannot find reference sequence file for $target. Skipping...\n";
   next;
  }

  foreach my $id (@queries) {
   my $start                = $queries{$id}{'lower_bound'};
   my $end                  = $queries{$id}{'upper_bound'};
   my $genome_region_length = $scaffold_obj->length();
   if ( $genome_region_length > ( $end + $padding ) ) {
    $end += $padding;
   }
   else {
    $end = $genome_region_length;
   }
   if ( $start <= $padding ) {
    $start = 1;
   }
   else {
    $start -= $padding;
   }
   my $genome_region = $scaffold_obj->subseq( $start, $end );
   my $comp_q = $id;
   $comp_q =~ s/\W+/_/g;
   $comp_q =~ s/__+/_/g;
   $comp_q =~ s/_$//g;
   $comp_q =~ s/^_//g;

#        my $header= "Processing reference target $target with query $comp_q\n";
   my $query_file = $base_query_file . "_$comp_q";
   next if -s $query_file;
   open( QUERY, ">$query_file" );
   print QUERY ">$id\n" . $queries{$id}{'seq'} . "\n" if !$coding_only;
   print QUERY ">$id\n" . $queries{$id}{'orf'}{'seq'} . "\n" if $coding_only;
   close QUERY;

   my $target_file = "$query_file.scaffold";
   open( SCAFFOLD, ">$target_file" );
   print SCAFFOLD ">$target" . ":$start-$end\n$genome_region\n";
   close SCAFFOLD;

### -D -M control memory usage!
   #"echo '$header' > $query_file.exonerate_results && $exonerate"
   my $common =
       "$exonerate "
     . ' --ryo "RYOAP_START\nRYOAP_STATS_D\tAlignment length\tidentical\tsimilar\tmismatch\tidentical%%\tsimilar%%\nRYOAP_STATS\t%s\t%et\t%ei\t%es\t%em\t%pi\t%ps\nRYOAP_CODING_QUERY\t%qi\t%qab\t%qae\n>%qi\n%qas\nRYOAP_CODING_GENOME\t%ti\t%tcb\t%tce\n>%ti\n%tcs\nRYOAP_END\n" '
     . " --targettype dna -D 8192 -M 8192 --showtargetgff y -n 1 -S n --maxintron $max_intron  ";
   my $cmd;
   if ($same_species) {
    $cmd =
      $coding_only
      ? " -s 2000 -Q dna -m coding2genome --percent 90"
      : " -s 2000 -Q dna -m cdna2genome --annotation $fasta_file.annotations.txt";
    $cmd =
" -s 1000 -Q protein --exhaustive 1 -m protein2genome:bestfit --percent 70"
      if $is_protein;
    $cmd .= " --hspfilter 100 --geneseed 250 ";
   }
   else {

    # I DON'T KNOW WHAT GOOD SCORE/PERCENT VALUES ARE HERE
    $cmd =
      $coding_only
      ? " -s 200 -Q dna -m coding2genome --percent 30"
      : " -s 200 -Q dna -m cdna2genome --annotation $fasta_file.annotations.txt";
    $cmd =
" -s 100 -Q protein --exhaustive 1 -m protein2genome:bestfit --percent 20 "
      if $is_protein;
   }
   $cmd .= ' --refine region ' if !$no_refine && !$is_protein;
   $cmd .= " --softmasktarget y " if $softmasked;
   $cmd .= " $query_file $target_file ";
   print CMD (
    #"test ! -e $query_file.exonerate_results  && " .

                     $common
      . $cmd . " > $query_file.exonerate_results\n"
   );
  }
 }
 close(CMD);
}

sub splitfasta(){
my @files;
my $file2split=shift;
my $outdir = shift;
my $how_many_in_a_file = shift; 
return unless $file2split && -s $file2split && $outdir && $how_many_in_a_file;
mkdir($outdir) unless -d $outdir;   
my ($flag);
my $filecount=0;
my $seqcount=0;

open (FILE,$file2split);
while (my $line=<FILE>){
        if ($line=~/^\s*$/){next;}      #empty line
        elsif ($line=~/^>(\S+)/){
                $seqcount++;
                if (!$flag){
                        $filecount++;
                        my $outfile=$file2split;
                        $outfile.="_".$filecount;
                        open (OUT, ">$outdir/$outfile")||die("Cannot open $outdir/$outfile");
                        push(@files,"$outfile");
                        $flag=1;
                }
                elsif ($seqcount>=$how_many_in_a_file){
                        $seqcount=0;
                        $filecount++;
                        my $outfile=$file2split;
                        $outfile.="_".$filecount;
                        close (OUT);
                        open (OUT, ">$outdir/$outfile");
                        push(@files,"$outfile");
                 }
                 print OUT $line;
         }
        else { print OUT $line; }
}
close (FILE);
close (OUT);
return \@files;
}
                        
