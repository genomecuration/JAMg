#!/usr/bin/env perl

=pod

=head1 NAME

run_exonerate

=head1 USAGE

Run exonerate
 

 @ Stage 1:
 *  -in|fasta:s  => Fasta file (cDNA or protein)
 *  -reference:s => Fasta with reference sequences (genome)
  -protein      => File is a protein (otherwise assumed to be cDNA)
  -minorf:i     => minimum ORF size for cDNA
  -annotation   => Define ORF in cDNA input (rather that find longest)
  -blast_hash:s => Blast HASH file from analyze_blast
  -only_coding  => For cDNA: Use only the coding sequence (coding2genome model). Otherwise (default) use cdna2genome (includes UTR)
  -exclude:s    => Exclude file (one ID per line) or comma separated string of IDs 
  -threads      => Threads/CPUs to use
  -eval_cutoff:s=> Evalue cutoff for BLAST (def 1e-80)
  -bit_cutoff:i => Bit score cutoff for BLAST (def 100 for protein tblastn, 500 for nucleotide megablast)
  -prep_only:i  => Set to 1 or more. Do BLAST and prepare commands only. If number >1 is given, then split cmds to that many subfiles (e.g. to run with Parafly)
  -norefine     => Don't refine region in exonerate. Prevents segmentation faults in some searches
  -aat:s        => AAT FILTER output files (one or more). Do not use blast hash
  -score_dps:i  => Minimum DDS/DPS score (def 100)
  -intron_max   => Def 70000
  -same_species => FASTA file is from same species as referene genome
  -softmask     => genome has been softmasked with repeatmasker
  -local_protein => Don't run an exhaustive global alignment for the protein (use it for those that crash)
  -separate      => Only run exonerate on part of the assembly, padded with -padding. Very fast!
  -padding      => genome padding for exonerate alignment (2000)

 @ Stage 2:
  -postprocess:s{,}' => The result files of exonerate if run with -separate in orted to reconstruct co-ordinates of the GFF (and only the GFF)
  -postscore:i  => require at least this score before reporting in post-processing (def. 500) 
 
 Example cmd with existing AAT filter output:
  $ run_exonerate.pl -in arab.pep -ref arab.genomicSeq -thread 4 -same -protein -aat *filter -separ
  # this will run it locally with 4 threads. Use -prep_only to just create the command file but not run it.
 NB: Protein exhaustive alignments sometimes fail. So re-run with local alignment (good output is not overwritten).
  $ run_exonerate.pl -in arab.pep -ref arab.genomicSeq -thread 4 -same -protein -aat *filter -separ -local_protein
  # post-process results to correct co-ordinates.
  $ run_exonerate.pl -postprocess arab.pep_queries/*results > exonerate.result

 NB: Protein option works well. Exonerate often has issues with cDNA.

=head1 DEPENDENCIES

 * exonerate (obviously)
 * emboss
 * NCBI's BLAST C++ 
 
=head1 SUGGESTIONS
  
If you use RepeatMasker then:
$ RepeatMasker -e abblast -pa 35 -s -species insecta -no_is -gccalc -gff  -xsmall 
 
=head1 AUTHORS

 Alexie Papanicolaou

        CSIRO Ecosystem Sciences
        alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

Copyright 2012-2014 the Commonwealth Scientific and Industrial Research Organization. 
See LICENSE file for license info
It is provided "as is" without warranty of any kind.

=head1 BUGS & LIMITATIONS

None known, please report.

=cut

use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use File::Basename;
use Getopt::Long;
use List::Util 'shuffle';

use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";
use CdbTools;

my (
     $fasta_file,                    $reference_file,
     $blast_file,                    $coding_only,
     $help,                          $exclude_file,
     $is_protein,                    $prep_only,
     $no_refine,                     @aat_filter_files,
     $same_species,                  $softmasked,
     @post_process_similarity_files, $separate,
     @post_process_result_files,     $not_exhaustive,
     $annotation_file,               $paths_need_to_be_made
);

#variables
my $minorf              = 290;       #minimum orf size in bp
my $threads             = 4;
my $exonerate_postscore = 500;
my $eval_cut            = '1e-80';
my $bit_cut             = 500;
my $dps_min_score       = 100;
my $max_intron          = 70000;
my $padding             = 2000;

#globals
my (%get_id_seq_from_fasta_hash);
my ( $cdbfasta_exec, $cdbyank_exec, $parafly_exec, $exonerate ) =
  &check_program( 'cdbfasta', 'cdbyank', 'ParaFly', 'exonerate' );
my ( $makeblastdb_exec, $blastn_exec, $tblastn_exec ) =
  &check_program( 'makeblastdb', 'blastn', 'tblastn' );
my ($getorf_exec) = &check_program('getorf');

pod2usage $! unless &GetOptions(
         'help'             => \$help,
         'annotations:s'    => \$annotation_file,
         'minorf:i'         => \$minorf,
         'in|fasta:s'       => \$fasta_file,
         'reference:s'      => \$reference_file,
         'blast_hash:s'     => \$blast_file,
         'only_coding'      => \$coding_only,
         'exclude:s'        => \$exclude_file,
         'threads:i'        => \$threads,
         'protein'          => \$is_protein,
         'eval_cutoff:s'    => \$eval_cut,
         'bit_cutoff:i'     => \$bit_cut,
         'prep_only:i'      => \$prep_only,
         'norefine'         => \$no_refine,
         'aat|filter:s{,}'  => \@aat_filter_files,
         'score_dps:i'      => \$dps_min_score,
         'intron_max:i'     => \$max_intron,
         'same_species'     => \$same_species,
         'softmask'         => \$softmasked,
         'padding:i'        => \$padding,
         'separate'         => \$separate,
         'postprocess:s{,}' => \@post_process_result_files,
         'similarity:s{,}'  => \@post_process_similarity_files,    #from bioperl
         'local_protein'    => \$not_exhaustive,
         'postscore:i'      => \$exonerate_postscore,
         'from_blast'       => \$paths_need_to_be_made,
);
pod2usage if $help;

if (@post_process_result_files) {
 my $file_number = scalar(@post_process_result_files);
 warn "Post-processing $file_number exonerate result files...\n";
 &process_exonerate( \@post_process_result_files );
 exit(0);
}
elsif (@post_process_similarity_files) {
 my $file_number = scalar(@post_process_similarity_files);
 warn "Post-processing $file_number exonerate result files...\n";
 &process_exonerate_similarity( \@post_process_similarity_files );
 exit(0);
}
pod2usage
  unless $fasta_file && $reference_file && -s $fasta_file && -s $reference_file;
$bit_cut = 100 if $is_protein && $bit_cut == 500;
my $minaa = int( $minorf / 3 );

# make indexes
print "Indexing...\n";
&process_cmd("$cdbfasta_exec $fasta_file");
&process_cmd("$cdbfasta_exec $reference_file");
my $contig_seq_hash  = &read_fasta($fasta_file);


my ( %queries, %exclude_ids );

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
 my $blast_file     = $base . '_vs_' . $reference_base;
 $blast_file .= '.tblastn' if $is_protein;
 $blast_file .= '.blastn'  if !$is_protein;
 if ( !-s $blast_file ) {
  print "Blasting input FASTA vs reference\n";
  &process_cmd(
"$makeblastdb_exec -in $reference_file -dbtype nucl -title $reference_file -parse_seqids -out $reference_file"
  ) unless -s $reference_file . '.nin';
  &process_cmd(
"$blastn_exec -outfmt 6 -max_target_seqs 1  -evalue $eval_cut -task megablast -num_threads $threads -dust no -query $fasta_file -db $reference_file -out $blast_file"
  ) if !$is_protein;
  &process_cmd(
"$tblastn_exec -outfmt 6 -max_target_seqs 1 -evalue $eval_cut -num_threads $threads -seg no -query $fasta_file -db $reference_file -out $blast_file"
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

print "Preprocessing for exonerate run...\n";

############################
my $not_used;

#clear out
my $query_dir = basename($fasta_file) . '_queries';
if ( -d $query_dir ) {
 print "Cleaning uncompleted output files from previous runs\n";
 my @todelete1 = glob( $fasta_file . "_queries/*results" );
 if (@todelete1) {
  foreach my $fasta_file (@todelete1) {
   my $check = `tail -n 1 $fasta_file`;
   unlink($fasta_file)
     unless $check =~ /completed exonerate analysis\s*$/;
  }
  my @todelete2 = glob( $fasta_file . "_queries/*queries" );
  foreach my $fasta_file (@todelete2) {
   unlink($fasta_file) unless -s $fasta_file . '.exonerate_results';
  }
 }
}
else {
 mkdir($query_dir);
}

if ( !$is_protein ) {
 if ($annotation_file) {
  open( IN, $annotation_file );
  while ( my $ln = <IN> ) {
   my @data = split( /\s+/, $ln );
   if ( $data[3] < $minorf ) {
    warn $data[0] . " is smaller than $minorf bp. Skipping\n";
    next;
   }
   $queries{ $data[0] }{'orf'}{'strand'} = $data[1];
   $queries{ $data[0] }{'orf'}{'start'}  = $data[2];
   $queries{ $data[0] }{'orf'}{'size'}   = $data[3];

  }
  close IN;
 }
 else {
  unless ( -s "$fasta_file.orf" ) {
   print "Preparing cDNA file ORFs\n";
   &process_cmd(
       "$getorf_exec $fasta_file $fasta_file.orf.u  -minsize $minorf -find 3" );
   &convert_fasta_to_single_line("$fasta_file.orf.u");
   rename( "$fasta_file.orf.u", "$fasta_file.orf.orf" );
   die "Cannot translate $fasta_file\n"
     unless ( -s "$fasta_file.orf" );
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
    if (    exists $queries{$id}
         && exists $queries{$id}{'orf'}{'size'} )
    {

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
  $annotation_file = "$fasta_file.annotations";
  open( ANNOT, ">$annotation_file" );

  foreach my $id ( keys %queries ) {
   print ANNOT "$id "
     . $queries{$id}{'orf'}{'strand'} . ' '
     . $queries{$id}{'orf'}{'start'} . ' '
     . $queries{$id}{'orf'}{'size'} . "\n"
     if $queries{$id}{'orf'}{'size'};
  }
  close ANNOT;
 }
}
else {

 #is protein, add stop codon at the end
 #min aa = $minaa
 foreach my $id ( keys %queries ) {
  if ( $queries{$id}{'length'} < $minaa ) {
   warn "$id is smaller than $minaa amino acids. Skipping\n";
   delete $queries{$id};
   $not_used++;
   next;
  }
  $queries{$id}{'seq'} .= '*'
    unless !$queries{$id}{'seq'} || $queries{$id}{'seq'} =~ /\*$/;
  # convert U to X
  $queries{$id}{'seq'} =~ tr/U/X/;
 }
}

#exonerate efficiency
my (%reference_hash);
open( NOTUSED, ">$fasta_file.annotations.notused" );
foreach my $id ( keys %queries ) {
 unless (    $queries{$id}{'hit'}
          && $queries{$id}{'seq'}
          && ( $is_protein || exists $queries{$id}{'orf'} ) )
 {

  #warn "Query $id: no hit, or ORF found\n";
  print NOTUSED "Query $id: no hit, or ORF found\n";
  my $err = Dumper $queries{$id};
  print NOTUSED $err . "\n";
  $not_used++;
  delete $queries{$id};
  next;
 }
 push( @{ $reference_hash{ $queries{$id}{'hit'} } }, $id );
}
close NOTUSED;
die "No references have been matched.\n"
  unless scalar( keys %reference_hash ) > 0;
warn "$not_used queries had no hit or ORF and will not be searched\n"
  if $not_used;
my $reference_dir = basename($reference_file) . "_dir/";
unless ( -d $reference_dir ) {
 print "Splitting reference file\n";
 my $files_ref = &splitfasta( $reference_file, $reference_dir, 1 );
 die "Cannot split $reference_file\n"
   unless $files_ref && scalar( @$files_ref > 0 );
}
print "Preparing/processing files\n";
if ($separate) {
 &separate_exonerate();
}
else {
 &scaffold_exonerate();
}
unlink("run_exonerate_commands.cmd") unless -s "run_exonerate_commands.cmd";
warn "No commands to run\n"          unless -s "run_exonerate_commands.cmd";
exit(0)                              unless -s "run_exonerate_commands.cmd";

open( SHUF,    "run_exonerate_commands.cmd" );
open( SHUFOUT, ">run_exonerate_commands.cmd." );
my @shuf = <SHUF>;
close SHUF;
print SHUFOUT shuffle(@shuf);
close SHUFOUT;
rename( "run_exonerate_commands.cmd.", "run_exonerate_commands.cmd" );

if ($prep_only) {
 if ( $prep_only > 1 ) {
  my $lns = `wc -l < run_exonerate_commands.cmd`;
  chomp($lns);
  my $file_lns = int( ( $lns + 10 ) / $prep_only );
  system("rm -f run_exonerate_commands.cmd.? run_exonerate_commands.cmd.??");
  system(
"split -d -l $file_lns -a 3 run_exonerate_commands.cmd run_exonerate_commands.cmd."
  );
  warn "Data prepared, see run_exonerate_commands.cmd.*\n";
 }
 else {
  warn "Data prepared, see run_exonerate_commands.cmd\n";
 }
 exit();
}
print "Starting exonerate with $threads runs versus "
  . scalar( keys %reference_hash )
  . " reference targets in $fasta_file"
  . "_queries/\n";
&process_cmd(
"$parafly_exec -shuffle -CPU $threads -c run_exonerate_commands.cmd -failed_cmds run_exonerate_commands.cmd.failed -v"
) if -s "run_exonerate_commands.cmd";
print "Finished!\n" if -s "run_exonerate_commands.cmd";
print "Something went wrong or no exonerates jobs to run run...\n"
  if !-s "run_exonerate_commands.cmd";
######################################3
sub parse_analyze_blast() {
 
 # not to be used
 my $blast_file = shift;
 print "Parsing BLAST tabular file\n";
 open( BLAST, $blast_file );
 while ( my @data = split( "\t", <BLAST> ) ) {
  chomp( $data[-1] );
  my $query_id = $data[0];
  next if exists $exclude_ids{$query_id};
  next if $data[-1] < $dps_min_score;
  next
    if $queries{$query_id}{'score'}
     && $data[-1] < $queries{$query_id}{'score'};
  my $seq = &get_id_seq_from_fasta( $query_id, $fasta_file );
  if ( !$seq ) {
   #warn "Sequence not found for $query_id. Skipping...\n";
   next;
  }
  $queries{$query_id}{'score'}  = $data[-1];
  $queries{$query_id}{'seq'}    = $seq;
  $queries{$query_id}{'length'} = length($seq);
  $queries{$query_id}{'hit'}    = $data[1];
  $queries{$query_id}{'lower_bound'} =
    $data[8] < $data[9] ? $data[8] : $data[9];
  $queries{$query_id}{'upper_bound'} =
    $data[8] < $data[8] ? $data[9] : $data[8];

 }
 close BLAST;
}
#################3
sub parse_aat_filter() {
 my $file = shift || die;

 my $basename = basename($file);
 open( IN, $file ) || die "Cannot open $file";
 $basename =~ /^(.+)\.\w+\.filter/;    #TODO GROSS ASSUMPTION!!
 my $ref = $1 || die;

 my ( %hash, %top_score );

 #  warn "Assuming file's reference is $ref\n";

 #   $dstart, $dend, $score, $astart, $aend, $orient, $zero1, $zero2, $acc
 # want $acc, and $dstart, $dend
 my $header = <IN>;
 while ( my $ln = <IN> ) {
  $ln =~
    /^\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)/;
  my ( $dstart, $dend, $score, $astart, $aend, $orient, $zero1, $zero2, $query )
    = ( $1, $2, $3, $4, $5, $6, $7, $8, $9 );
  die $ln unless $query;
  next if exists $exclude_ids{$query};
  next if $score < $dps_min_score;

  if ( $dend < $dstart ) {
   my $t = $dstart;
   $dstart = $dend;
   $dend   = $t;
  }
  if ($paths_need_to_be_made) {
   # this still doesn't make paths.
   $hash{$query}{$ref}{'score'}  = $score;
   $hash{$query}{$ref}{'dstart'} = $dstart
     if !$hash{$query}{$ref}{'dstart'}
      || $hash{$query}{$ref}{'dstart'} > $dstart;
   $hash{$query}{$ref}{'dend'} = $dend
     if !$hash{$query}{$ref}{'dend'} || $hash{$query}{$ref}{'dend'} < $dend;
   if ( !$top_score{$query} ) {
    $top_score{$query}{'ref'}   = $ref;
    $top_score{$query}{'score'} = $score;
   }
   else {
    if ( $top_score{$query}{'score'} < $score ) {
     $top_score{$query}{'ref'}   = $ref;
     $top_score{$query}{'score'} = $score;
    }
   }
  }
  else {
   next
     if (    $top_score{$query}{'ref'}
          && $score < $top_score{$query}{'score'} );
   $hash{$query}{$ref}{'score'}  = $score;
   $hash{$query}{$ref}{'dstart'} = $dstart;
   $hash{$query}{$ref}{'dend'}   = $dend;
   $top_score{$query}{'ref'}   = $ref;
   $top_score{$query}{'score'} = $score;
  }

 }
 close IN;
# print "Preparing hits from $file...\n";
 foreach my $query ( keys %hash ) {
  my $seq       = $contig_seq_hash->{$query};
  if ( !$seq ) {
   #warn "Sequence not found for $query. Skipping...\n";
   $exclude_ids{$query} = 1;
   next;
  }

  $queries{$query}{'seq'}         = $seq;
  $queries{$query}{'length'}      = length($seq);
  $queries{$query}{'hit'}         = $top_score{$query}{'ref'};
  $queries{$query}{'score'}       = $top_score{$query}{'score'};
  $queries{$query}{'lower_bound'} = $hash{$query}{$top_score{$query}{'ref'}}{'dstart'};
  $queries{$query}{'upper_bound'} = $hash{$query}{$top_score{$query}{'ref'}}{'dend'};
 }
}

sub parse_aat_core() {

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
   next
     if $queries{$query}{'score'}
      && $score < $queries{$query}{'score'};
   my $strand = 1;

   if ( $s > $e ) {
    $strand = -1;
    my $t = $e;
    $e = $s;
    $s = $t;
   }
   my $seq = &get_id_seq_from_fasta( $query, $fasta_file );
   if ( !$seq ) {
    #warn "Sequence not found for $query. Skipping...\n";
    $exclude_ids{$query} = 1;
    next;
   }
   $queries{$query}{'seq'}         = $seq;
   $queries{$query}{'length'}      = length($seq);
   $queries{$query}{'hit'}         = $ref;
   $queries{$query}{'score'}       = $score;
   $queries{$query}{'lower_bound'} = $s;
   $queries{$query}{'upper_bound'} = $e;
  }
 }
 close IN;
}

sub scaffold_exonerate() {
 open( CMD, ">run_exonerate_commands.cmd" );
 my $do_annotations = $annotation_file
   && -s $annotation_file ? "--annotation $annotation_file" : ' ';

 #only process those with ORF
 foreach my $scaffold ( shuffle keys %reference_hash ) {
  my @queries    = @{ $reference_hash{$scaffold} };
  my $query_file = $fasta_file . '_queries/' . $scaffold . '.queries';
  open( QUERY, ">$query_file" );
  foreach my $id (@queries) {
   if ($coding_only) {
    next
      if !$queries{$id}{'orf'}{'seq'}
       || $queries{$id}{'orf'}{'seq'} =~ /^\s*$/;
    print QUERY ">$id\n" . $queries{$id}{'orf'}{'seq'} . "\n"
      if $coding_only;
   }
   else {
    next
      if !$queries{$id}{'seq'} || $queries{$id}{'seq'} =~ /^\s*$/;
    print QUERY ">$id\n" . $queries{$id}{'seq'} . "\n";
   }
  }
  close QUERY;
  my $scaffold_file = $reference_dir . "/$scaffold";
### -D -M control memory usage!
  #"echo '$header' > $query_file.exonerate_results && $exonerate"
  my $common =
      "$exonerate "
    . ' --ryo "RYOAP_START\nRYOAP_STATS_D\tAlignment length\tidentical\tsimilar\tmismatch\tidentical%%\tsimilar%%\nRYOAP_STATS\t%s\t%et\t%ei\t%es\t%em\t%pi\t%ps\nRYOAP_CODING_QUERY\t%qi\t%qab\t%qae\n>%qi\n%qas\nRYOAP_CODING_GENOME\t%ti\t%tcb\t%tce\n>%ti\n%tcs\nRYOAP_END\n" '
    . " --targettype dna -D 8192 -M 8192 --showquerygff y --showtargetgff y -n 1 -S n --maxintron $max_intron  ";
  my $cmd;
  if ($same_species) {
   $cmd = $coding_only
     ? " -s 2000 -Q dna -m coding2genome  "    # removed --percent 90
     : " -s 2000 -Q dna -m cdna2genome  $do_annotations "
     ;                                         # removed --percent 70
   $cmd = " -s 1000 -Q protein --exhaustive 1 -m protein2genome:bestfit "
     if $is_protein && !$not_exhaustive;       # removed --percent 70
   $cmd = " -s 1000 -Q protein -m protein2genome "
     if $is_protein && $not_exhaustive;        # removed --percent 70
   $cmd .= " --hspfilter 100 --geneseed 250 ";
  }
  else {

   # I DON'T KNOW WHAT GOOD SCORE/PERCENT VALUES ARE HERE
   $cmd =
     $coding_only
     ? " -s 200 -Q dna -m coding2genome --percent 30"
     : " -s 200 -Q dna -m cdna2genome --percent 20 $do_annotations ";
   $cmd =
     " -s 100 -Q protein --exhaustive 1 -m protein2genome:bestfit --percent 30 "
     if $is_protein && !$not_exhaustive;
   $cmd = " -s 100 -Q protein -m protein2genome --percent 30 "
     if $is_protein && $not_exhaustive;
  }
  $cmd .= ' --refine region ' if !$no_refine && !$is_protein;
  $cmd .= " --softmasktarget y " if $softmasked;
  $cmd .= " $query_file $scaffold_file ";
  print CMD (

     $common 
     . $cmd
     . " > $query_file.exonerate_results\n"
  );
 }
 close(CMD);
}

sub separate_exonerate() {
 open( CMD, ">run_exonerate_commands.cmd" );
 my $do_annotations = $annotation_file
   && -s $annotation_file ? "--annotation $annotation_file" : ' ';
 my $reference_base = basename($reference_file);
 my $base           = basename($fasta_file);

 #only process those with ORF
 foreach my $scaffold ( sort keys %reference_hash ) {
  my @queries         = @{ $reference_hash{$scaffold} };
  my $base_query_file = $base . '_queries/' . $scaffold . '.queries';
  my $scaffold_seq    = &get_id_seq_from_fasta( $scaffold, $reference_file );
  if ( !$scaffold_seq ) {
   warn "Cannot find reference sequence file for $scaffold. Skipping...\n";
   next;
  }
  foreach my $id ( shuffle @queries ) {
   my $start = $queries{$id}{'lower_bound'};
   my $end   = $queries{$id}{'upper_bound'};
   if ($end < $start){
    my $t = $end;
    $end = $start;
    $start = $t;
   }
   if ( !$start ) {
    die "No lower bound for $id";

    #next;
   }
   if ( !$end ) {
    die "No upper bound for $id";

    #next;
   }
   my $scaffold_length = length($scaffold_seq);
   if ( $scaffold_length > ( $end + $padding ) ) {
    $end += $padding;
   }
   else {
    $end = $scaffold_length;
   }
   if ( $start <= $padding ) {
    $start = 1;
   }
   else {
    $start -= $padding;
   }
   my $genome_region_length = $end - $start + 1;     
   my $genome_region = substr( $scaffold_seq, $start - 1, $genome_region_length );
   next if !$genome_region_length;
   my $comp_q = $id;
   $comp_q =~ s/\W+/_/g;
   $comp_q =~ s/__+/_/g;
   $comp_q =~ s/_$//g;
   $comp_q =~ s/^_//g;

#        my $header= "Processing reference target $scaffold with query $comp_q\n";
   my $query_file = $base_query_file . "_$comp_q";
   if ( !-s $query_file ) {
    open( QUERY, ">$query_file" ) ||die ($!);
    print QUERY ">$id\n" . &wrap_text($queries{$id}{'seq'})
      if !$coding_only;
    print QUERY ">$id\n" . &wrap_text($queries{$id}{'orf'}{'seq'})
      if $coding_only;
    close QUERY;
   }
   my $scaffold_file = "$query_file.scaffold";
   if ( !-s $scaffold_file ) {
    open( SCAFFOLD, ">$scaffold_file" ) ||die($!);
    print SCAFFOLD ">$scaffold" . ":$start-$end\n".&wrap_text($genome_region);
    close SCAFFOLD;
   }
### -D -M control memory usage!
   #"echo '$header' > $query_file.exonerate_results && $exonerate"
   my $common =
       "$exonerate "
     . ' --ryo "RYOAP_START\nRYOAP_STATS_D\tAlignment length\tidentical\tsimilar\tmismatch\tidentical%%\tsimilar%%\nRYOAP_STATS\t%s\t%et\t%ei\t%es\t%em\t%pi\t%ps\nRYOAP_CODING_QUERY\t%qi\t%qab\t%qae\n>%qi\n%qas\nRYOAP_CODING_GENOME\t%ti\t%tcb\t%tce\n>%ti\n%tcs\nRYOAP_END\n" '
     . " --targettype dna -D 8192 -M 8192 --showtargetgff y -n 1 -S n --maxintron $max_intron  ";
   my $cmd;
   if ($same_species) {
    $cmd = $coding_only
      ? " -s 2000 -Q dna -m coding2genome "    ## removed --percent 90
      : " -s 2000 -Q dna -m cdna2genome $do_annotations ";    # --percent 70
    $cmd =
" -s 100 -Q protein --exhaustive 1 -m protein2genome:bestfit --percent 20 "
      if $is_protein && !$not_exhaustive;                     #
    $cmd = " -s 100 -Q protein -m protein2genome --percent 20 "
      if $is_protein && $not_exhaustive;
    $cmd .= " --hspfilter 100 --geneseed 250 ";
   }
   else {

    # I DON'T KNOW WHAT GOOD SCORE/PERCENT VALUES ARE HERE

    $cmd =
      $coding_only
      ? " -s 200 -Q dna -m coding2genome --percent 30"
      : " -s 200 -Q dna -m cdna2genome $do_annotations ";
    $cmd =
" -s 100 -Q protein --exhaustive 1 -m protein2genome:bestfit --percent 20 "
      if $is_protein && !$not_exhaustive;
    $cmd = " -s 100 -Q protein -m protein2genome --percent 20 "
      if $is_protein && $not_exhaustive;
   }
   $cmd .= ' --refine region ' if !$no_refine && !$is_protein;
   $cmd .= " --softmasktarget y " if $softmasked;
   $cmd .= " $query_file $scaffold_file ";
   print CMD (
       $common 
      . $cmd
      . " > $query_file.exonerate_results\n"
   );
  }
 }
 close(CMD);
}

sub splitfasta() {
 my @files;
 my $file2split         = shift;
 my $outdir             = shift;
 my $how_many_in_a_file = shift;
 return
   unless $file2split && -s $file2split && $outdir && $how_many_in_a_file;
 mkdir($outdir) unless -d $outdir;
 my $filecount;
 my $seqcount = int(0);
 open( FILE, $file2split );
 my $orig_sep = $/;
 $/ = ">";
 <FILE>;

 while ( my $record = <FILE> ) {
  my @lines = split( "\n", $record );
  my $id = shift @lines;
  $id =~ /^(\S+)/;
  $id = $1 || die "Cannot find ID for a sequence in $file2split";
  chomp(@lines);
  my $seq = join( '', @lines );
  $seq =~ s/>$//;
  $seq =~ s/\s+//g;
  next unless $seq;
  $seq = &wrap_text($seq);
  $seqcount++;

  if (   !$filecount
       || $how_many_in_a_file == 1
       || $seqcount > $how_many_in_a_file )
  {
   $seqcount = int(0);
   $filecount++;
   my $outfile =
     $how_many_in_a_file == 1 ? $id : $file2split . "_" . $filecount;
   $outfile = $outdir . '/' . $outfile;
   close(OUT);
   open( OUT, ">$outfile" ) || die("Cannot open $outfile");
   push( @files, $outfile );
  }
  print OUT ">$id\n$seq\n";
 }
 close(FILE);
 close(OUT);
 $/ = $orig_sep;
 return \@files;
}

sub process_exonerate() {

 #this fixes the co-ords
 my $exonerate_result_files = shift;
 my $gene_counter           = int(0);
 my $file_counter           = int(0);
 my $all_files_number       = scalar(@$exonerate_result_files);
 print "##gff-version 3\n";
 foreach my $exonerate_file (@$exonerate_result_files) {
  if ( !-s $exonerate_file ) {
   warn "File $exonerate_file not found! Skipping...\n";
   next;
  }
  open( EXONERATE, $exonerate_file );
  $file_counter++;
  print STDERR
"Processed $file_counter/$all_files_number files. Found $gene_counter genes....\r"
    if $file_counter =~ /00$/;
  while ( my $ln = <EXONERATE> ) {
   next if $ln =~ /^\s*$/;
   if ( $ln =~ /^# --- START OF GFF DUMP/ ) {
    my ( $gene_id, $gene_start, $gene_end );
    my $exon_counter   = 0;
    my $intron_counter = 1;
    my ( $mRNA, $mRNA_id, $current_score );

    while ( my $ln2 = <EXONERATE> ) {
     if ( $ln2 =~ /^# --- END OF GFF DUMP/ ) {
      print "\n";
      $gene_counter++;
      last;
     }
     elsif ( $ln2 =~ /^#/ ) {
      next;
     }
     my @data = split( "\t", $ln2 );
     if ( $data[8] ) {
      next if $data[2] eq 'similarity';
      my $offset = int(0);
      if ( $data[0] =~ s/:(\d+)-(\d+)$// ) {
       $offset = $1 - 1;
       $data[3] += $offset;
       $data[4] += $offset;
      }
      if ( $data[2] eq 'gene' ) {
       $data[8] =~ /sequence\s+(\S+)/;
       $gene_id    = $1;
       $gene_start = $data[3];
       $gene_end   = $data[4];
       undef($mRNA);    #reset
       undef($mRNA_id);

       #pasa
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
       else {
        $data[8] = "ID=$gene_id\n";
       }
      }
      elsif ( $data[2] eq 'cds' ) {
       $data[2] = 'CDS';
       $data[8] = "ID=cds.$gene_id;Parent=$gene_id\n";
      }
      elsif ( $data[2] eq 'exon' ) {
       $exon_counter++;
       $data[8] = "ID=$gene_id.exon.$exon_counter;Parent=$gene_id\n";
      }
      elsif ( $data[2] eq 'intron' ) {
       $data[8] = "ID=$gene_id.intron.$intron_counter;Parent=$gene_id\n";
      }
      else {

       # should we checking for splice site here as part of the golden checks?
       next;
      }
      $data[8] =~ s/\s+\;\s+/\;/g;
      if ( $data[2] eq 'gene' ) {
       $current_score = $data[5];
       next if ( $current_score < $exonerate_postscore );
       print join( "\t", @data );
       print $mRNA if ($mRNA);
       $gene_counter++;
      }
      else {
       next if ( $current_score < $exonerate_postscore );
       print join( "\t", @data );
      }
     }
    }

   }
  }
  close EXONERATE;
 }
 print STDERR
"Processed $file_counter/$all_files_number files. Found $gene_counter genes....\n";
}

sub process_exonerate_similarity() {

 #this fixes the co-ords
 my $exonerate_result_files = shift;
 my $gene_counter           = int(0);
 my $file_counter           = int(0);
 my $all_files_number       = scalar(@$exonerate_result_files);
 foreach my $exonerate_file (@$exonerate_result_files) {
  if ( !-s $exonerate_file ) {
   warn "File $exonerate_file not found! Skipping...\n";
   next;
  }
  open( EXONERATE, $exonerate_file );
  $file_counter++;
  print STDERR
"Processed $file_counter/$all_files_number files. Found $gene_counter genes....\r";
  while ( my $ln = <EXONERATE> ) {
   next if $ln =~ /^\s*$/;
   my ( $gene_id, $gene_start, $gene_end );
   my $exon_counter   = 0;
   my $intron_counter = 1;
   my ( $mRNA, $mRNA_id );
   my @data = split( "\t", $ln );
   if ( $data[8] ) {
    next unless $data[2] eq 'similarity';
    my $offset = int(0);
    if ( $data[0] =~ s/:(\d+)-(\d+)$// ) {
     $offset = $1 - 1;
     $data[3] += $offset;
     $data[4] += $offset;
    }
    print join( "\t", @data );
   }
  }
  close EXONERATE;
 }
 print STDERR
"Processed $file_counter/$all_files_number files. Found $gene_counter genes....\r";
 print STDERR "\n";
}

sub get_id_seq_from_fasta {
 my ( $acc, $fasta_db ) = @_;
 my $seq;
 if ( !$get_id_seq_from_fasta_hash{$fasta_db}{$acc} ) {
  $seq = `$cdbyank_exec $fasta_db.cidx -a '$acc'`;
  chomp($seq);
  $get_id_seq_from_fasta_hash{$fasta_db}{$acc} = $seq;
 }
 else {
  $seq = $get_id_seq_from_fasta_hash{$fasta_db}{$acc};
 }
 unless ($seq) {
  warn "FATAL: couldn't retrieve seq for $acc from $fasta_db\n";
  return;
 }
 my @x = split( /\n/, $seq );
 my $id = shift @x;
 $seq = join( "", @x );
 $seq = uc($seq);
 $seq =~ s/\s//g;
 return $seq;
}

sub check_program() {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  die "Error, path to required $prog cannot be found\n"
    unless $path =~ /^\//;
  chomp($path);
  $path = readlink($path) if -l $path;
  push( @paths, $path );
 }
 return @paths;
}

sub process_cmd {
 my ($cmd) = @_;
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  die "Error, cmd died with ret $ret\n";
 }
 return $ret;
}

sub convert_fasta_to_single_line() {
 my $file = shift || die;
 open( IN, $file ) || die;
 open( OUT, ">$file.t" );
 my $orig_sep = $/;
 $/ = ">";
 <IN>;
 while ( my $record = <IN> ) {
  my @lines = split( "\n", $record );
  my $id    = shift @lines;
  my $seq   = join( '', @lines );
  $seq =~ s/\s+//g;
  $seq =~ s/>$//;
  next unless $seq;
  print OUT ">$id\n$seq\n";
 }
 close OUT;
 close IN;
 unlink($file);
 rename( "$file.t", $file );
}

sub wrap_text() {
 my $string      = shift;
 my $wrap_length = shift;
 $wrap_length = 120 if !$wrap_length;
 $string =~ s/(.{0,$wrap_length})/$1\n/g;
 $string =~ s/\n{2,}/\n/;
 return $string;
}


sub read_fasta(){
 my $fasta = shift;
 my %hash;
 my $orig_sep = $/; 
 $/ = '>';
 open (IN,$fasta) || die($!);
 while (my $record = <IN>){
  chomp($record);
  next unless $record; 
  my @lines = split("\n",$record);
  my $id = shift (@lines);
  my $seq = join('',@lines);
  $seq=~s/\s+//g;
  if ($id && $seq && $id=~/^(\S+)/){
   $hash{$1} = $seq;
  }
 }
 close IN;
 $/ = $orig_sep;
 return \%hash;
}
