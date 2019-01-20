#!/usr/bin/env perl

# TODO : reverse transcriptase and transposons from uniref
# change to refseq_insect?

=head1 NAME

prepare_domain_exon_annotation.pl

=head1 USAGE

Mandatory

 -fasta|genome|in :s   => FASTA file of genome. If the file GENOMEFASTA.masked is present, then RepeatMasker is not run
 -engine          :s   => How to run hhblits: none, local, localmpi, PBS or cluster (def. local)
 -transposon_db   :s   => HHblits transposon database (provided)
 -uniprot_db      :s   => HHblits Uniprot database (see ftp://toolkit.genzentrum.lmu.de/pub/HH-suite/databases/hhsuite_dbs)
 -hosts           :s   => Only for -engine mpi: a definition for which hosts to use in the format hostname1:number_of_cpus-hostname2:number_of_cpus, e.g. localhost:5-remote:5
 -repeat_taxon    :s   => A species category for RepeatMasker such as insecta, primates, plants, vertebrates etc 

Optional

 -minsize         :i   => Minimum number of nucleotides without a stop codon to define an exon (def. 100bp)
 -circular             => If genome is a circular molecule (bacteria, mtDNA etc)
 -repeatoptions        => Any options to pass on to /all/ of the RepeatMasker runs using key=value notation. Pass multiple options delimited by :colon: (e.g. -repeatoptions nopost:frag=1000000)
 
 -rep_cpus        :i   => Number of CPUs to use for /each/ of the 5 Repeatmasking steps in parallel (def. 2, i.e. 10 CPUs across five steps)
 -mpi_cpus        :i   => Number of MPI threads (or CPUs for local and nodes for cluster) to use (if applicable; def to 2). Careful of memory usage if local!
 -scratch         :s   => If engine is local or MPI, a 'local' scratch directory to copy databases to each node, e.g. /dev/shm if there is enough space
 -no_uniprot           => Don't search for Uniprot hits. Useful if you want to conduct the transposon search separately from the Uniprot (e.g. different engine/computing environment)
 -no_transposon        => Don't search for transposon hits. See above.
 -min_exons_hints :i   => Minimum number of hits for a scaffold with same Uniprot ID before flagging sequence having a valid hit (def 2). Really useful for small scaffolds or Augustus training. Not useful otherwise but 2 is a good minimum

 -help                 => This help text and some more info
 -verbose              => Print out every command before it is executed.

 -only_repeat	       => Only do RepeatMasking and then exit

 -only_parse      :s   => Just parse this HHblits output and then exit. Useful if ran outside this script and want to produce .hints/gff files.
 
=head1 NOTES
            
Requires EMBOSS tools (getorf), RepeatMasker and HHblits (installed and in path).

-hhblits_cpus can be as high as your I/O would allow. E.g. for MPI 40 or even 200... 

For local, it shouldn't be more than the number of CPUs on your local box. for cluster, it should be the number of nodes you'd like to use

If a ffindex-powered HHblits run stops prematurely, you can restart by using the original input file and specifying transposon or uniprot, e.g.:

 ffindex_gather.sh masked.exons.aa.trim transposon 
 mv masked.exons.aa.trim.db.idx masked.exons.aa.trim.db.idx.orig
 ln -s masked.exons.aa.trim.db.idx.notdone masked.exons.aa.trim.db.idx

Now the index will tell ffindex to only process what hasn't been done. Repeat the above procedure to capture the final output to a single file

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use List::Util 'shuffle';
use POSIX qw(ceil);
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
use Fasta_reader;
use Thread_helper;
use Data::Dumper;
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/:$RealBin:$RealBin/../3rd_party/RepeatMasker:$RealBin/../3rd_party/hhsuite/bin";
$ENV{HHLIB} =  "$RealBin/../3rd_party/hhsuite";

my (
     $genome,          $circular,          
     $help,            $verbose,$only_repeat,
     $scratch_dir,     $no_uniprot_search,
     $no_transposon_search, $only_parse
);
my $repeat_taxon;
my $mpi_host_string = '';
my $repeatmasker_options = '';
my $minsize       = 100;
my $repeatcpus          = 2;
my $hhblits_cpus  = 10;
my $engine        = 'local';
my $transposon_db = $RealBin. "/../databases/hhblits/transposons";
my $uniprot_db =  $RealBin . "/../databases/hhblits/refseq_plant";
my $min_exons_before_reporting = 2;
my $gc_cutoff = 0.4;
my ( $getorf_exec, $repeatmasker_exec, $bedtools_exec ) = &check_program( 'getorf', 'RepeatMasker', 'bedtools' );
my ( $hhblits_exec, $ffindex_apply_exec, $ffindex_from_fasta_exec, $segmasker_exec, $ffindex_get_exec ) = &check_program( 'hhblits', 'ffindex_apply', 'ffindex_from_fasta', 'segmasker', 'ffindex_get' ); 


GetOptions(
            'fasta|genome|in:s' => \$genome,
            'minsize:i'         => \$minsize,
            'circular'          => \$circular,
            'rep_cpus:i'      => \$repeatcpus,
            'repeatoptions:s'   => \$repeatmasker_options,
            'engine:s'          => \$engine,
            'mpi_cpus:i'        => \$hhblits_cpus,
            'transposon_db:s'   => \$transposon_db,
            'uniprot_db:s'      => \$uniprot_db,
            'hosts:s'           => \$mpi_host_string,
            'verbose|debug'     => \$verbose,
            'help'              => \$help,
            'scratch:s'         => \$scratch_dir,
            'no_uniprot'        => \$no_uniprot_search,
            'no_transposon'     => \$no_transposon_search,
            'only_repeat'       => \$only_repeat,
            'only_parse:s'      => \$only_parse,
            'min_exons_hints:i' => \$min_exons_before_reporting,
	    'gc_cutoff:f'       => \$gc_cutoff,
	    'repeat_taxon:s'    => \$repeat_taxon
);

pod2usage( -verbose => 2 ) if $help;

die pod2usage "No genome FASTA provided\n" unless $genome && -s $genome;

die "-genome must be a full path\n" unless $genome=~/^\//;

unless ($only_repeat){
	die "The Transposon DB (eg. -transposon_db or $transposon_db) has not been prepared. See ".$RealBin . "/../databases/hhblits/README\n" unless $no_transposon_search || -s $transposon_db."_hhm_db" || -s $transposon_db."_hhm.ffdata";
	die "The UniProt DB (e.g. -uniprot_db or $uniprot_db) has not been prepared. See ".$RealBin . "/../databases/hhblits/README\n" unless $no_uniprot_search || -s $uniprot_db."_hhm_db" || -s $uniprot_db."_hhm.ffdata";
}

if ($only_parse){
 die "Cannot find HHBlits output file $only_parse\n" unless -s $only_parse;
 my $parse_results = &parse_hhr(  $only_parse, 70, 1e-3, 1e-6, 100, 50, 30 );
 print "Completed, see $parse_results\n";
 exit(0);
}
mkdir($scratch_dir) if $scratch_dir && !-s $scratch_dir && !-d $scratch_dir;
$ENV{'TMPDIR'} = $scratch_dir if $scratch_dir && -d $scratch_dir;

die "-mpi_cpus needs to be larger than 1 (not $hhblits_cpus)\n" if $hhblits_cpus < 2;
$engine = lc($engine);
die pod2usage "Engine must be local, localmpi or PBS\n"
  unless (    $engine =~ /local/
           || $engine =~ /mpi/
           || $engine =~ /pbs/
           || $engine =~ /none/
           || $engine =~ /cluster/ );

die pod2usage "-repeat_taxon is mandatory\n\n" unless $repeat_taxon;

my $genome_name = basename($genome);
$genome_name=~s/[^\w\.\-]+//g;
print "Using $genome_name as the name for $genome\n";
&check_for_iupac_violation($genome);


if (-s $genome.'.hardmasked' && !-s $genome . '.masked'){
	print "Found $genome.hardmasked. Using it as a masked file\n";
	symlink($genome.'.hardmasked',$genome . '.masked');
}

if (-s $genome . '.masked'){
  print "Found masked file $genome.masked. Skipping repeatmasking\n";
}else{
  &do_repeat_masking($repeatmasker_options);
}

$genome .= '.masked';
die "Could not find masked genome $genome.\n" unless -s $genome;
if ($only_repeat){
	print "User stop requested after RepeatMasking step.\n";
	exit;
}


my $exons = basename($genome).".exons";
my $getorf_options .= $circular ? '-circular' : '';

&process_cmd(
       "$getorf_exec -sequence $genome -outseq $exons.aa -minsize 200 -find 0 ") # nucleotide length
  unless -s $exons . '.aa';
#&process_cmd(       "$getorf_exec -sequence $genome -outseq $exons.nt -minsize 200 -find 2 ")  unless -s $exons . '.nt';

unless ( -s $exons . '.aa.trim'  ) {
 print "Post-processing...\n";
 &process_cmd("$segmasker_exec -outfmt fasta -in $exons.aa -out $exons.aa.seg") unless -s "$exons.aa.seg";
 my $hash_to_keep = &trim_X("$exons.aa.seg");
 rename("$exons.aa.seg.trim","$exons.aa.trim");
 #&trim_id( "$exons.nt", $hash_to_keep );
}
if ( $engine =~ /none/ ) {
 print "Engine is $engine. No HHblits related files prepared. Exiting\n";
 exit(0);
}


print "Preparing HHblits files for $engine\n";
# change the following so it happens per database
if ( $engine =~ /pbs/ ) {
 &prepare_pbs();
}
elsif ( $engine =~ /localmpi/ ) {
 &prepare_localmpi();
}
elsif ( $engine =~ /mpi/ ) {
 &prepare_mpi($mpi_host_string);
}
elsif ( $engine =~ /cluster/ ) {
 &prepare_cluster();
}
else {
 &prepare_local();
}

###############################################################
sub mpi_version() {
 my $exec  = shift;
 my @check = `$exec --help 2>&1`;
 my $version;
 foreach my $ln (@check) {
  if    ( $ln =~ /mpiexec/ ) { $version = 'mpich2';  last; }
  elsif ( $ln =~ /mpirun/ )  { $version = 'openmpi'; last; }
 }
 die "Can't tell which MPI version you are using, MPICH or OpenMPI"
   if !$version;
 print "Found MPI $version\n";
 $ENV{LD_LIBRARY_PATH} .= ":$RealBin/../lib64";
 return $version;
}

sub check_for_mpd($$) {
 my $option  = shift;
 my $version = shift;
 if ( $version eq 'mpich2' ) {
  if ( $option =~ /^\d+$/ ) {
   &process_cmd("mpdboot --ncpus=$option");
  }
  elsif ( -s $option ) {
   my @nodes = `cut -f 1 -d ':' $option|sort -u`;
   chomp(@nodes);
   my $nodes_size = scalar(@nodes);
   &process_cmd("mpdboot --file=$option -n $nodes_size");
  }
  else {
   die "I don't know what the MPD option $option is...\n";
  }
 }
 else {

  # nothing needs to be done for openmpi?
 }
}

sub check_program() {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  die "Error, path to required $prog cannot be found in your path\nPATH is:\n".$ENV{PATH}."\n"
    unless $path =~ /^\//;
  chomp($path);
  push( @paths, $path );
 }
 return @paths;
}

sub process_cmd {
 my ($cmd) = @_;
 print "CMD: $cmd\n" if $verbose;
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  die "Error, cmd died with ret $ret\n";
 }
 return $ret;
}

sub trim_X($) {
 my $file     = shift;
 my $fcounter = int(0);
 my $scounter = int(0);
 my %hash;
 my $minaa    = int( $minsize / 3 );
 my $orig_sep = $/;
 $/ = '>';
 open( IN,   $file );
 open( TRIM, '>' . $file . '.trim' );
 my $discard = <IN>;

 while ( my $record = <IN> ) {
  $scounter++;
  chomp($record);
  my @lines = split( "\n", $record );
  my $id    = shift @lines;
  my $seq   = join( '', @lines );
  $seq =~ s/\s*//g;
  next if length($seq) < $minaa;
  $seq =~ s/[a-z]/X/g;
  my $seq_temp = $seq;
  my $xs = ( $seq_temp =~ tr/X// );

  if ($xs) {
   next if ( $xs / length($seq) ) > 0.3;
   next if length($seq_temp) < ( $minaa * 0.9 );
  }
  print TRIM ">$id\n$seq\n";
  $fcounter++;
  $hash{$id}++;
 }
 close IN;
 close TRIM;
 $/ = $orig_sep;
 print "$file: Found $scounter exons and $fcounter passing criteria\n";
 return \%hash;
}

sub trim_id() {
 my $file     = shift;
 my $hash_ref = shift;

 my $orig_sep = $/;
 $/ = '>';
 open( IN,   $file );
 open( TRIM, '>' . $file . '.trim' );
 my $discard = <IN>;
 while ( my $record = <IN> ) {
  chomp($record);
  my @lines = split( "\n", $record );
  my $id = shift @lines;
  next unless $hash_ref->{$id};
  my $seq = join( '', @lines );
  $seq =~ s/\s*//g;
  print TRIM ">$id\n$seq\n";
 }
 close IN;
 close TRIM;
 $/ = $orig_sep;

}

sub prepare_pbs() {
 my ( $mpirun_exec, $ffindex_apply_mpi_exec ) =
   &check_program( 'mpirun', 'ffindex_apply_mpi' );

 unless ($no_transposon_search) {

  open( SCRIPT, ">hhblits_mpi_transposon.sh" );
  print SCRIPT "#!/bin/bash
NUMBERSPROCESSES=$hhblits_cpus
PROTEIN_FILE=$exons.aa.trim

if [ ! -e \$PROTEIN_FILE.db ]; then
 $ffindex_from_fasta_exec -s \$PROTEIN_FILE.db \$PROTEIN_FILE.db.idx \$PROTEIN_FILE
 mv \$PROTEIN_FILE.db.idx \$PROTEIN_FILE.db.idx.orig ; cp \$PROTEIN_FILE.db.idx.orig \$PROTEIN_FILE.db.idx.orig.notdone; ln -s \$PROTEIN_FILE.db.idx.orig.notdone \$PROTEIN_FILE.db.idx
fi

qsub -l select=\$NUMBERSPROCESSES:ncpus=1:mpiprocs=1:mem=4gb:NodeType=any -l walltime=12:00:00 -V -r n -N hbtransposons -- \$PWD/hhblits_mpi_transposon.pbs \$PROTEIN_FILE.db \$NUMBERSPROCESSES
";
  close SCRIPT;
  open( SCRIPT, ">hhblits_mpi_transposon.pbs" );
  print SCRIPT "#!/bin/bash
MPIRUN_EXEC=$mpirun_exec
MPIRUN_ARGS=\"-gmca mpi_warn_on_fork 0 -cpus-per-proc 1 -np \$2 -machinefile workers.\$PBS_JOBID.mpi\"
DB=$transposon_db

export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$RealBin/../lib64
cd \$PBS_O_WORKDIR
cat \${PBS_NODEFILE} > workers.\$PBS_JOBID.mpi
\$MPIRUN_EXEC \$MPIRUN_ARGS $ffindex_apply_mpi_exec \\
  -d \"\$1\".transposon.db \\
  -i \"\$1\".transposon.db.idx \\
  \$1 \\
  \$1.idx \\
  -- $hhblits_exec -maxmem 3 -d \$DB -n 1 -mact 0.5 -cpu 1 -i stdin -o stdout -e 1e-5 -E 1E-05 -id 80 -p 80 -z 0 -b 0 -v 0 -B 3 -Z 3 2>/dev/null
";
  close SCRIPT;

 }
 unless ($no_uniprot_search) {
  open( SCRIPT, ">hhblits_mpi_uniprot.sh" );
  print SCRIPT "#!/bin/bash
NUMBERSPROCESSES=$hhblits_cpus
PROTEIN_FILE=$exons.aa.trim

if [ ! -e \$PROTEIN_FILE.db ]; then
 $ffindex_from_fasta_exec -s \$PROTEIN_FILE.db \$PROTEIN_FILE.db.idx \$PROTEIN_FILE
 mv \$PROTEIN_FILE.db.idx \$PROTEIN_FILE.db.idx.orig ; cp \$PROTEIN_FILE.db.idx.orig \$PROTEIN_FILE.db.idx.orig.notdone; ln -s \$PROTEIN_FILE.db.idx.orig.notdone \$PROTEIN_FILE.db.idx
fi

qsub -l select=\$NUMBERSPROCESSES:ncpus=1:mpiprocs=1:mem=7gb:NodeType=any -l walltime=12:00:00 -V -r n -N hbuniprot -- \$PWD/hhblits_mpi_uniprot.pbs \$PROTEIN_FILE.db \$NUMBERSPROCESSES
";
  close SCRIPT;

  open( SCRIPT, ">hhblits_mpi_uniprot.pbs" );
  print SCRIPT "#!/bin/bash
MPIRUN_EXEC=$mpirun_exec
MPIRUN_ARGS=\"-gmca mpi_warn_on_fork 0 -cpus-per-proc 1 -np \$2 -machinefile workers.\$PBS_JOBID.mpi\"
DB=$uniprot_db

export OMP_NUM_THREADS=1
export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$RealBin/../lib64
cd \$PBS_O_WORKDIR
cat \${PBS_NODEFILE} > workers.\$PBS_JOBID.mpi
\$MPIRUN_EXEC \$MPIRUN_ARGS $ffindex_apply_mpi_exec \\
  -d \"\$1\".uniprot.db \\
  -i \"\$1\".uniprot.db.idx \\
  \$1 \\
  \$1.idx \\
  -- $hhblits_exec -maxmem 5 -d \$DB -n 1 -mact 0.5 -cpu 1 -i stdin -o stdout -e 1e-5 -E 1E-05 -id 80 -p 80 -z 0 -b 0 -v 0 -B 3 -Z 3 2>/dev/null
";
  close SCRIPT;
 }
 chmod 0755, qw/hhblits_mpi_transposon.sh hhblits_mpi_transposon.pbs hhblits_mpi_uniprot.sh hhblits_mpi_uniprot.pbs/;
 print
"Wrote PBS scripts. See scripts hhblits_mpi* and $exons.aa.trim and transfer them to your PBS batch system. You will almost certainly need to edit the pbs/.sh scripts to make them work for your setup (e.g. your account string or load any modules)\n";

}

sub prepare_localmpi() {
 my ( $mpirun_exec, $ffindex_apply_mpi_exec ) =
   &check_program( 'mpirun', 'ffindex_apply_mpi' );
 my $mpi_version = &mpi_version($mpirun_exec);
 &check_for_mpd( $hhblits_cpus, $mpi_version );
 unless ($no_transposon_search || -s "hhr.$exons.aa.trim.db.transposon.db") {
  print "Transposon database...\n";
  &process_cmd("$ffindex_from_fasta_exec -s $exons.aa.trim.db $exons.aa.trim.db.idx $exons.aa.trim") if !-s "$exons.aa.trim.db";
  my $number_of_entries = `wc -l < $exons.aa.trim.db.idx`;
  chomp($number_of_entries);
  unlink("mpi_err.log");
  print "Processing $number_of_entries entries with $hhblits_cpus threads...\n";
  my $transposon_cmd ="$mpirun_exec -n $hhblits_cpus $ffindex_apply_mpi_exec -d $exons.aa.trim.db.transposon.db -i $exons.aa.trim.db.transposon.db.idx $exons.aa.trim.db $exons.aa.trim.db.idx  \\
 -- $hhblits_exec -maxmem 3 -d $transposon_db -n 1 -mact 0.5 -cpu 1 -i stdin -o stdout -e 1E-5 -E 1E-5 -id 80 -p 80 -z 0 -b 0 -B 3 -Z 3 -v 0 2>mpi_err.log
 ";
  &process_cmd($transposon_cmd) unless $number_of_entries == 0 || -s "$exons.aa.trim.db.transposon.db" ;
 }
 my $transposon_results = &parse_hhr( "$exons.aa.trim.db.transposon.db", 70, 1e-3, 1e-6, 100, 50, 30, 'yes' );
 my $noreps_fasta = &remove_transposons( "$exons.aa.trim", $transposon_results );

 unless ($no_uniprot_search || -s "hhr.$noreps_fasta.uniprot.db"){
  print "Uniprot database...\n";
  &process_cmd("$ffindex_from_fasta_exec -s $noreps_fasta.db $noreps_fasta.db.idx $noreps_fasta"
  ) if !-s "$noreps_fasta.db";
  my $number_of_entries = `wc -l < $noreps_fasta.db.idx`;
  chomp($number_of_entries);
  unlink("mpi_err.log");
  print
"Processing $number_of_entries entries with $hhblits_cpus threads for uniprot...\n";
  my $uniprot_cmd =
"$mpirun_exec -n $hhblits_cpus $ffindex_apply_mpi_exec -d $noreps_fasta.db.uniprot.db -i $noreps_fasta.db.uniprot.db.idx $noreps_fasta.db $noreps_fasta.db.idx  \\
 -- $hhblits_exec -maxmem 5 -d $uniprot_db -n 1 -mact 0.5 -cpu 1 -i stdin -o stdout -e 1E-5 -E 1E-5 -id 80 -p 80 -z 0 -b 0 -B 3 -Z 3 -v 0 2>mpi_err.log
 ";

  &process_cmd($uniprot_cmd)  unless $number_of_entries == 0 || -s "$noreps_fasta.db.uniprot.db";
  &parse_hhr( "$noreps_fasta.db.uniprot.db", 70, 1e-3, 1e-6, 100, 50, 30 ) unless  -s "hhr.$noreps_fasta.uniprot.db";
 }
}

sub prepare_mpi() {
 my $host_string = shift;
 pod2usage "For MPI engine I need a host definition with -hosts\n" if $engine && $engine eq 'mpi' && !$host_string;
 my %hash;
 my @hosts_defs = split( '-', $host_string );
 foreach my $host_str (@hosts_defs) {
  my ( $host, $nodes ) = split( ':', $host_str );
  die "Invalid MPI host definition\n" unless $nodes && $nodes =~ /^\d+$/;
  my @check = `ssh $host -t echo testing $host 2>&1`;
  die "Cannot connect to host $host\n" unless $check[0] =~ /^test/;
  $hash{$host} = $nodes;
 }
 my $workers_file = "machinefile.$$";
 my ( $mpirun_exec, $ffindex_apply_mpi_exec ) =
   &check_program( 'mpirun', 'ffindex_apply_mpi' );

 # which mpirun? mpich2 or openmpi?
 my $mpi_version = &mpi_version($mpirun_exec);

 my $cpu_count = int(0);
 open( OUT, ">$workers_file" );
 foreach my $worker (sort keys %hash ) {
  print OUT "$worker\n" x $hash{$worker} if $mpi_version eq 'openmpi';
  print OUT "$worker:" . $hash{$worker} . "\n" if $mpi_version eq 'mpich2';
  $cpu_count += $hash{$worker};
 }
 close OUT;
 die
   "CPU count $cpu_count from -hosts is not equal to -mpi_cpu $hhblits_cpus\n"
   unless $cpu_count == $hhblits_cpus;

 &check_for_mpd( $workers_file, $mpi_version );

 if ($scratch_dir) {
  foreach my $worker (sort keys %hash ) {
   unless ($no_transposon_search) {
    print "Copying $transposon_db to $worker scratch\n";
    &process_cmd("rsync -lua --perms $transposon_db* $worker:$scratch_dir/");
   }
   unless ($no_uniprot_search) {
    print "Copying $uniprot_db to $worker scratch\n";
    &process_cmd("rsync -lua --perms $uniprot_db* $worker:$scratch_dir/");
   }
  }
  $transposon_db = $scratch_dir . '/' . basename($transposon_db);
  $uniprot_db    = $scratch_dir . '/' . basename($uniprot_db);
 }
 unless ($no_transposon_search) {
  print "Transposon database...\n";
  &process_cmd(
"$ffindex_from_fasta_exec -s $exons.aa.trim.db $exons.aa.trim.db.idx $exons.aa.trim"
  ) if !-s "$exons.aa.trim.db";
  my $number_of_entries = `wc -l < $exons.aa.trim.db.idx`;
  chomp($number_of_entries);
  unlink("mpi_errors.log");
  print "Processing $number_of_entries entries with $hhblits_cpus threads for transposons...\n";
  &process_cmd("$mpirun_exec -machinefile $workers_file -n $hhblits_cpus $ffindex_apply_mpi_exec -d $exons.aa.trim.db.transposon.db -i $exons.aa.trim.db.transposon.db.idx $exons.aa.trim.db $exons.aa.trim.db.idx  \\
 -- $hhblits_exec -maxmem 3 -d $transposon_db -n 1 -mact 0.5 -cpu 1 -i stdin -o stdout -e 1E-5 -E 1E-5 -id 80 -p 80 -z 0 -b 0 -B 3 -Z 3 -v 0 2>> mpi_errors.log
 ") unless $number_of_entries == 0 || -s "$exons.aa.trim.db.transposon.db" ;
 }
 my $transposon_results = &parse_hhr( "$exons.aa.trim.db.transposon.db", 70, 1e-3, 1e-6, 100, 50, 30, 'yes' ) unless -s "hhr.$exons.aa.trim.transposon.db";
 
 my $noreps_fasta = &remove_transposons( "$exons.aa.trim", $transposon_results );

 unless ($no_uniprot_search) {
  print "Uniprot database...\n";
  &process_cmd(
"$ffindex_from_fasta_exec -s $noreps_fasta.db $noreps_fasta.db.idx $noreps_fasta"
  ) if !-s "$noreps_fasta.db";
  my $number_of_entries = `wc -l < $noreps_fasta.db.idx`;
  chomp($number_of_entries);
  unlink("mpi_errors.log");
  print "Processing $number_of_entries entries with $hhblits_cpus threads for uniprot...\n";

  &process_cmd(
"$mpirun_exec -machinefile $workers_file -n $hhblits_cpus $ffindex_apply_mpi_exec -d $noreps_fasta.db.uniprot.db -i $noreps_fasta.db.uniprot.db.idx $noreps_fasta.db $noreps_fasta.db.idx  \\
 -- $hhblits_exec -maxmem 5 -d $uniprot_db -n 1 -mact 0.5 -cpu 1 -i stdin -o stdout -e 1E-5 -E 1E-5 -id 80 -p 80 -z 0 -b 0 -B 3 -Z 3 -v 0 2>> mpi_errors.log
 "
  ) unless $number_of_entries == 0 || -s "$noreps_fasta.db.uniprot.db";
  &parse_hhr( "$noreps_fasta.db.uniprot.db", 70, 1e-3, 1e-6, 100, 50, 30 ) unless -s "hhr.$noreps_fasta.uniprot.db";
 }
 
 
 if ($scratch_dir ) {
  foreach my $worker (sort keys %hash ) {
   print "Removing $transposon_db from $worker scratch\n";
   &process_cmd("rm -f $worker:$transposon_db*");
   unless ($no_uniprot_search) {
    print "Removing $uniprot_db from $worker scratch\n";
    &process_cmd("rm -f  $worker:$uniprot_db*");
   }
  }
 }
}

sub remove_transposons() {
 my $fasta_file  = shift;
 my $result_file = shift;
 my $out_fasta   = "$fasta_file.norep";
 return $out_fasta if -s $out_fasta; 
 return $fasta_file if !$result_file || !-s $result_file;
 return $out_fasta if -s $out_fasta; 
 my %transposon_hits;
 open( IN, $result_file );
 while ( my $ln = <IN> ) {
  if ( $ln =~ /^(\S+)/ ) {
   $transposon_hits{$1} = 1;
  }
 }
 close IN;

 open( FASTA, $fasta_file );
 open( OUT,   ">$out_fasta" );
 my $orig_sep = $/;
 $/ = '>';
 my $disc = <FASTA>;
 while ( my $record = <FASTA> ) {
  chomp($record);
  my @lines = split( "\n", $record );
  my $id = shift(@lines);
  $id =~ /^(\S+)/;
  my $lid = $1;
  next if $transposon_hits{$1};
  print OUT ">$id\n" . join( "\n", @lines ) . "\n";
 }
 close FASTA;
 close OUT;
 $/ = $orig_sep;
 return $out_fasta;
}

sub remove_zero_bytes() {
 my $infile    = shift;
 my ($name, $path, $suffix) = fileparse($infile);
 my $outfile = $path ."hhr.$name$suffix";
 return $outfile if (-s $outfile);
 &process_cmd("cat $infile* | tr -d '\\000' > $outfile");
 system("rm -f $infile*");
 return $outfile;
}

sub prepare_local() {
 &cleanup_threaded_exit();
 my $fasta = "$exons.aa.trim";
 if ($scratch_dir) {
   unless ($no_transposon_search) {
    print "Copying $transposon_db to $scratch_dir\n";
    &process_cmd("rsync -lua --perms $transposon_db* $scratch_dir/");
   }
   unless ($no_uniprot_search) {
    print "Copying $uniprot_db to $scratch_dir\n";
    &process_cmd("rsync -lua --perms $uniprot_db* $scratch_dir/");
   }
  $transposon_db = $scratch_dir . '/' . basename($transposon_db);
  $uniprot_db    = $scratch_dir . '/' . basename($uniprot_db);
 }

 &process_cmd("$ffindex_from_fasta_exec -s $fasta.db $fasta.db.idx $fasta")  unless -s "$fasta.db";
 my $number_of_entries = `wc -l < $fasta.db.idx`;
 unlink("local_errors.log");
 chomp($number_of_entries);
 # transposons
 unless ( -s "$fasta.hhblits.transposon.cmds" || -s "hhr.$exons.aa.trim.transposon.db" || $no_transposon_search || $number_of_entries == 0 ) {
   open( CMD, ">$fasta.hhblits.transposon.cmds" );
   for ( my $i = 1 ; $i <= $number_of_entries ; $i++ ) {
    print CMD "$ffindex_get_exec -n $fasta.db $fasta.db.idx $i | $hhblits_exec -maxmem 3 -d $transposon_db -n 1 -mact 0.5 -cpu 1 -i stdin -o stdout -e 1E-5 -E 1E-5 -id 80 -p 80 -z 0 -b 0 -B 3 -Z 3 -v 0  >> $fasta.transposons.hhr 2>> local_errors.log\n"; 
   }
   close CMD;
  }

 unless ($no_transposon_search || !-s "$fasta.hhblits.transposon.cmds" ) {
   print "Processing transposon library\n";
   &just_run_my_commands("$fasta.hhblits.transposon.cmds");
 }

  my $transposon_results = &parse_hhr( "$fasta.transposons.hhr", 70, 1e-3, 1e-6, 100, 50, 30, 'yes' );
  $fasta = &remove_transposons( $fasta, $transposon_results );
  &process_cmd("$ffindex_from_fasta_exec -s $fasta.db $fasta.db.idx $fasta")  unless -s "$fasta.db";
  $number_of_entries = `wc -l < $fasta.db.idx`;
  chomp($number_of_entries);
  unlink("local_errors.log");

  # uniprot
  unless ( -s "$fasta.hhblits.uniprot.cmds" || -s "hhr.$exons.aa.trim.uniprot.db" || $no_uniprot_search || $number_of_entries == 0) {
   open( CMD, ">$fasta.hhblits.uniprot.cmds" );
   for ( my $i = 1 ; $i <= $number_of_entries ; $i++ ) {
    print CMD "$ffindex_get_exec -n $fasta.db $fasta.db.idx $i | $hhblits_exec -maxmem 5 -d $uniprot_db -n 1 -mact 0.5 -cpu 1 -i stdin -o stdout -e 1E-5 -E 1E-5 -id 80 -p 80 -z 0 -b 0 -B 3 -Z 3 -v 0 >> $fasta.uniprot.hhr 2>> local_errors.log\n";
   }
   close CMD;
  }

 unless ($no_uniprot_search || !-s "$fasta.hhblits.uniprot.cmds") {
  print "Processing UniProt library\n";
  &just_run_my_commands("$fasta.hhblits.uniprot.cmds");
 }
 &parse_hhr( "$fasta.uniprot.hhr", 70, 1e-3, 1e-6, 100, 50, 30 );
}

sub just_run_my_commands_helper(){
	my ($cmd,$failed_filehandle,$completed_filehandle) = @_;
	chomp($cmd);
	print "CMD: $cmd\n" if $verbose;
	my $ret = system($cmd);
	if ($ret && $ret != 256 ){
		print $failed_filehandle $cmd."\n" if $failed_filehandle;
	}else{
		print $completed_filehandle $cmd."\n" if $completed_filehandle;
	}
}

sub just_run_my_commands(){
 my $cmd_file = shift;
 return unless $cmd_file && -s $cmd_file;
 my (%completed,%failed);
 my $number_commands = int(0);
 my $cmd_count = int(0);
 my $failed_count = int(0);
 my $completed_count = int(0);
 if (-s $cmd_file.'.completed'){
	open (IN,$cmd_file.'.completed');
	while (my $cmd=<IN>){
		next if $cmd=~/^\s*$/;
		$completed{$cmd}++;
	}
	close IN;
 }
 open (CMDS,$cmd_file);
 while (my $cmd=<CMDS>){
	next if $cmd=~/^\s*$/ || $completed{$cmd};
	$number_commands++;
 }
 close CMDS;
 return unless $number_commands && $number_commands > 0;
 print "Processing with $hhblits_cpus CPUs...\n";
 my $thread_helper = new Thread_helper($hhblits_cpus);
 open (CMDS,$cmd_file);
 open (my $failed_fh,">$cmd_file.failed");
 open (my $completed_fh,">>$cmd_file.completed");
 while (my $cmd=<CMDS>){
	next if $cmd=~/^\s*$/ || $completed{$cmd};
	$cmd_count++;
        my $thread = threads->create('just_run_my_commands_helper', $cmd,$failed_fh,$completed_fh);
        $thread_helper->add_thread($thread);
        sleep(1) if ($cmd_count % 100 == 0);
        print "\r $cmd_count / $number_commands                 " if $verbose;
 }
 close CMDS;
 $thread_helper->wait_for_all_threads_to_complete();
 close $failed_fh;
 close $completed_fh;
 my @failed_threads = $thread_helper->get_failed_threads();
 if (@failed_threads || -s "$cmd_file.failed") {
  die "Error, " . scalar(@failed_threads) . " threads failed. Also see $cmd_file.failed\n";
  exit(1);
 }
}

sub prepare_cluster() {

# to be honest, I don't know what this should produce/ currently it produces commands like local but with a twist on number of CPUs.
 my $workdir = 'exons_hhsearch';
 if ( -d $workdir ) {
  warn
"Working directory $workdir already exists! I don't want to overwrite anything so I won't create any files for 'cluster' unless you delete it\n";
 # return;
 }
 else {
  mkdir($workdir);
 }
 my @fasta_files =
   &partition_transcript_db( "$exons.aa.trim", 'exons_hhsearch' );

 foreach my $fasta (@fasta_files) {
  &process_cmd("$ffindex_from_fasta_exec -s $fasta.db $fasta.db.idx $fasta");
  my $number_of_entries = `wc -l < $fasta.db.idx`;
  chomp($number_of_entries);
  unless ($no_transposon_search || $number_of_entries == 0) {
   open( CMD, ">$fasta.hhblits.transposon.cmds" );
   for ( my $i = 1 ; $i <= $number_of_entries ; $i++ ) {
    print CMD
"$ffindex_get_exec -n $fasta.db $fasta.db.idx $i | $hhblits_exec -maxmem 3 -d $transposon_db -n 1 -mact 0.5 -cpu 2 -i stdin -o stdout -e 1E-5 -E 1E-5 -id 80 -p 80 -z 0 -b 0 -B 3 -Z 3 -v 0 >> $fasta.transposons.hhr 2>/dev/null\n";
   }
   close CMD;
  }

  if ( -s "$fasta.transposon.hhr" ){
    my $transposon_results = &parse_hhr( "$fasta.transposon.hhr", 70, 1e-3, 1e-6, 100, 50, 30, 'yes' );
    $fasta = &remove_transposons( $fasta, $transposon_results );
    &process_cmd("$ffindex_from_fasta_exec -s $fasta.db $fasta.db.idx $fasta")  unless -s "$fasta.db";
    $number_of_entries = `wc -l < $fasta.db.idx`;
    chomp($number_of_entries);

    
    unless ($no_uniprot_search) {
      open( CMD, ">$fasta.hhblits.uniprot.cmds" );
      for ( my $i = 1 ; $i <= $number_of_entries ; $i++ ) {
	print CMD
	  "$ffindex_get_exec -n $fasta.db $fasta.db.idx $i | $hhblits_exec -maxmem 5 -d $uniprot_db -n 1 -mact 0.5 -cpu 2 -i stdin -o stdout -e 1E-5 -E 1E-5 -id 80 -p 80 -z 0 -b 0 -B 3 -Z 3 -v 0 >> $fasta.uniprot.hhr 2>/dev/null\n";
      }
      close CMD;
    }
  }
}
 print "Done, prepared for 2 CPUs per node. See $workdir/*cmds\n";
}

sub shuffle_file() {
 my $file = shift;
 open( IN, $file );
 my @lines = <IN>;
 close IN;
 open( OUT, ">$file." );
 print OUT shuffle @lines;
 close OUT;
 unlink($file);
 rename( $file . '.', $file )

}

sub partition_transcript_db {
 my $transcript_db  = shift;
 my $workdir        = shift;
 my $number_of_peps = `grep '^>' $transcript_db | wc -l `;
 chomp $number_of_peps;

 my $seqs_per_partition = ceil( $number_of_peps / ($hhblits_cpus) );
 $seqs_per_partition = 1 if $seqs_per_partition < 1;
 $seqs_per_partition = $seqs_per_partition < 5000 ? $seqs_per_partition : 5000;
 my @files;
 my $fasta_reader      = new Fasta_reader($transcript_db);
 my $partition_counter = 0;
 my $counter           = 0;
 my $ofh;

 while ( my $seq_obj = $fasta_reader->next() ) {
  my $fasta_entry = $seq_obj->get_FASTA_format();
  $fasta_entry =~ s/[\*\s]+$//;    #strip stop codon/empty space
  $fasta_entry .= "\n";
  if ( $counter % $seqs_per_partition == 0 ) {
   close $ofh if $ofh;
   $partition_counter++;
   my $outfile = "$workdir/partition.$counter.fa";
   push( @files, $outfile );
   open( $ofh, ">$outfile" ) or die "Error, cannot write to outfile: $outfile";
  }
  print $ofh $fasta_entry;
  $counter++;
 }
 close $ofh if $ofh;
 return (@files);
}

sub cleanup_threaded_exit() {
 my $process_group = $$;
 $SIG{'INT'} = sub {
  warn "Interrupt for $process_group\n";
  system("kill -9 -$process_group");
 };
 $SIG{'KILL'} = sub {
  warn "Interrupt for $process_group\n";
  system("kill -9 -$process_group");
 };
}

sub parse_hhr() {

 # (70,1e-3,1e-6,100,50,30);
 my ( $infile, $homology_prob_cut, $eval_cut, $pval_cut, $score_cut,
      $align_col_cut, $template_aln_size_cut, $is_repeat )
   = @_;
 my ($name, $path, $suffix) = fileparse($infile);
 my $outfile = $path . "hhr.$name$suffix.results";
 return $outfile if (-s $outfile);  
 return if !$infile || !-s $infile;
 $infile = &remove_zero_bytes( $infile);
 print "Post-processing $infile\n";   
 my $min_filesize = 500;
 my ( $qcounter, %reporting_hash, %hit_counter );

 die "Can't find $infile or it is too small\n" unless -s $infile && ( -s $infile ) >= $min_filesize;

 open( IN, $infile ) || die($!);

 my ($query,$scaffold_id,$reverse,$start,$stop,%scaffold_hits);

 while ( my $ln = <IN> ) {
  if ( $ln =~ /^\W*Query\s+(\S+)/ ) {
   $qcounter++;
   $query = $1;	# original name
   $scaffold_id    = $query;
   $scaffold_id =~ s/_\d+$//; # strip ORF identifier to make it a scaffold id
   $reverse = ($ln =~ /REVERSE SENSE/) ? 1 : 0;
   $ln =~ /\[(\d+)\s\-\s(\d+)\]/;
   $start = $1 && $1 =~ /^(\d+)$/ ? $1 : int(0);
   $stop  = $2 && $2 =~ /^(\d+)$/ ? $2 : int(0);
   next;
  }
  elsif ( $ln =~ /^\s*No Hit/ ) {
   while ( my $ln2 = <IN> ) {
    last if $ln2 =~ /^\s*$/;
    last if $ln2 =~ /^Done/;
    $ln2 =~ /^\s*(\d+)/;
    my $hit_number = $1;
    next unless $hit_number == 1;
    my ( $hit_desc, $hit_data, $hit_id, $strand );
    $hit_desc = substr( $ln2, 4, 31 );
    $hit_data = substr( $ln2, 35 );
    $hit_desc =~ /^(\S+)\s*(.*)/;
    $hit_id   = $1;
    $hit_desc = $2;

    if ($hit_desc) {
     $hit_desc =~ s/[\s\.]+$//;
     $hit_desc =~ s/\s+/_/g;
    }
    chomp($hit_data);
    $hit_data =~ s/^\s+//;
    my @data = split( /\s+/, $hit_data );
    my ( $prob, $evalue, $pvalue, $score, $structure_score, $alignment_length )
      = ( $data[0], $data[1], $data[2], $data[3], $data[4], $data[5] );
    $data[6] =~ /(\d+)\-(\d+)/;
    my $aa_start = $1;
    my $aa_stop  = $2;
    $data[7] =~ /(\d+)\-(\d+)/;
    my $hit_start = $1;
    my $hit_stop  = $2;

    if ( $data[7] =~ s/\((\d+)\)// ) {
     $data[8] = $1;
    }
    else {
     $data[8] =~ s/[\(\)]//g;
    }
    my $template_size     = $data[8];
    my $template_aln_size = abs( $hit_start - $hit_stop ) + 1;

    next if $homology_prob_cut > $prob;
    next if $eval_cut && $eval_cut < $evalue;
    next if $pval_cut && $pval_cut < $pvalue;
    next if $score_cut && $score_cut > $score;
    next if $alignment_length < $align_col_cut;
    next if $template_aln_size < $template_aln_size_cut;
    my ( $gff_start, $gff_end );
    if ( !$reverse ) {
     $strand = '+';
     $gff_start = $start + ( 3 * $aa_start ) - 1;
     $gff_end   = $start + ( 3 * $aa_stop ) - 1;
    }
    else {
     $strand = '-';
     $gff_start = $start - ( 3 * $aa_start ) + 1;
     $gff_end   = $start - ( 3 * $aa_stop ) + 1;
    }
    my $src  = $is_repeat ? 'RM'           : 'HU';
    my $type = $is_repeat ? 'nonexonpart' : 'CDSpart';
    my $prio = $is_repeat ? 6             : 5;
    my $uid  = "$hit_id.s$hit_start.e$hit_stop";
    $hit_counter{$uid}++;
    $uid .= '.n' . $hit_counter{$uid};

    my $name = $uid;
    $name .= " ($hit_desc)" if $hit_desc;
    my %hash = (
       'scaffold_id'=>$scaffold_id, 'type' => $type,'gff_start'=>$gff_start,'gff_end'=>$gff_end,'score'=>$score,'strand'=>$strand,'src'=>$src,'prio'=>$prio,'uid'=>$uid,'hit_start'=>$hit_start,'hit_stop'=>$hit_stop,'name'=>$name,'evalue'=>$evalue
    );
    $scaffold_hits{$scaffold_id}{$hit_id}++;
    push(@{$reporting_hash{$query}{$hit_id}},\%hash);
    last;    # top hit
   }
  }
 }
 close IN;

 open( OUT,     ">$outfile" );
 open( GLIMMER, ">$outfile.glimmer" );
 open( GENEID,  ">$outfile.geneid" );
 open( GFF3,    ">$outfile.gff3" );
 open( HINTS,   ">$outfile.hints" );
 foreach my $query_id (sort keys %reporting_hash) {
	my $scaffold_id = $query_id;
	$scaffold_id=~s/_\d+$//;
        foreach my $hit_id (keys %{$scaffold_hits{$scaffold_id}}){
                my $number_of_times = $scaffold_hits{$scaffold_id}{$hit_id};
                if ($is_repeat || $number_of_times >= $min_exons_before_reporting){
                  print OUT $query_id . "\n";
                  last;
                }
        }
  foreach my $hit_id (keys %{$reporting_hash{$query_id}}){
	my $scaffold_id = $query_id;
	$scaffold_id=~s/_\d+$//;
        my $number_of_times = $scaffold_hits{$scaffold_id}{$hit_id};
        if ($is_repeat || $number_of_times >= $min_exons_before_reporting){
             foreach my $hash_ref (@{$reporting_hash{$query_id}{$hit_id}}){
                my ($id,$type,$gff_start,$gff_end,$score,$strand,$src,$prio,$uid,$name,$hit_start,$hit_stop,$evalue) = (
       $hash_ref->{'scaffold_id'},  $hash_ref->{'type'},$hash_ref->{'gff_start'},$hash_ref->{'gff_end'},$hash_ref->{'score'},$hash_ref->{'strand'},$hash_ref->{'src'},$hash_ref->{'prio'},$hash_ref->{'uid'},$hash_ref->{'name'},$hash_ref->{'hit_start'},$hash_ref->{'hit_stop'},$hash_ref->{'evalue'}
                        );
			# In GFF3 (such as hints etc START < STOP
			my $gff3_start = $gff_start;
			my $gff3_stop = $gff_end;
				if ($gff3_start > $gff3_stop){
					my $t = $gff3_start;
					$gff3_start = $gff3_stop;
					$gff3_stop = $t;
			}

                     print HINTS "$id\thhblits\t$type\t$gff3_start\t$gff3_stop\t$score\t$strand\t.\tsrc=$src;grp=$hit_id;pri=$prio\n";
                     print GFF3  "$id\thhblits\tprotein_match\t$gff3_start\t$gff3_stop\t$score\t$strand\t.\tID=$uid;Name=$name;Target=$hit_id $hit_start $hit_stop\n";
                     print GENEID "$id\thhblits\tsr\t$gff_start\t$gff_end\t$score\t$strand\t.\n";
                     if ($strand eq '-') {
                        print GLIMMER "$id $gff_end $gff_start $score $evalue\n\n";
                     }else {
                        print GLIMMER "$id $gff_start $gff_end $score $evalue\n\n";
                     }
             }
        }
  }
 }
 close OUT;
 close GLIMMER;
 close HINTS;
 close GFF3;
 close GENEID;
 if ( -s "$outfile.gff3" ) {
  system("sort -nk 4,4 $outfile.gff3| sort -s -k 1,1 > $outfile.gff3.sorted");
  rename( "$outfile.gff3.sorted", "$outfile.gff3" );
  system("sort -nk4,4 $outfile.hints|sort -s -k1,1 > $outfile.hints. ");
  rename( "$outfile.hints.", "$outfile.hints" );
 }

 return $outfile;

}


sub do_repeat_masking(){
  # we run them separately
  # we combine to make a redundant GFF and fasta

  # these options are run in every step.
  my $repeatmasker_options = shift;
  my $frag;#=5000000;
  print "Masking repeats with $repeatcpus CPUs (cf $genome.repeatmasking.log)..\n";
  if ($repeatmasker_options){
    my @option_array = split(":",$repeatmasker_options);
    $repeatmasker_options = '';
    foreach my $option (@option_array){
      my @key_value = split("=",$option);
      $repeatmasker_options.=' -'.$key_value[0].' '.$key_value[1] if $key_value[1];
      $repeatmasker_options.=' -'.$key_value[0] if !$key_value[1];
      $frag = $key_value[1] if $key_value[0] eq 'frag';
    }
    print "Will use RepeatMasker options: $repeatmasker_options\n";
  } 
  $frag=5000000 if !$frag;
  my $cwd = `pwd`;chomp($cwd);
  mkdir("simple-only") unless -d "simple-only";
  mkdir("rna-specific") unless -d "rna-specific";
  mkdir("species-specific") unless -d "species-specific";
  mkdir("general") unless -d "general";
  mkdir("contamination-check") unless -d "contamination-check";
  my (@threads_submitted);
  my $thread_helper = new Thread_helper(5);
  my $cmd = "$repeatmasker_exec $repeatmasker_options -e ncbi -gff -pa $repeatcpus -qq -s -excln -xsmall -gccalc -frag $frag ";
  #1 "contamination-check"
	unless (-s "contamination-check/genome.fasta.cat.gz"){
		my $local_cmd = "cd $cwd/contamination-check && rm -f genome.fasta && ln -s $genome genome.fasta && $cmd -is_only -species $repeat_taxon genome.fasta 2>&1 > repeatmasking.log && cd $cwd";
  		my $thread = threads->create('just_run_my_commands_helper', $local_cmd, undef, undef);
	        $thread_helper->add_thread($thread);
		push(@threads_submitted,$thread);
		sleep(10);
	}

  #2 "simple-only"
	 unless (-s "simple-only/genome.fasta.cat.gz"){
		my $local_cmd = "cd $cwd/simple-only && rm -f genome.fasta && ln -s $genome genome.fasta && $cmd -species $repeat_taxon -no_is -noint -norna genome.fasta 2>&1 > repeatmasking.log && cd $cwd/" ;
	  	my $thread = threads->create('just_run_my_commands_helper', $local_cmd, undef, undef);
        	$thread_helper->add_thread($thread);
		push(@threads_submitted,$thread);
		sleep(10);
	}

  #3 "general"
	unless (-s "general/genome.fasta.cat.gz"){
		my $local_cmd = "cd $cwd/general && rm -f genome.fasta && ln -s $genome genome.fasta && $cmd -species $repeat_taxon -no_is -nolow genome.fasta 2>&1 > repeatmasking.log && cd $cwd";
 	 	my $thread = threads->create('just_run_my_commands_helper', $local_cmd, undef, undef);
        	$thread_helper->add_thread($thread);
		push(@threads_submitted,$thread);
		sleep(10);
	}

  #4 "rna-specific"
	unless (-s "rna-specific/genome.fasta.cat.gz"){
		my $local_cmd = "cd $cwd/rna-specific && rm -f genome.fasta && ln -s $genome genome.fasta && $cmd -no_is -nolow -lib $RealBin/../databases/repeats/rnammer-SILVA.classified.nr95.renamed.fasta genome.fasta 2>&1 > repeatmasking.log && cd $cwd";
 	 	my $thread = threads->create('just_run_my_commands_helper', $local_cmd, undef, undef);
        	$thread_helper->add_thread($thread);
		push(@threads_submitted,$thread);
 		sleep(10);
	}

  #5 RepeatModeller for "species-specific"
	unless (-s "$cwd/species-specific/$genome_name.rm.nal" || -s "$cwd/species-specific/$genome_name.rm.nin"){
	  	&process_cmd("cd $cwd/species-specific && $RealBin/../3rd_party/RepeatModeler/BuildDatabase -name $genome_name.rm -engine ncbi $genome >/dev/null && cd $cwd");
		sleep(10);
	}
	unless (-s "$genome.consensi.fa.classified"){
		my $local_cmd = "cd $cwd/species-specific && rm -rf RM_* genome.fasta && ln -s $genome genome.fasta && $RealBin/../3rd_party/RepeatModeler/RepeatModeler -engine ncbi -database $genome_name.rm -pa $repeatcpus 2>&1 > repeatmodelling.log && ln -s RM_*/consensi.fa.classified $genome.consensi.fa.classified && cd $cwd/";
  		my $thread = threads->create('just_run_my_commands_helper', $local_cmd, undef, undef);
	        $thread_helper->add_thread($thread);
		push(@threads_submitted,$thread);
 		sleep(10);
	}

   # wait for all to finish
	while (@threads_submitted && scalar(@threads_submitted)>0){
	   sleep(600);
   	   for (my $i=0;$i<(@threads_submitted);$i++){
		my $thread = $threads_submitted[$i] || next;
	        next if ($thread->is_running);
       		my $thread_id = $thread->tid;
		if ($thread->is_joinable){
			$thread->join();
		        warn "GOOD: thread $thread_id completed\n";
			splice(@threads_submitted,$i,1)
		}
	        if (my $error = $thread->error()) {
		        warn "ERROR: thread $thread_id exited with error $error\n";
		}
    	   }
	}
   # "species-specific"
	my $local_cmd = "cd $cwd/species-specific && $cmd -no_is -nolow -lib $genome.consensi.fa.classified genome.fasta 2>&1 > repeatmasking.log && cd ..";
	if (!-s "$genome.consensi.fa.classified"){
		warn "RepeatModeller failed!";
	}else{
		&process_cmd($local_cmd) unless -s "species-specific/genome.fasta.cat.gz";
	}


  # final file
  my $repeat_gff_file = $genome.".out.gff";
  &process_cmd("cat simple-only/genome.fasta.out.gff rna-specific/genome.fasta.out.gff species-specific/genome.fasta.out.gff general/genome.fasta.out.gff contamination-check/genome.fasta.out.gff | sort -nk1,5|uniq > $repeat_gff_file");

  die "Repeatmasking failed to produce $repeat_gff_file\n" unless -s $repeat_gff_file;
  my $repeat_hints_file = $repeat_gff_file;
  $repeat_hints_file =~s/.out.gff$/.repeatmasked.hints/;
  open (IN,$repeat_gff_file);
  open (OUT,">$repeat_hints_file");
  while (my $ln=<IN>){
	next if $ln=~/^#/;
	my @data = split("\t",$ln);
	$data[2] = 'nonexonpart';
	$data[8] = 'src=RM;pri=6';
	print OUT join("\t",@data)."\n";
  }
  close OUT;
  close IN;

  #produce $genome.masked
  &process_cmd("$bedtools_exec maskfasta -fi $genome -fo $genome.softmasked -bed $repeat_gff_file -soft");
  &process_cmd("$bedtools_exec maskfasta -fi $genome -fo $genome.hardmasked -bed $repeat_gff_file");
  symlink("$genome.hardmasked","$genome.masked");

  print "RepeatMasking completed and hints file produced ($repeat_gff_file)\n";
  return $repeat_gff_file;
}


sub check_for_iupac_violation(){
	my $fasta = shift;
	return if -f $fasta.'.checked';
	print "Checking FASTA $fasta for non ATCGN characters...\n";
	my $orig_sep = $/;
	$/ = ">";
	open (IN,$fasta);
	while (my $rec = <IN>){
		chomp($rec); next unless $rec;
		my @lines = split("\n",$rec);
		my $id = shift(@lines);
		my $seq = join('',@lines);
		$seq=~s/\s+//g;
		die "FASTA file $fasta has non-ATCGN bases\n" if $seq=~/[^ATCGN]/;
	}
	close IN;
	print "None found :-)\n";
	open(T,">$fasta.checked");
	print T "Checked\n";
	close T;
}
