#!/usr/bin/env perl

=pod

=head1 TODO
 
 new trimmomatic method about quality trimming (non-stringent vs stringent)

 add blue

 add kraken

 multiple files

 remove options not needed

 add quality sanitization (too high values for sanger)

 add length sanitization (too long reads 251 bp. last base is wrong)

=head1 NAME

 preprocess_illumina_gdna.pl

=head1 USAGE

 Preprocess Illumina data and create archives. This involves trimming based on quality, converting to FASTA, getting rid of Ns etc.
 Simply give one or more FASTQ files either .bz2, .gz or uncompressed

 Needs pbzip2
    
    -debug           => Verbose output. Prints out commands before they are run. Tracks memory for native deduplication

 Options
    -paired          => If 2 files have been provided, then treat them as a pair.
    -stop_qc         => Stop after QC of untrimmed file (e.g. in order to specify -trim_5 or -qtrim)
    -backup          => If bz2 files provided, then re-compress them using parallel-bzip2

 Screening:
    -adapters        => Illumina adapters FASTA (default provided)
    -noadapter       => Do not search for adapters
    -deduplicate :s  => Perform deduplication. There are two approaches. Native and Allpaths. For Native provide an integer to use as the length of the 5' seed (give '1' to use default of 16 for short reads and 32 for long reads). To use Allpaths provide a read name prefix this will overwrite the original readname. Allpaths used 16 bp as the 5' seed. Both approaches need lots of memory.

 These happen after any adapter trimming (in this order)
    -trim_5      :i     => Trim these many bases from the 5'. Happens before quality trimming but after adapter trimming (def 0)
    -max_keep    :i     => After any -trim5 then trim 3' end so that it is no longer than these many bases. Have seen erroneous 251th base in 250 bp sequencing (def automatic to closest whole 10 b.p decrement - 100, 150, 170 etc - if not user specified)
    -min_length  :i     => Discard sequences shorter than this (after quality trimming). Defaults to 32. Increase to 50-80 if you plan to use if it for alignments
    -qtrim       :i     => Trim 3' so that mean quality is that much in the phred scale (def. 5)
    -no_average_quality => Do not do any average quality trimming, just -trim_5, -max_keep and -min_length
    

 Quality is auto-calculated based on the minimum number found but
    -sanger          => Force sanger Quality scores for fastq (otherwise autodetect)
    -illumina        => Force Illumina 1.3-1.7 qual scores (autodetect)
    -casava18        => Force input as Fastq from Casava 1.8 (autodetect)
    -convert_fastq   => Convert to Sanger FASTQ flavour if Illumina 1.3 format is detected
    -dofasta         => Create FASTA file

Alexie tip: For RNA-Seq,            I check the FASTQC report of the processed data but do not trim the beginning low complexity regions (hexamer priming) as some tests with TrinityRNAseq did not show improvement (the opposite in fact).

=head1 DEPENDECIES

 * pbzip2/gzip

=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Digest::MD5 qw/md5/;
use BSD::Resource;
use FindBin qw/$RealBin/;
$ENV{PATH} .= ":$RealBin:$RealBin/3rd_party/FastQC/:$RealBin/3rd_party/allpaths/";

my ( $pbzip_exec, $fastqc_exec,$java_exec ) = &check_program( 'pbzip2', 'fastqc','java' );

my (
     $is_sanger,   $do_fasta,     $genome_size,  $use_meryl,
     $is_illumina, $delete_fastb, $is_cdna,      $no_preprocess,
     $is_casava,      @user_labels,  @user_bowties, $convert_fastq,
     $is_paired,   $trim_5,       $stop_qc,      $no_screen,
     $backup_bz2,  $debug,        $is_gdna,      $nohuman, 
     $noadapters, $max_keep_3, $no_qc, $mate_pair,$do_deduplicate, $max_length,
     $no_av_quality
);
my $cwd = `pwd`;
chomp($cwd);
my $kmer_ram   = int(0);
my $cpus       = 4;
my $qtrim      = 5;
my $min_length = 32;
my $slide_window  = 8 ; 
my $slide_quality  = 8 ; 

# edit these if you use it often with the same variables
my $trimmomatic_exec = $RealBin . "/3rd_party/Trimmomatic-0.36/trimmomatic-0.36.jar";
my $rDNA_db          = $RealBin . '/dbs/' . 'rDNA_nt_inv.fsa_nr';      #bowtie2
my $contam_db = $RealBin . '/dbs/' . 'ecoli_pseudomonas.fsa.masked.nr'; #bowtie2
my $human_db  = $RealBin . '/dbs/' . 'human_genome.fasta';              #bowtie2
my $phix_db   = $RealBin . '/dbs/' . 'phage_phiX174';                   #bowtie2
my $adapters_db = $RealBin . '/dbs/' . 'illumina_PE_2008_adapters.fsa'; #fasta

GetOptions(
            'debug'              => \$debug,
            'rDNA:s'             => \$rDNA_db,
            'contam:s'           => \$contam_db,
            'human:s'            => \$human_db,
            'nohuman'            => \$nohuman,
            'adaptors|adapters:s' => \$adapters_db,
            'noscreen|no_screen' => \$no_screen,
            'noqc|no_qc' => \$no_qc,
            'sanger'             => \$is_sanger,
            'dofasta'            => \$do_fasta,
            'genome_size:s'      => \$genome_size,
            'kmer_ram:i'         => \$kmer_ram,
            'meryl'              => \$use_meryl,
            'illumina'           => \$is_illumina,
            'delete_fastb'       => \$delete_fastb,
            'cdna'               => \$is_cdna,
            'gdna'               => \$is_gdna,
            'no_preprocess'      => \$no_preprocess,
            'casava18'           => \$is_casava,
            'user_bowtie:s{,}'   => \@user_bowties,
            'user_label:s{,}'    => \@user_labels,
            'convert_fastq'      => \$convert_fastq,
            'paired'             => \$is_paired,
            'trim_5:i'           => \$trim_5,
            'max_keep:i'         => \$max_keep_3,
            'stop_qc'            => \$stop_qc,
            'qtrim:i'            => \$qtrim,
            'backup'             => \$backup_bz2,
            'trimmomatic:s'      => \$trimmomatic_exec,
            'noadapters'         => \$noadapters,
            'min_length:i'       => \$min_length,
	    'slide_quality:i'    => \$slide_quality,
	    'slide_window:i'    => \$slide_window,
            'no_average_quality' => \$no_av_quality,
	    'mate_pair'=> \$mate_pair,
	    'deduplicate:s' => \$do_deduplicate
) || pod2usage();
my @files = @ARGV;

pod2usage "Trimmomatic is now required.\n" if !$trimmomatic_exec;

pod2usage
  "Cannot use paired option without specifying two (and only two) files\n"
  if ( $is_paired && scalar(@files) != 2 );
if ( !$is_paired && scalar(@files) == 2 ) {
 warn
"You provided two files. If these are paired/mated files, then you should also specify -paired\n";
 sleep(3);
}

for ( my $i = 0 ; $i < scalar(@user_bowties) ; $i++ ) {
 my $user_bowtie = $user_bowties[$i];
 die "It doesn't seem that $user_bowtie is formatted for Bowtie2\n" unless -s $user_bowtie.'.1.bt2';
 my $user_label  = $user_labels[$i];
 pod2usage "Must specify a label if also using a custom bowtie file\n"
   if $user_bowtie && !$user_label;
}
pod2usage "No files given\n" unless @files;

undef($adapters_db) if $noadapters;

#setrlimit( RLIMIT_VMEM, $kmer_ram * 1000 * 1000 , $kmer_ram * 1024 * 1024 )  if $kmer_ram;

my ( %files_to_delete_master );

for ( my $i = 0 ; $i < @files ; $i++ ) {

# warn "Requested -kmer_ram is less than file size. May cause problems.\n" if ($kmer_ram && ($kmer_ram * 1024 * 1024 * 1024) < -s $file;
 my ( $file_is_sanger, $file_is_illumina, $file_is_casava, $file_format );
 my $file = $files[$i];
 warn "Cannot find $file" if !-s $file;
 sleep(1) if !-s $file;
 next if $file =~ /trim.filtered/;
 my $out = $file;
 $out =~ s/.bz2$//;
 $out =~ s/.gz$//;
 $out .= '.trimmomatic' if ( $is_paired && $trimmomatic_exec && $adapters_db );
 $out .= '.trim.filtered.fastq';
 print "\n#############################\nProcessing filename $file\n";

 if (    ( -f $out && ( -s $out ) > 104857600 )
      || ( -f $out . '.bz2' && ( -s $out . '.bz2' ) > 10485760 ) )
 {
  warn "$out seems to exist already. Using already uncompressed file\n";
  &process_cmd("$pbzip_exec -dkp8 $out.bz2");
  $file =~ s/.bz2$// if ( $file =~ /.bz2$/ );
  $file =~ s/.gz$//  if ( $file =~ /.gz$/ );
  $files[$i] = $file;
 }
 else {
  if ( $file =~ /.bz2$/ && $backup_bz2 ) {
   print "Uncompressing $file\n";
   &process_cmd("$pbzip_exec -dp4 $file")
     ; # baylor files are compressed with single threaded bzip2 so make a $pbzip_exec backup
   $file =~ s/.bz2$//;
  }
  elsif ( $file =~ /.bz2$/ ) {
   print "Uncompressing $file\n";
   &process_cmd("$pbzip_exec -kdp4 $file");
   $file =~ s/.bz2$//;
  }
  elsif ( $file =~ /.gz$/ ) {
   print "Uncompressing $file\n";
   $file =~ s/.gz$//;
   &process_cmd("gunzip  $file.gz ");
  }
  $files_to_delete_master{$file} = 1;

  $files[$i] = $file;
   my $check = &check_fastq_format($file);
   $file_is_sanger   = 1 if $check eq 'sanger';
   $file_is_casava   = 1 if $check eq 'casava';
   $file_is_illumina = 1 if $check eq 'illumina';
  if ($file_is_casava) {
   $file_is_sanger = 1;
   undef($file_is_illumina);
   print "Converting from CASAVA 1.8 to Sanger header format\n"  if $convert_fastq;
   &process_cmd(
             'sed --in-place \'s/^@\([^ ]*\) \([0-9]\).*/@\1\/\2/\' ' . $file ) if $convert_fastq;
  }
  $files_to_delete_master{$file} = 1;
  my $fastqc_basename = $file;$fastqc_basename=~s/\.[^\.\-\_]+$//;$fastqc_basename.='_fastqc'; # probably
  unless (-s $fastqc_basename . ".zip" || !$fastqc_exec || $no_qc){
    system("$fastqc_exec --noextract --nogroup -q $file");
  }
 }
}

if ($stop_qc) {
 print "User asked to stop after QC\n";
 exit(0);
}

###############################
# TRIMMING
###############################
if (!$no_preprocess){
if ( $is_paired && $trimmomatic_exec) {
 my $file1 = $files[0];
 my $file2 = $files[1];
 print "Pre-processing $file1 and $file2\n";
 my $check1 = &check_fastq_format($file1);
 my $check2 = &check_fastq_format($file2);
 my $cmd = "java -jar $trimmomatic_exec PE -threads $cpus -phred33 $file1 $file2 $file1.trimmomatic $file1.trimmomatic.unpaired $file2.trimmomatic $file2.trimmomatic.unpaired MINLEN:32 ";
 $cmd .= " ILLUMINACLIP:$adapters_db:2:40:15 " if $adapters_db ;
 $cmd .= " HEADCROP:$trim_5 " if $trim_5;
 $cmd .= " CROP:$max_keep_3 " if $max_keep_3 && $max_keep_3 >0;
 $cmd .= " LEADING:4 TRAILING:$qtrim SLIDINGWINDOW:$slide_window:$slide_quality " unless $no_av_quality;
 if ( $check1 eq $check2 && ( $check1 eq 'illumina' ) ) {
  $cmd =~ s/phred33/phred64/;
  $cmd .= " TOPHRED33 ";
 }
 $cmd .= " MINLEN:$min_length 2>&1";
 &process_cmd($cmd) unless -s "$file1.trimmomatic" && -s "$file2.trimmomatic";
 die "Something bad happened... one of the files are empty/missing...\n" unless -s "$file1.trimmomatic" && -s "$file2.trimmomatic";
 $files_to_delete_master{$file1}                 = 1;
 $files_to_delete_master{ $file1 . '.trimmomatic.unpaired' } = 1;
 $files_to_delete_master{$file2}                 = 1;
 $files_to_delete_master{ $file2 . '.trimmomatic.unpaired' } = 1;

 $files[0] .= '.trimmomatic';
 $files[1] .= '.trimmomatic';
}else {
 for ( my $i = 0 ; $i < @files ; $i++ ) {
  my $file = $files[$i];
  print "Pre-processing $file\n";
  my $check = &check_fastq_format($file);
  my $cmd = "java -jar $trimmomatic_exec SE -threads $cpus -phred33 $file $file.trimmomatic MINLEN:32 ";
  $cmd .= " ILLUMINACLIP:$adapters_db:2:40:15 " if  $adapters_db;
  $cmd .= " HEADCROP:$trim_5 " if $trim_5;
  $cmd .= " CROP:$max_keep_3 " if $max_keep_3 && $max_keep_3 >0;
  $cmd .= " LEADING:4 TRAILING:$qtrim SLIDINGWINDOW:$slide_window:$slide_quality " unless $no_av_quality;

  if ( $check && $check eq 'illumina' ) {
   $cmd =~ s/phred33/phred64/;
   $cmd .= " TOPHRED33 ";
  }
  $cmd .= " MINLEN:$min_length 2>&1";
  &process_cmd($cmd) unless -s "$file.trimmomatic";
  die "Something bad happened... one of the files are empty/missing...\n" unless -s "$file.trimmomatic";
  $files_to_delete_master{$file} = 1;
  $files[$i] .= '.trimmomatic';
 }
}
}

print "Stage 1 completed (".$files[0].")\n";
##########################

&remove_dodgy_reads($files[0],$files[1]) if $is_paired && $do_deduplicate;

##########################
foreach my $file (@files){
  my $fastqc_basename = $file;$fastqc_basename.='_fastqc'; # probably
  system("$fastqc_exec --noextract --nogroup -q $file") unless -s $fastqc_basename . ".zip" || !$fastqc_exec;
  $files_to_delete_master{$fastqc_basename.".html"} = 1;
}

##########################

for ( my $i = 0 ; $i < @files ; $i++ ) {
 my $file = $files[$i];
 print "Post-processing $file\n";
 &fastq_to_fasta("$file")  if $do_fasta;
 &do_genome_kmers("$file") if $kmer_ram > 0;
 $files_to_delete_master{$file} = 1;
}

print "Completed. Compressing/cleaning up...\n";

&process_cmd(   "$pbzip_exec -fvp$cpus "
              . join( " ", ( keys %files_to_delete_master ) )
              . " 2>/dev/null" )
  if %files_to_delete_master && scalar( keys %files_to_delete_master ) > 0;


##################################################################################################


sub check_fastq_format() {
 my $file = shift;
 
 die "No file or does not exist\n" unless $file && -s $file;
 return 'sanger'  if $is_sanger;
 return 'casava' if $is_casava;
 return 'illumina' if $is_illumina;
 my $max_seqs = 100000;
 my ($min_number,$max_number); 
 my ( @line, $l, $number, $counter );
 my $id; # store for future

 open( FQ, $file );
 for (my $counter=0;$counter<$max_seqs;$counter++){
        $id = <FQ> || last;
        my $seq = <FQ>;
        my $qid    = <FQ>;
        my $qual   = <FQ>;
  die "File is not a fastq file:\n$id\n" unless ( $id =~ /^@/ );
  chomp($qual);
  $max_length = length($qual) if !$max_length || $max_length < length($qual);
  @line = split( '', $qual );    # divide in chars
  for ( my $i = 0 ; $i < length($qual) ; $i++ ) {    # for each char
   $number = ord( $line[$i] );    # get the number represented by the ascii char
   $min_number = $number if (!$min_number || $number < $min_number);
   $max_number = $number if (!$max_number || $number > $max_number);
   # check if it is sanger or illumina/solexa, based on the ASCII image at http://en.wikipedia.org/wiki/FASTQ_format#Encoding
  }
 }
 close FQ;

 # use $max_length to determine if last base should be cut.
 if (!$max_keep_3){
	$max_keep_3 = $max_length;
        if ($max_keep_3 > 100){
		while ($max_keep_3 % 10 != 0 ){
			$max_keep_3--;
        	}
	}
 }


 if ( $min_number >= 75 ) {    # if solexa/illumina
    print "This file is solexa/illumina format ($min_number,$max_number) with max length $max_length\n";
    close FQ;
    return 'illumina';
  }
 elsif ($min_number < 50 ) {      # if sanger
  if ($id =~ /(\S+)\s*(\S*)/){
          my $description = $2;
          if ( $description && $description =~ /(\d)\:[A-Z]\:/ ) {
           print "This file is in Sanger quality but CASAVA 1.8 header format ($min_number,$max_number) with max length $max_length\n";
           return 'casava';
   }
  }
  print "This file is sanger format ($min_number,$max_number) with max length $max_length\n";    # print result to terminal and die
  return 'sanger';
 }
 die "Cannot determine fastq format ($min_number,$max_number) with max length $max_length\n";
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
 print "COMPLETED: $cmd ($ret)\n" if $debug;
 return;
}

sub illumina2sanger() {
 my $file = shift || die;
 my $outfile = $file . '.sanger';
 print "Using existing file $file.sanger\n" if -s $outfile;
 return $outfile if -s $outfile;
 open( IN,  $file );
 open( OUT, ">$outfile" );
 my $count = int(0);
 while ( my $id = <IN> ) {
  my $seq = <IN>;
  my $id2 = <IN>;
  my $qlt = <IN>;
  $count++;
  if ( !$seq ) {
   next;
  }
  if ( length($seq) ne length($qlt) ) {
   die "$file: Length of sequence and quality for id $id differs ("
     . length($seq) . " vs "
     . length($qlt)
     . "! sequence $count, line "
     . ( $count * 4 ) . "\n";
  }
  chomp($qlt);
  $qlt =~ tr/\x40-\xff\x00-\x3f/\x21-\xe0\x21/;
  print OUT $id . $seq . "+\n" . $qlt . "\n";
 }
 return $outfile;
}

sub check_program() {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  pod2usage "Error, path to a required program ($prog) cannot be found\n\n"
    unless $path =~ /^\//;
  chomp($path);
  $path = readlink($path) if -l $path;
  push( @paths, $path );
 }
 return @paths;
}


sub fastq_to_fasta() {
 my $infile = shift || die("No input file");
 my $output = $infile . ".fsa";
 return if ( -s $output || $output . '.bz2' );
 print "Preparing FASTA\n";
 open( IN,  $infile );
 open( OUT, ">$output" );
 my $count = 1;
 while ( my $ln = <IN> ) {
  next if $ln =~ /^\s*$/;
  if ( $ln =~ /^@/ ) {
   my $seq = <IN>;
   $ln =~ s/^@/>/;
   print OUT $ln . $seq;
   my $discard = <IN> . <IN>;
  }
  else {
   die "Invalid line $count found:\n$ln\n";
  }
  $count += 4;
 }
 close IN;
 close OUT;
 system("$pbzip_exec -p4 $output");
}

sub fastq_keep_trim3(){
    my ($file,$max_size)=@_;
    open (IN,$file) || die $!;
    open (OUT,">$file.x")||die $!;
    while (my $id=<IN>){
        my $seq = <IN>;
        my $qid = <IN>;
        my $qlt = <IN>;
        
        if (length($seq) > ($max_size+1)){
            chomp($seq);chomp($qlt);
            $seq = substr($seq,0,$max_size)."\n";
            $qlt = substr($qlt,0,$max_size)."\n";
            
        }
        print OUT $id.$seq.$qid.$qlt;
    }
    

    close IN;
    close OUT;
    rename("$file.x",$file);
}

sub get_path(){
    my $file = shift ||die;
    my $path = &check_program($file);
    die $! unless $path;
    $path = `which $file`;chomp($path);
    die unless $path;
    return $path;
}

sub remove_dodgy_reads(){
  if ($do_deduplicate=~/^\d+$/){
  remove_dodgy_reads_native(@_);
  }else{
      &remove_dodgy_reads_allpaths(@_);
  }
}


sub total_quality(){
  # get total quality
  my $total_q = int(0);
  foreach my $string (@_){
    chomp($string);
    my @array = split('',$string);
    for (my $i=0;$i<scalar(@array);$i++){
       $total_q += ord($array[$i]);
    }
  }
  return $total_q;
}

sub remove_dodgy_reads_native(){
    print "Removing duplicate PCR fragments...\n";
    # using the allpaths technique
    # assume paired; $max_length has the length
    my ($file1,$file2)=@_;
    my %hash;
    die "Sequences are too short for reliable deduplication\n" if $max_length < 16;
    my $size_search = $max_length < 60 ? 16 : 32;
    if ($do_deduplicate=~/^\d+$/ && $do_deduplicate != 1){
        $size_search = $do_deduplicate;
    }
    die "Search hash it too short ($size_search)\n" if $size_search < 8;
    print "Hashing files using K=$size_search\n";
    my $total_seqs=int(0);
    my $ignored_seqs = int(0);
    open (FILE1,$file1);
    open (FILE2,$file2);
    $|=1;
    while (my $sid1=<FILE1>){
        $total_seqs++;
        if ($total_seqs =~/00000$/){
            print "Processed $total_seqs   ";
            &get_memory_usage() if $debug;
            print "\r";
        }
        my $seq1 = <FILE1>;
        my $qid1 = <FILE1>;
        my $qlt1 = <FILE1>;
        
        my $sid2 = <FILE2>;
        my $seq2 = <FILE2>;
        my $qid2 = <FILE2>;
        my $qlt2 = <FILE2>;
        die "FATAL: Number of lines for second file ($file2) is not the same as first file ($file1)\n" unless $qlt2;
        
        next if length($seq1) < $size_search || length($seq2) < $size_search;
        my $id; # common ID
        if ($sid1=~/^(\S+)/){
            $id = $1;
            $id =~s/\/[12]$//;
        }
        
        my $md5_1 = md5(substr($seq1,$size_search));
        my $md5_2 = md5(substr($seq2,$size_search));
        if (($md5_1 cmp $md5_2) == 0){
            $ignored_seqs++;
#            warn "WARNING: possibly identical sequences ($id). Skipping.\n $seq1 $seq2\n";
            next;
        }
        elsif (($md5_1 cmp $md5_2) == 1){
            my $t = $md5_1;
            $md5_1 = $md5_2;
            $md5_2 = $t;
        }
        my $total_q = &total_quality($qlt1,$qlt2);


        # this could be more efficient if we used lists instead of hashes.
        
        # this could go into an in-memory/file sqlite
        if ($hash{$md5_1}{$md5_2}){
            if ($total_q > $hash{$md5_1}{$md5_2}{'q'}){
                $hash{$md5_1}{$md5_2}{'i'} = $id;
                $hash{$md5_1}{$md5_2}{'q'} = $total_q;
            }
        }else{
          $hash{$md5_1}{$md5_2}{'i'} = $id;
          $hash{$md5_1}{$md5_2}{'q'} = $total_q;
        }
    }
    close FILE1;
    close FILE2;
    $|=0;
    print "\nFound $ignored_seqs pairs with identical seed for forward and reverse sequences.\nDeduplicating files...\n";
    my $counter=int(0);
    
    open (FILE1,$file1);
    open (FILE2,$file2);
    open (OUT1,">$file1.dedup");
    open (OUT2,">$file2.dedup");
    while (my $sid1=<FILE1>){
        my $seq1 = <FILE1>;
        my $qid1 = <FILE1>;
        my $qlt1 = <FILE1>;
        
        my $sid2 = <FILE2>;
        my $seq2 = <FILE2>;
        my $qid2 = <FILE2>;
        my $qlt2 = <FILE2>;
        
        
        next if (length($seq1) < $size_search || length($seq2) < $size_search);
        my $id;
        if ($sid1=~/^(\S+)/){
            $id = $1;
            $id =~s/\/[12]$//;
        }
        
        my $md5_1 = md5(substr($seq1,$size_search));
        my $md5_2 = md5(substr($seq2,$size_search));
        if (($md5_2 cmp $md5_1) == 1){
            my $t = $md5_1;
            $md5_1 = $md5_2;
            $md5_2 = $t;
        }
          
        next unless ($hash{$md5_1}{$md5_2} && $hash{$md5_1}{$md5_2}{'i'} eq $id) || ($hash{$md5_2}{$md5_1} && $hash{$md5_2}{$md5_1}{'i'} eq $id) ;
        print OUT1 $sid1.$seq1.$qid1.$qlt1;
        print OUT2 $sid2.$seq2.$qid2.$qlt2;
        $counter++;
    }
    close FILE1;
    close FILE2;
    close OUT1;
    close OUT2;
    print "Deduplication done. Kept $counter sequences from a total of $total_seqs sequences\n";
    rename("$file1.dedup",$file1);
    rename("$file2.dedup",$file2);
}


sub remove_dodgy_reads_allpaths(){
    print "Removing duplicate PCR fragments...\n";
    my $FastbQualbToFastq_exec = &get_path('FastbQualbToFastq');
    my $FastqToFastbQualb_exec = &get_path('FastqToFastbQualb');
    my $RemoveDodgyReads_exec = &get_path('RemoveDodgyReads');
    my $MergePairedFastbs_exec  = &get_path('MergePairedFastbs');
    # sadly using allpaths for deduplication... one day write it's code to use fastq.
    my ($file1,$file2)=@_;
    die unless $file1 && $file2 && -s $file1 && -s $file2;    
    my $cmd = "$FastqToFastbQualb_exec WRITE_QUALB=True ";

    &process_cmd($cmd." FASTQ=$file1 OUT_HEAD=$file1");
    &process_cmd($cmd." FASTQ=$file2 OUT_HEAD=$file2");
    &process_cmd($MergePairedFastbs_exec. " HEAD1_IN=$file1 HEAD2_IN=$file2 HEAD_OUT=tmp WRITE_PAIRS=True");
    die unless -s "tmp.pairs";
    unlink("$file1.fastb","$file1.qualb");
    unlink("$file2.fastb","$file2.qualb");
    
    &process_cmd($RemoveDodgyReads_exec." REMOVE_POLY_A=False NUM_THREADS=$cpus IN_HEAD=tmp OUT_HEAD=clean");
    unlink("tmp.fastb","tmp.qualb","tmp.pairs");
    die unless -s "clean.readtrack";
    
    
    &process_cmd("$FastbQualbToFastq_exec HEAD_IN=clean HEAD_OUT=clean PAIRED=True PHRED_OFFSET=33 PICARD_NAMING_SCHEME=True NAMING_PREFIX=$do_deduplicate");
    unlink("clean.fastb","clean.qualb","clean.pairs","clean.readtrack");
    die "Deduplication failed\n" unless -s "clean.A.fastq" && -s "clean.B.fastq";
    rename("clean.A.fastq",$file1);
    rename("clean.B.fastq",$file2);
}

sub get_memory_usage(){
    my $rusage = getrusage();
    print "\tMemory usage: ".$rusage->maxrss."k\n";
}

