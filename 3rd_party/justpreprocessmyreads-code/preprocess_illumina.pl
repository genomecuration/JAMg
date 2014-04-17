#!/usr/bin/env perl

=pod


=head1 NAME

 preprocess_illumina_gdna.pl

=head1 USAGE

 Preprocess Illumina data and create archives. This involves trimming based on quality, converting to FASTA, getting rid of Ns etc.
 Simply give one or more FASTQ files either .bz2, .gz or uncompressed

 Requires ALLPATHS-LG or meryl for kmer-hist

 Needs pbzip2
    
    debug           => Verbose output. Prints out commands before they are run

 Options
    paired          => If 2 files have been provided, then treat them as a pair.
    cdna            => Input is cDNA
    gdna            => Input is gDNA
    no_preprocess   => Do not do any trimming    
    stop_qc         => Stop after QC of untrimmed file (e.g. in order to specify -trim_5 or -qtrim)
    backup          => If bz2 files provided, then re-compress them using parallel-bzip2 (e.g. if source is Baylor)
    no_screen       => Do not screen for contaminants using bowtie2

 Screening options:
    rDNA        :s  => rDNA bowtie2 database
    contam      :s  => contaminants bowtie2 database
    human       :s  => human genome bowtie2 database (unless you're sequencing humans!)
    nohuman         => Don't search the human genome
    adapters        => Illumina adapters FASTA if also -trimmomatic (default provided)
    noadaptor       => Do not search for adaptors
    trimmomatic     => Trimmomatic jar file
    user_bowtie :s  => 1+ bowtie2 database array for optional custom screening. Must match -user_label
    user_label  :s  => 1+ this label will be used for custom screening. Must match -user_bowtie
    
 These happen after any adaptor trimming (in this order)
    trim_5      :i  => Trim these many bases from the 5' (def 0)
    max_keep    :i  => Trim 3' end so that it is no longer than these many bases (def 0)
    qtrim       :i  => Trim 3' so that mean quality is that much in the phred scale (def. 5)
    max_length  :i  => Discard sequences shorter than this (after quality trimming). Defaults to 32. Increase to 50-80 if you plan to use if it for alignments

 These are not really used:
    genome_size  :s => Give genome size to estimate coverage (can use mb or gb as suffix for megabases or gigabases)
    kmer_ram:i      => To create Fastb and produce a 25-mer histogram, then tell me how many GB of RAM to use (0 switches it off)
    threads:i       => Threads for building kmer table using allpaths (def 4). Meryl will only use 1 thread. 
                      NB: The program also uses 2-4 threads for parallel stats creation & compressing/backup of files but 
                      these ought to have finished by the time kmer building starts (but this is not checked)
    meryl           => force the use of meryl (keeps 25mer database)
    sanger          => Force sanger Quality scores for fastq (otherwise autodetect)
    illumina        => Force Illumina 1.3-1.7 qual scores (autodetect)
    casava18        => Force input as Fastq from Casava 1.8 (autodetect)
    noconvert_fastq => Don't convert to Sanger FASTQ flavour if Illumina 1.3 format is detected
    dofasta         => Create FASTA file

Alexie tip: For RNA-Seq, I check the FASTQC report of the processed data but do not trim the beginning low complexity regions (hexamer priming) as some tests with TrinityRNAseq did not show improvement (the opposite in fact).

I have seen Trimmomatic to inadvertedly shuffle reads when threads are used. So it runs with a single thread here.

=head1 DEPENDECIES

 * pbzip2/gzip
 * fastqc
 * build_illumina_mates.pl (author: Alexie)
 * samtools
 * bowtie2
 * AllpathsLG or meryl (optional)
 * R
 * sed, ps2pdf

=cut

use strict;
use warnings;
use Data::Dumper;
use threads;
use Getopt::Long;
use Pod::Usage;
use BSD::Resource;
use FindBin qw/$RealBin/;
$ENV{PATH} .= ":$RealBin:$RealBin/3rd_party/FastQC/";

my ( $pbzip_exec, $fastqc_exec ) = &check_program( 'pbzip2', 'fastqc' );

my (
     $is_sanger,   $do_fasta,     $genome_size,  $use_meryl,
     $is_illumina, $delete_fastb, $is_cdna,      $no_preprocess,
     $casava,      @user_labels,  @user_bowties, $noconvert_fastq,
     $is_paired,   $trim_5,       $stop_qc,      $no_screen,
     $backup_bz2,  $debug,        $is_gdna,      $nohuman,
     $noadaptors, $trim_3
);
my $cwd = `pwd`;
chomp($cwd);
my $kmer_ram   = int(0);
my $cpus       = 4;
my $qtrim      = 5;
my $max_length = 32;

# edit these if you use it often with the same variables
my $trimmomatic_exec = $RealBin . "/3rd_party/trimmomatic-0.32.jar";
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
            'adapters:s'         => \$adapters_db,
            'noscreen|no_screen' => \$no_screen,
            'sanger'             => \$is_sanger,
            'dofasta'            => \$do_fasta,
            'genome_size:s'      => \$genome_size,
            'kmer_ram:i'         => \$kmer_ram,
            'cpu|threads:i'      => \$cpus,
            'meryl'              => \$use_meryl,
            'illumina'           => \$is_illumina,
            'delete_fastb'       => \$delete_fastb,
            'cdna'               => \$is_cdna,
            'gdna'               => \$is_gdna,
            'no_preprocess'      => \$no_preprocess,
            'casava18'           => \$casava,
            'user_bowtie:s{,}'   => \@user_bowties,
            'user_label:s{,}'    => \@user_labels,
            'noconvert_fastq'    => \$noconvert_fastq,
            'paired'             => \$is_paired,
            'trim_5:i'           => \$trim_5,
            'max_keep:i'           => \$trim_3,
            'stop_qc'            => \$stop_qc,
            'qtrim:i'            => \$qtrim,
            'backup'             => \$backup_bz2,
            'trimmomatic:s'      => \$trimmomatic_exec,
            'noadaptors'         => \$noadaptors,
            'max_length:i'       => \$max_length
);
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

undef($adapters_db) if $noadaptors;

#setrlimit( RLIMIT_VMEM, $kmer_ram * 1000 * 1000 , $kmer_ram * 1024 * 1024 )  if $kmer_ram;

my ( @threads, %files_to_delete_master );

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
  if ( !$casava && !$is_sanger && !$is_illumina ) {
   my $check = &check_fastq_format($file);
   $file_is_sanger   = 1 if $check eq 'sanger';
   $file_is_casava   = 1 if $check eq 'casava';
   $file_is_illumina = 1 if $check eq 'illumina';
  }
  else {
   $file_is_sanger   = 1 if $is_sanger;
   $file_is_casava   = 1 if $casava;
   $file_is_illumina = 1 if $is_illumina;
  }
  if ($file_is_casava) {
   $file_is_sanger = 1;
   undef($file_is_illumina);
   print "Converting from CASAVA 1.8 to Sanger header format\n";
   &process_cmd(
             'sed --in-place \'s/^@\([^ ]*\) \([0-9]\).*/@\1\/\2/\' ' . $file );
  }
  $files_to_delete_master{$file} = 1;
  push( @threads, &create_fork("$fastqc_exec --noextract --nogroup -q $file") )
    unless -s $file . "_fastqc.zip" || !$fastqc_exec;
 }
}

if ($stop_qc) {
 foreach my $thread (@threads) {
  $thread->join();
 }
 print "User asked to stop after QC\n";
 exit(0);
}

###############################
if ( $is_paired && $trimmomatic_exec ) {
 my $file1 = $files[0];
 my $file2 = $files[1];
 print "Pre-processing $file1 and $file2\n";
 my $check1 = &check_fastq_format($file1);
 my $check2 = &check_fastq_format($file2);
 my $cmd = "java -jar $trimmomatic_exec PE -threads 1 -phred33 $file1 $file2 $file1.trimmomatic $file1.unpaired $file2.trimmomatic $file2.unpaired MINLEN:32 ";
 $cmd .= " ILLUMINACLIP:$adapters_db:2:40:15 " if $adapters_db ;
 $cmd .= " HEADCROP:$trim_5 " if $trim_5;
 $cmd .= " CROP:$trim_3 " if $trim_3;
 $cmd .= " LEADING:4 TRAILING:$qtrim SLIDINGWINDOW:8:10 ";
 if ( $check1 eq $check2 && ( $check1 eq 'illumina' ) ) {
  $cmd =~ s/phred33/phred64/;
  $cmd .= " TOPHRED33 ";
 }
 $cmd .= " MINLEN:$max_length 2>&1";
 &process_cmd($cmd) unless -s "$file1.trimmomatic" && -s "$file2.trimmomatic";
 die "Something bad happened... one of the files are empty/missing...\n" unless -s "$file1.trimmomatic" && -s "$file2.trimmomatic";
 $files_to_delete_master{$file1}                 = 1;
 $files_to_delete_master{ $file1 . '.unpaired' } = 1;
 $files_to_delete_master{$file2}                 = 1;
 $files_to_delete_master{ $file2 . '.unpaired' } = 1;

 $files[0] .= '.trimmomatic';
 $files[1] .= '.trimmomatic';
 push(
       @threads,
       &create_fork(
                    "$fastqc_exec --noextract --nogroup -q $file1.trimmomatic ")
 ) unless -f "$file1.trimmomatic_fastqc.zip" || !$fastqc_exec;
 push(
       @threads,
       &create_fork(
                    "$fastqc_exec --noextract --nogroup -q $file2.trimmomatic ")
 ) unless -f "$file2.trimmomatic_fastqc.zip" || !$fastqc_exec;
}
else {
 for ( my $i = 0 ; $i < @files ; $i++ ) {
  my $file = $files[$i];
  print "Pre-processing $file\n";
  my $check = &check_fastq_format($file);
  my $cmd = "java -jar $trimmomatic_exec SE -threads 1 -phred33 $file $file.trimmomatic MINLEN:32 ";
  $cmd .= " ILLUMINACLIP:$adapters_db:2:40:15 " if  $adapters_db;
  $cmd .= " HEADCROP:$trim_5 " if $trim_5;
  $cmd .= " CROP:$trim_3 " if $trim_3;
  $cmd .= " LEADING:4 TRAILING:$qtrim SLIDINGWINDOW:8:10 ";

  if ( $check && $check eq 'illumina' ) {
   $cmd =~ s/phred33/phred64/;
   $cmd .= " TOPHRED33 ";
  }
  $cmd .= " MINLEN:$max_length 2>&1";
  &process_cmd($cmd) unless -s "$file.trimmomatic";
  die "Something bad happened... one of the files are empty/missing...\n" unless -s "$file.trimmomatic";
  $files_to_delete_master{$file} = 1;
  $files_to_delete_master{ $file . '.unpaired' } = 1;
  $files[$i] .= '.trimmomatic';
  push(
        @threads,
        &create_fork(
                     "$fastqc_exec --noextract --nogroup -q $file.trimmomatic ")
  ) unless -f "$file.trimmomatic_fastqc.zip" || !$fastqc_exec;
 }
}

print "Stage 1 completed\n";
##########################

for ( my $i = 0 ; $i < @files ; $i++ ) {
 my $file = $files[$i];
 print "Post-processing $file\n";
 &fastq_to_fasta("$file")  if $do_fasta;
 &do_genome_kmers("$file") if $kmer_ram > 0;
 &screen_bowtie2( $files[$i] ) unless $no_screen || !-s $files[$i];
 $files_to_delete_master{$file} = 1;
}

print "Completed. Compressing/cleaning up...\n";
foreach my $thread (@threads) {
 $thread->join();
}

&process_cmd(   "$pbzip_exec -fvp$cpus "
              . join( " ", ( keys %files_to_delete_master ) )
              . " 2>/dev/null" )
  if %files_to_delete_master && scalar( keys %files_to_delete_master ) > 0;

##################################################################################################
sub screen_bowtie2() {
 my $file = shift;
 my ( $bowtie2_exec, $samtools_exec ) = &check_program( 'bowtie2', 'samtools' );
 print "Screening for contaminants\n";
 if ( $is_cdna && !-s $file . "_vs_rDNA.bam" && $rDNA_db ) {
  print "\tInvertrebrate rRNA\n";

#  &process_cmd("kanga -T $cpus -i $file -I $rDNA_db -r2 -R300 -M5 -F "  . $file . "_vs_rDNA.log -o " . $file . "_vs_rDNA.sam" );
  &process_cmd(
"$bowtie2_exec --threads $cpus --fast --no-unal -x $rDNA_db -U $file 2> $file"
     . "_vs_rDNA.log |$samtools_exec view -S -F4 -b -o $file"
     . "_vs_rDNA.bam - > /dev/null" )
    unless -s $file . "_vs_rDNA.log";
  system(   "$samtools_exec view $file"
          . "_vs_rDNA.bam |cut -f 2|sort -n |uniq -c >> $file"
          . "_vs_rDNA.log " );
  sleep(3);
 }
 print "\tEcoli/Pseudomonas\n";
 if ( !-s $file . "_vs_contam.bam" && $contam_db ) {

#  &process_cmd( "kanga -T $cpus -i $file -I $contam_db -r2 -R300 -M5 -F "    . $file   . "_vs_contam.log -o "   . $file   . "_vs_contam.sam" );
  &process_cmd(
"$bowtie2_exec --threads $cpus --fast --no-unal -x $contam_db -U $file 2> $file"
     . "_vs_contam.log |$samtools_exec view -S -F4 -b -o $file"
     . "_vs_contam.bam - > /dev/null" )
    unless -s $file . "_vs_contam.log";
  system(   "$samtools_exec view $file"
          . "_vs_contam.bam |cut -f 2|sort -n |uniq -c >> $file"
          . "_vs_contam.log " );
  sleep(3);
 }

 print "\tPhiX 174 / Illumina control\n";
 if ( !-s $file . "_vs_phix.bam" && $phix_db ) {

#  &process_cmd(  "kanga -T $cpus -i $file -I $phix_db -r2 -R300 -M5 -F " . $file    . "_vs_phix.log -o "    . $file    . "_vs_phix.sam" );
  &process_cmd(
"$bowtie2_exec --threads $cpus --fast --no-unal -x $phix_db -U $file 2> $file"
     . "_vs_phix.log |$samtools_exec view -S -F4 -b -o $file"
     . "_vs_phix.bam - > /dev/null" )
    unless -s $file . "_vs_phix.log";
  system(   "$samtools_exec view $file"
          . "_vs_phix.bam |cut -f 2|sort -n |uniq -c >> $file"
          . "_vs_phix.log " );
  sleep(3);
 }

 if ( !-s $file . "_vs_human.bam" && $human_db && !$nohuman ) {
  print "\tHuman genome\n";

#  &process_cmd(  "kanga -T $cpus -i $file -I $human_db -r0 -R300 -M5 -F " . $file    . "_vs_human.log -o "    . $file    . "_vs_human.sam" );
  &process_cmd(
"$bowtie2_exec --threads $cpus --fast --no-unal -x $human_db -U $file 2> $file"
     . "_vs_human.log |$samtools_exec view -S -F4 -b -o  $file"
     . "_vs_human.bam - > /dev/null" )
    unless -s $file . "_vs_human.log";
  system(   "$samtools_exec view $file"
          . "_vs_human.bam |cut -f 2|sort -n |uniq -c >> $file"
          . "_vs_human.log " );
  sleep(3);
 }

 for ( my $i = 0 ; $i < scalar(@user_bowties) ; $i++ ) {
  my $user_bowtie = $user_bowties[$i];
  my $user_label  = $user_labels[$i];
  if ( $user_bowtie && -s $user_bowtie.'.1.bt2') {
   print "\tUser specified bowtie file $user_bowtie (labelled $user_label)\n";
   if ( !-s $file . "_vs_$user_label.bam" ) {

#   &process_cmd(   "kanga -T $cpus -i $file -I $user_bowtie -r2 -R300 -M5 -F "  . $file . "_vs_$user_label.log -o " . $file . "_vs_$user_label.sam" );
    &process_cmd(
"$bowtie2_exec --threads $cpus --fast --no-unal -x $user_bowtie -U $file 2> $file"
       . "_vs_$user_label.log | $samtools_exec view -S -F4 -b -o  $file"
       . "_vs_$user_label.bam - > /dev/null" )
      unless -s $file . "_vs_$user_label.log";
    system(   "$samtools_exec view $file"
            . "_vs_$user_label.bam |cut -f 2|sort -n |uniq -c >> $file"
            . "_vs_$user_label.log " );
    sleep(3);
   }
  }
 }
}

sub check_fastq_format() {
 my $file = shift;
 die "No file or does not exist\n" unless $file && -s $file;
 open( FQ, $file );
 my $max_lines = 100000;
 my $id        = <FQ>;
 die "File is not a fastq file:\n$id\n" unless ( $id =~ /^@/ );
 my ( @line, $l, $number, $counter );
 while ( my $seq = <FQ> ) {
  $counter++;
  my $qid    = <FQ>;
  my $qual   = <FQ>;
  my $nextid = <FQ>;
  @line = split( //, $qual );    # divide in chars
  for ( my $i = 0 ; $i <= $#line ; $i++ ) {    # for each char
   $number = ord( $line[$i] );    # get the number represented by the ascii char
     # check if it is sanger or illumina/solexa, based on the ASCII image at http://en.wikipedia.org/wiki/FASTQ_format#Encoding
   if ( $number > 75 ) {    # if solexa/illumina
    print "This file is solexa/illumina format\n"
      ;                     # print result to terminal and die
    close FQ;
    return 'illumina';
   }
  }
  last if $counter >= $max_lines;
 }
 close FQ;
 if ( $number < 59 ) {      # if sanger
  $id =~ /(\S+)\s*(\S*)/;
  my $description = $2;
  if ( $description && $description =~ /(\d)\:[A-Z]\:/ ) {
   print "This file is in Sanger quality but CASAVA 1.8 header format\n";
   return 'casava';
  }
  print "This file is sanger format\n";    # print result to terminal and die
  return 'sanger';
 }
 die "Cannot determine fastq format\n";
}

sub prepare_r_histogram() {
 my $file = shift;
 return unless $file && -s $file;
 open( R, ">$file.R" );
 if ( $file =~ /spec/ ) {
  print R 'data.25mer<-read.table(sep="\t",header=T,file="' . $file . '");';
  print R "\n" . 'postscript(file="' . $file . '.ps");';
  print R "\n"
    . 'plot(data.25mer$X1.kmer_frequency[1:10],data.25mer$X2.num_distinct_kmers[1:10],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$X1.kmer_frequency[1:20],data.25mer$X2.num_distinct_kmers[1:20],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$X1.kmer_frequency[1:30],data.25mer$X2.num_distinct_kmers[1:30],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$X1.kmer_frequency[3:30],data.25mer$X2.num_distinct_kmers[3:30],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$X1.kmer_frequency[5:30],data.25mer$X2.num_distinct_kmers[5:30],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$X1.kmer_frequency[3:50],data.25mer$X2.num_distinct_kmers[3:50],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$X1.kmer_frequency[3:80],data.25mer$X2.num_distinct_kmers[3:80],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$X1.kmer_frequency[3:100],data.25mer$X2.num_distinct_kmers[3:100],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$X1.kmer_frequency[3:150],data.25mer$X2.num_distinct_kmers[3:150],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$X1.kmer_frequency[3:200],data.25mer$X2.num_distinct_kmers[3:200],type="l",col="red");';
 }
 else {
  print R 'data.25mer<-read.table(sep="\t",header=F,file="' . $file . '");';
  print R "\n" . 'postscript(file="' . $file . '.ps");';
  print R "\n"
    . 'plot(data.25mer$V1[1:10],data.25mer$V2[1:10],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$V1[1:20],data.25mer$V2[1:20],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$V1[1:30],data.25mer$V2[1:30],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$V1[3:30],data.25mer$V2[3:30],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$V1[5:30],data.25mer$V2[5:30],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$V1[3:50],data.25mer$V2[3:50],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$V1[3:80],data.25mer$V2[3:80],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$V1[3:100],data.25mer$V2[3:100],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$V1[3:150],data.25mer$V2[3:150],type="l",col="red");';
  print R "\n"
    . 'plot(data.25mer$V1[3:200],data.25mer$V2[3:200],type="l",col="red");';
 }
 print R "\n" . 'dev.off();';
 close R;
 return "$file.R";
}

sub process_cmd {
 my ( $cmd, $dir ) = @_;
 print "CMD: $cmd\n" if $debug;
 chdir($dir) if $dir;
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  foreach my $thread (@threads) {
   $thread->join();
  }
  die "Error, cmd died with ret $ret\n";
 }
 chdir($cwd) if $dir;
 return;
}

sub create_fork() {
 my $cmd = shift;
 my $thread = threads->create( 'process_cmd', $cmd );
 sleep(1);
 return $thread;
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

sub do_genome_kmers() {
 my $infile = shift;
 die if !$infile || !-s $infile;
 my $format        = &check_fastq_format($infile);
 my $have_allpaths = `which FastqToFastbQualb`;
 my $have_meryl    = `which meryl`;
 if ( $have_allpaths && !$use_meryl ) {
  print "Using allpathsLG to estimate kmer distribution for K=25\n";
  my $is_mp_lib = ( "$infile" =~ /kb_/ ) ? int(1) : int(0);
  if ( $is_mp_lib == 1 ) {
   print
"Warning: I think this is an MP outie library. Creating FASTB as MP, will not be useful to allpaths if this is a PE innie library\n";
   if ( $format eq 'sanger' ) {
    &process_cmd(
"FastqToFastbQualb WRITE_QUALB=True REVERSE_READS=True FASTQ=$infile OUT_HEAD=$infile"
    ) unless -s "$infile.fastb";
   }
   else {
    &process_cmd(
"FastqToFastbQualb WRITE_QUALB=True REVERSE_READS=True PHRED_64=True FASTQ=$infile OUT_HEAD=$infile"
    ) unless -s "$infile.fastb";
   }

  }
  else {
   print
"Warning: Creating FASTB as PE data (innie), will not be useful to allpaths if this is a MP outie library\n";
   if ( $format eq 'sanger' ) {
    &process_cmd(
            "FastqToFastbQualb WRITE_QUALB=True FASTQ=$infile OUT_HEAD=$infile")
      unless -s "$infile.fastb";
   }
   else {
    &process_cmd(
"FastqToFastbQualb WRITE_QUALB=True PHRED_64=True FASTQ=$infile OUT_HEAD=$infile"
    ) unless -s "$infile.fastb";
   }
  }
  &process_cmd(
"KmerSpectrum K=25 HEAD=$infile MAX_MEMORY_GB=$kmer_ram NUM_THREADS=$cpus G_ESTIMATE=True PLOIDY=2 2>&1 | tee $infile.kmer.log"
  ) unless -s "$infile.25mer.kspec";
  if ( -s "$infile.25mer.kspec" ) {
   &process_cmd( 'search_replace.pl -s "^ " -i ' . "$infile.25mer.kspec" );
   &process_cmd(
              'search_replace.pl -s "^# " -d -i ' . "$infile.25mer.kspec.out" );
   &process_cmd(
         'search_replace.pl -s " " -r "	" -d -i ' . "$infile.25mer.kspec.out" );
   my $r_script = &prepare_r_histogram("$infile.25mer.kspec.out");
   &process_cmd("R --no-save  --no-restore -q -f $r_script ") if $r_script;
   &process_cmd("ps2pdf $infile.25mer.kspec.out.ps");
  }
  unlink("$infile.fastb") if $delete_fastb;
 }
 elsif ($have_meryl) {
  my $meryl_ram = ( $kmer_ram * 1024 ) . 'MB';
  print "Using meryl to estimate kmer distribution for K=25\n";
  my $ram = $meryl_ram;
  &process_cmd("meryl -B -C -m 25 -v -threads 2 -s $infile -o $infile")
    unless -s "$infile.mcdat";
  &process_cmd("meryl -Dh -s $infile > $infile.25mers")
    unless -s "$infile.25mers";
  my $r_script = &prepare_r_histogram("$infile.25mers");
  &process_cmd("R --no-save  --no-restore -q -f $r_script ") if $r_script;
  &process_cmd("ps2pdf $infile.25mers.ps");
 }
 else {
  print
"Warning: neither allpaths-LG nor meryl found in path. Will not create a kmer histogram\n";
 }
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

sub samtools_version_check() {
 my $samtools_exec = shift;
 my @version_lines = `$samtools_exec 2>&1`;
 foreach my $ln (@version_lines) {
  if ( $ln =~ /^Version:\s+\d+\.(\d+).(\d+)/i ) {
   die "Samtools version 1.19+ is needed\n" unless $1 >= 1 && $2 >= 19;

   #print "Good: Samtools version 1.19+ found\n";
  }
 }
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
