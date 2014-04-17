#!/usr/bin/env perl

use FindBin;
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);
use Data::Dumper;
use List::Util qw (min max);
use File::Basename;

use lib ("$FindBin::RealBin/PerlLib");

use POSIX qw(ceil);
use Gene_obj;
use Nuc_translator;
use Fasta_reader;
use Longest_orf;



my $UTIL_DIR = "$FindBin::RealBin/util";

my $help;
my $workdir;
my $verbose;
my ($reuse,$pfam_out);
my $CPU = 2;

my $usage =  <<_EOH_;

######################################## Options ###################################################################################
#
###############
# ** Required:
###############

# --pep <string>                         peptide files
#
# --pfam_db <string>                 /path/to/pfam_db.hmm to search 
#                                        using hmmscan (which should be accessible via your PATH setting)
#
################
# ** Optional:
################ 
#
#
# --reuse                                If this option is given, any existing files are not overwritten but reused
#
# 
# --pfam_out|o <string>                    You can also pre-run the pfam searches if --reuse is set. In that case, 
#                                        --pfam_out is the output of hhmscan --tblout using --noali --cut_nc --acc --notextw
#
# --prepare_pfam                         Prepare data for PFAM search and then quit (for running PFAM on HPC/computing cluster)
#
# --workdir                              Force temporary output directory to this directory (e.g. if --reuse is needed)
#
#
# -h                                     print this option menu and quit
# -v                                     verbose
#
# --CPU <int>                            number of threads to use; (default: 2)
#
# --MPI                                  use MPI (via ffindex_apply_mpi)
#
# --quiet                                send stderr to /dev/null
#
#
####################################################################################################################################

_EOH_

    ;


my $QUIET = 0;
my $pep_file;
my $search_pfam;
my $prepare_pfam_only = 0;
my $MPI = 0;

&GetOptions( 'pep=s' => \$pep_file,
             'h' => \$help,
             'v' => \$verbose,
             'CPU=i' => \$CPU,
             'pfam_db=s' => \$search_pfam,
             'reuse' => \$reuse,
             'workdir:s' => \$workdir,
             'pfam_out|o=s' => \$pfam_out,
             'prepare_pfam' => \$prepare_pfam_only,
             'quiet' => \$QUIET,
             'MPI' => \$MPI,
             );



if ($help) {
    die $usage;
}

if (@ARGV) {
    die "Error, don't understand options: @ARGV";
}

$|++;

our $SEE = $verbose;

unless ($pep_file && $search_pfam) {
    die "$usage\n";
}

&check_for_pfam_execs if ($search_pfam);   

$workdir = "transdecoder.tmp.$$" unless $workdir;
mkdir($workdir) unless -d $workdir;
die "Error, cannot mkdir $workdir" unless -d $workdir;

unless ($pfam_out) {
    $pfam_out = basename($pep_file) . ".transdecoder.pfam.dat";
}


main: {


    #   my $pfam_cmd = "hmmscan -o /dev/null --cpu 1 --noali --cut_nc --acc --notextw --tblout $workdir/INFILE.pfam.out --domtblout $workdir/INFILE.pfam.domtbl.out $search_pfam - ";
    # this is what we could use for hhblits need $hhblits_exec and $hhblits_db (PFAM)
    # my $hhblits_cmd = "$hhblits_exec -i INFILE -d $hhblits_db -oa3m BASENAME.a3m -n 2 -mact 0.5 >/dev/null "

    
    my $parafly_cmd_file = &multithread($pep_file);
    
    my $parafly_cmd = "ParaFly -CPU $CPU -c $parafly_cmd_file --failed $parafly_cmd_file.failed ";
    if ($QUIET) { 
        $parafly_cmd .= " -v ";
    }
    else {
        # a little more verbose
        $parafly_cmd .= " -vv ";
    }
    
    if ($prepare_pfam_only) {
        
        
        print "We have prepared the $parafly_cmd_file command file for you to run hmmscan separately (e.g. on a cluster).\n";
        print "Example for a single node with $CPU CPUs :\n\t$parafly_cmd\n\n";
        print "After your PFAM searches are complete, then concatanate all the out.db files using this command into $pep_file.pfam.out.\n";
        print "cat $workdir/*.out.db" 
              . '|tr -d \'\000\' ' 
              . "|grep -v '^#' > $pep_file.pfam.out\n\n";
        print "Then in order to restart transdecoder use the following command (along with any other options you want):\n";
        print "\t$0 --pep $pep_file --search_pfam $search_pfam --pfam_out $pep_file.pfam.out --reuse --workdir $workdir\n\n";
        
        exit(0);
    }
    
    print "Processing with PFAM HMM searches...\n";
    
    unless ($reuse && -s $pfam_out){
                
        &process_cmd("$parafly_cmd");
        
        if (-s "$parafly_cmd_file.failed"){
            die "Some sequences failed to be searched against PFAM. Please resolve the situation (see $parafly_cmd_file), delete $parafly_cmd_file and use --rerun and --workdir $workdir to resume.\n";
        }
        
        if ($MPI) {

            &process_cmd("cat $workdir/*.out.db" . '|tr -d "\000" | grep -v "^#" > ' . $pfam_out);
        
        
            # lets not purge these just yet...  debugging a clobber problem.
            
            #&process_cmd("find $workdir -name '*.out.db' -delete");
            #&process_cmd("find $workdir -name '*.out.db.idx' -delete");
            #&process_cmd("find $workdir -name 'parafly.sh*' -delete");
           
            unless (-s $pfam_out) {
                die "Error, pfam results were not properly concatenated into file: $pfam_out";
            }
        
            ## just capturing the domain table
            &process_cmd("ln -s $pfam_out $pfam_out.domtbl");
            
            
            # parse PFAM output and split into DOM and TBL files. DOM file is reused by trinnotate
            #open (my $fh, $pfam_out) or die "Error, cannot open file: $pfam_out";
            #open (DOMOUT,">$pfam_out.domtbl");
            #open (TBLOUT,">$pfam_out.tbl");
            #while (my $ln=<$fh>) {
            #    next if $ln=~/^\#/;
            #    my @x = split(/\s+/,$ln);
            #    next unless $x[2];
            #    # if this is the output of the domtbl part, then skip it and write it into another file.
            #    if ($x[2]=~/^\d+$/ && $x[5]=~/^\d+$/){
            #        print DOMOUT $ln;
            #        next;
            #    }
            #    print TBLOUT $ln;
            #}
            #close $fh;
            #close DOMOUT;
            #close TBLOUT;
        }
        else {
            ## Parse regular pfam table output files:
            &process_cmd("find $workdir/ -name '*.fa.tbl' -exec cat {} \\\; | egrep -v '^\#' > $pfam_out.tbl");
            &process_cmd("find $workdir/ -name '*.fa.domtbl' -exec cat {} \\\; | egrep -v '^\#' > $pfam_out"); # just the domain hits.
            
            &process_cmd("ln -s $pfam_out $pfam_out.domtbl"); # so just like under MPI
            
        }
    }
    
    print STDERR "PFAM SEARCH DONE.\n\n";
    
    exit(0);

}

####
sub process_cmd {
	my ($cmd) = @_;

	print "CMD: $cmd\n";
	my $ret = system($cmd);

	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}
	
	return;

}

sub index_fasta(){
   # this where ffindex would really help
   my $fasta_file = shift;
   &process_cmd("ffindex_from_fasta -s $fasta_file.db $fasta_file.db.idx $fasta_file") unless -s "$fasta_file.db";
   my $sequence_number = `wc -l < $fasta_file.db.idx`;
   chomp($sequence_number);
   die unless $sequence_number;
   return $sequence_number;
}


# MPIrun might require expert users here.... so just use parafly and be done with it (experts can run mpirun directly)
sub multithread(){
    my ($protein_file) = @_;
    print STDERR "Partitioning fasta file $protein_file\n";
    my @fasta_files = &partition_transcript_db($protein_file);
    my $cmd_file = "$workdir/parafly.sh";
    unlink($cmd_file.'.completed');
    unlink($cmd_file.'.failed');
    open(OUT,">$cmd_file") || die; 
    foreach my $fasta_file (@fasta_files){

        if ($MPI) {

            #my $cmd = "hmmscan -o /dev/null --cpu 1 --noali --cut_nc --acc --notextw --tblout /dev/stdout --domtblout /dev/stdout $search_pfam - ";    
            # JUST CAPTURING THE DOMAIN TABLE UNDER MPI MODE
            my $cmd = "hmmscan -o /dev/null --cpu 1 --noali --cut_nc --acc --notextw  --domtblout /dev/stdout $search_pfam - ";    

            my $sequence_number = &index_fasta($fasta_file);
            my $pfam_out = "$fasta_file.out";
            my $ffidx = "ffindex_apply_mpi -d $pfam_out.db -i $pfam_out.idx $fasta_file.db $fasta_file.db.idx -- ";
            unlink($pfam_out);
            unlink($pfam_out.'.db');
            unlink($pfam_out.'.idx');
            
            if ($QUIET) {
                $ffidx .= " 2>/dev/null";
            }
            
            print OUT $ffidx . " $cmd \n";
        }
        else {
            ## use ParaFly / hmmscan w/o ffindex
            my $cmd = "hmmscan -o /dev/null --cpu 1 --noali --cut_nc --acc --notextw --tblout $fasta_file.tbl --domtblout $fasta_file.domtbl $search_pfam $fasta_file ";
            print OUT "$cmd\n";
        }
        
    }
    close OUT;
    return "$cmd_file";
}


sub partition_transcript_db {
    my $transcript_db = shift;
    my $number_of_peps = `grep '>' $transcript_db | wc -l `;
    chomp $number_of_peps;
    
    my $seqs_per_partition = ceil($number_of_peps/$CPU);
    $seqs_per_partition = 1 if $seqs_per_partition < 1;
    $seqs_per_partition = $seqs_per_partition < 5000 ?  $seqs_per_partition : 5000  ;
    my @files;
    my $fasta_reader = new Fasta_reader($transcript_db);
    my $partition_counter = 0;
    my $counter = 0;
    my $ofh;
    while (my $seq_obj = $fasta_reader->next()) {
            my $fasta_entry = $seq_obj->get_FASTA_format();
	    $fasta_entry=~s/[\*\s]+$//; #strip stop codon/empty space
	    $fasta_entry.="\n";
            if ($counter % $seqs_per_partition == 0) {
                close $ofh if $ofh;
                $partition_counter++;
                my $outfile = "$workdir/partition.$counter.fa";
                open ($ofh, ">$outfile") or die "Error, cannot write to outfile: $outfile";
                push (@files, $outfile);
            }
            print $ofh $fasta_entry;
            $counter++;
    }
    close $ofh if $ofh;
    return(@files);
}

sub check_for_pfam_execs(){
    
    $ENV{PATH} .= ":$UTIL_DIR/bin"; # now can find 3rd party tools in PATH setting
    
    $ENV{LD_LIBRARY_PATH} .= ":$UTIL_DIR/lib64/";
    
    die "Error, cannot locate pfam database at: $search_pfam"  unless (-s $search_pfam);
    my @utils = qw(hmmscan ParaFly);

    if ($MPI) {
        push (@utils, qw(ffindex_apply_mpi ffindex_from_fasta));
    }
    
    my $missing_flag = 0;
    foreach my $util (@utils) {
        my $path = `which $util`;
        unless ($path =~ /\w/) {
            print STDERR "ERROR, cannot locate tool: $util in PATH setting.\n";
            $missing_flag++;
        }
    }

    die "Fatal, cannot find $missing_flag tools required.\n" if ($missing_flag);
    
}



