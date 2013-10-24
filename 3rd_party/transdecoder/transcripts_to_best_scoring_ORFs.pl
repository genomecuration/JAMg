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

### cd-hit-est : http://www.bioinformatics.org/downloads/index.php/cd-hit/cd-hit-v4.5.4-2011-03-07.tgz
my $cd_hit_est_exec = `which cd-hit-est`;
chomp($cd_hit_est_exec) if $cd_hit_est_exec;

my ($transcripts_file,$train_file,$prepare_pfam_only);

my $min_prot_length = 100;
my $genetic_code;

my $top_ORFs_train = 500;

my $TOP_STRAND_ONLY = 0;

my $help;
my $workdir;
my $verbose;
my $search_pfam = "";
my ($reuse,$pfam_out);
my $CPU = 2;
my $RETAIN_LONG_ORFS = 900;


my $usage =  <<_EOH_;

######################################## Options ###################################################################################
#
# ** Required:
#
# -t <string>                            transcripts.fasta
#
# ** Optional:
# 
# --reuse                                If this option is given, any existing files are not overwritten but reused
#
# --train <string>                       FASTA file with ORFs to train Markov Mod for protein identification; otherwise 
#                                        longest non-redundant ORFs used
#
# -m <int>                               minimum protein length (default: 100)
#
# --search_pfam <string>                 /path/to/pfam_db.hmm to search 
#                                        using hmmscan (which should be accessible via your PATH setting)
# 
# --pfam_out <string>                    You can also pre-run the pfam searches if --reuse is set. In that case, 
#                                        --pfam_out is the output of hhmscan --tblout using --noali --cut_nc --acc --notextw
#
# --prepare_pfam                         Prepare data for PFAM search and then quit (for running PFAM on HPC/computing cluster)
#
# --cd_hit_est <string>                  Optionally and only if it is not in your path, the full path to the cd-hit-est executable
#
# --workdir                              Force temporary output directory to this directory (e.g. if --reuse is needed)
#
# -G <string>                            genetic code (default: universal, options: Euplotes, Tetrahymena, Candida, Acetabularia)
#
#
# -h                                     print this option menu and quit
# -v                                     verbose
#
# -S                                     strand-specific (only analyzes top strand)
# -T <int>                               If no --train, top longest ORFs to train Markov Model (hexamer stats) (default: 500)
#
# --CPU <int>                            number of threads to use; (default: 2)
#
# --retain_long_orfs <int>               retain all ORFs found that are of minimum length in nucleotides (default: 900 bp => 300aa)
#
# --quiet                                send stderr to /dev/null
#
####################################################################################################################################

_EOH_

    ;


my $QUIET = 0;

&GetOptions( 't=s' => \$transcripts_file,
             'train:s' => \$train_file,
             'm=i' => \$min_prot_length,
             'G=s' => \$genetic_code,
             'h' => \$help,
             'v' => \$verbose,
             'S' => \$TOP_STRAND_ONLY, 
             'T=i' => \$top_ORFs_train,
             'CPU=i' => \$CPU,
             'search_pfam=s' => \$search_pfam,
             'reuse' => \$reuse,
             'workdir:s' => \$workdir,
             'pfam_out:s' => \$pfam_out,
             'retain_long_orfs=i' => \$RETAIN_LONG_ORFS,
             'cd_hit_est=s' => \$cd_hit_est_exec,
             'prepare_pfam' => \$prepare_pfam_only,
             '--quiet' => \$QUIET,
             );



if ($help) {
    die $usage;
}

if (@ARGV) {
    die "Error, don't understand options: @ARGV";
}

$|++;

our $SEE = $verbose;

unless ($transcripts_file) {
    die "$usage\n";
}




&check_for_pfam_execs if ($search_pfam);   


$workdir = "transdecoder.tmp.$$" unless $workdir;
mkdir($workdir) unless -d $workdir;
die "Error, cannot mkdir $workdir" unless -d $workdir;

if ($genetic_code) {
    &Nuc_translator::use_specified_genetic_code($genetic_code);
}

my $number_of_peps = int(0);
my $prefix = "$workdir/longest_orfs";
my $cds_file = "$prefix.cds";
my $gff3_file = "$prefix.gff3";
my $pep_file = "$prefix.pep";
die "PFAM out file not found: $pfam_out\n" if $pfam_out && !-s $pfam_out;
$pfam_out = basename($transcripts_file) . ".transdecoder.pfam.dat" unless $pfam_out && -s $pfam_out;

my (%orf_lengths,$cmd);

if ($reuse && -s $pep_file && -s $cds_file && -s $gff3_file) {
    $number_of_peps = `grep -c "^>" $pep_file` ;
    chomp($number_of_peps);
    open (IN,$cds_file);
    while (my $ln=<IN>){
        next unless $ln=~/^>(\S+).+len:(\d+)/;
        $orf_lengths{$1}=$2 if $1 && $2;
    }
    close IN;
}
else {
	open (PEP, ">$pep_file") or die $!;
	open (CDS, ">$cds_file") or die $!; 
	open (GFF, ">$gff3_file") or die $!;
	
	
	my $counter = 0;
	
	my $fasta_reader = new Fasta_reader($transcripts_file);
	while (my $seq_obj = $fasta_reader->next()) {
		
		my $acc = $seq_obj->get_accession();
		my $sequence = $seq_obj->get_sequence();
		
		my $longest_orf_finder = new Longest_orf();
		$longest_orf_finder->allow_5prime_partials();
		$longest_orf_finder->allow_3prime_partials();
		
	    if ($TOP_STRAND_ONLY) {
			$longest_orf_finder->forward_strand_only();
		}
		
		my @orf_structs = $longest_orf_finder->capture_all_ORFs($sequence);
		
		@orf_structs = reverse sort {$a->{length}<=>$b->{length}} @orf_structs;
		
        while (@orf_structs) {
            my $orf = shift @orf_structs;
            
            my $start = $orf->{start};
            my $stop = $orf->{stop};
            
            my $length = int((abs($start-$stop)+1)/3); #int($orf->{length}/3);
            my $orient = $orf->{orient};
            my $protein = $orf->{protein};            
            
            ##################################
            # adjust for boundary conditions, since starts and stops run off the ends of the sequences at partial codons
            #################################
            
            # adjust at 3' end
            if ($stop > length($sequence)) {
                $stop -= 3;
            }
            if ($start > length($sequence)) {
                $start -= 3;
            }
            
            # adjust at 5' end
            if ($stop < 1) {
                $stop += 3;
            }
            if ($start < 1) {
                $start += 3;
            }
            
                        
            if ($length < $min_prot_length) { next; }
            
            my $cds_coords_href = { $start => $stop };
            my $exon_coords_href = ($start < $stop) ? { 1 => length($sequence) } : { length($sequence) => 1 };
            
            my $gene_obj = new Gene_obj();
            
            $counter++;
            $gene_obj->populate_gene_object($cds_coords_href, $exon_coords_href);
            $gene_obj->{asmbl_id} = $acc;
            
            my $model_id = "$acc|m.$counter";
            my $gene_id = "$acc|g.$counter";
            
            $gene_obj->{TU_feat_name} = $gene_id;
            $gene_obj->{Model_feat_name} = $model_id;

            
            my $cds = $gene_obj->create_CDS_sequence(\$sequence);
            
            unless ($cds) {
                die "Error, no CDS for gene: " . Dumper($cds_coords_href) . Dumper($exon_coords_href);
            }

            my $got_start = 0;
            my $got_stop = 0;
            if ($protein =~ /^M/) {
                $got_start = 1;
            } 
            if ($protein =~ /\*$/) {
                $got_stop = 1;
            }
            
            my $prot_type = "";
            if ($got_start && $got_stop) {
                $prot_type = "complete";
            } elsif ($got_start) {
                $prot_type = "3prime_partial";
            } elsif ($got_stop) {
                $prot_type = "5prime_partial";
            } else {
                $prot_type = "internal";
            }
            
            $gene_obj->{com_name} = "ORF $gene_id $model_id type:$prot_type len:$length ($orient)";            
            
            print PEP ">$model_id $gene_id type:$prot_type len:$length $acc:$start-$stop($orient)\n$protein\n";
            $number_of_peps++;
            
            print CDS ">$model_id $gene_id type:$prot_type len:$length\n$cds\n";
            
            print GFF $gene_obj->to_GFF3_format(source => "transdecoder") . "\n";
            

            $orf_lengths{$model_id} = length($cds);
            


        }
	}

    close PEP;
    close CDS;
    close GFF;
    
}

## Train a Markov model based on user-provided file or longest candidate CDS sequences, score all candidates, and select the final set.

my $top_cds_file = $train_file && -s $train_file ? $train_file : "$cds_file.top_${top_ORFs_train}_longest";
if (!-s $top_cds_file) {
    # get longest entries
    my $cmd = "$UTIL_DIR/get_top_longest_fasta_entries.pl $cds_file $top_ORFs_train > $top_cds_file";
    
    unless ($reuse && -s $top_cds_file){
        if ($cd_hit_est_exec){
            # to speed things up only check for redundancy up to 4x the number of entries we want
            my $red_num = $top_ORFs_train * 4 ;
            &process_cmd("$UTIL_DIR/get_top_longest_fasta_entries.pl $cds_file $red_num > $workdir/redundant_top");
            &process_cmd("$cd_hit_est_exec -r 1 -i $workdir/redundant_top -o $workdir/redundant_top.nr90 -M 0 -T $CPU >/dev/null 2>/dev/null");
            &process_cmd("$UTIL_DIR/get_top_longest_fasta_entries.pl $workdir/redundant_top.nr90 $top_ORFs_train > $top_cds_file");
            unlink("$workdir/redundant_top");
            unlink("$workdir/redundant_top.nr90");
            unlink("$workdir/redundant_top.nr90.bak.clstr");
        }
        else {
            &process_cmd($cmd);
        }
    }
}

$cmd = "$UTIL_DIR/compute_base_probs.pl $transcripts_file $TOP_STRAND_ONLY > $workdir/base_freqs.dat";
&process_cmd($cmd) unless $reuse && -s "$workdir/base_freqs.dat";


# get hexamer scores
#$cmd = "$UTIL_DIR/seq_n_background_to_logliklihood_vals.pl $top_cds_file $transcripts_file.random > hexamer.scores";
#&process_cmd($cmd) unless ($reuse && -s "hexamer.scores");

$cmd = "$UTIL_DIR/seq_n_baseprobs_to_logliklihood_vals.pl $top_cds_file $workdir/base_freqs.dat > $workdir/hexamer.scores";
&process_cmd($cmd) unless $reuse && -s "$workdir/hexamer.scores";


# score all cds entries
$cmd = "$UTIL_DIR/score_CDS_liklihood_all_6_frames.pl $cds_file $workdir/hexamer.scores > $cds_file.scores";
&process_cmd($cmd) unless ($reuse && -s "$cds_file.scores");

# run pfam
my %has_pfam_hit;
if ($search_pfam) {
    my $pfam_cmd = "hmmscan -o /dev/null --cpu 1 --noali --cut_nc --acc --notextw --tblout /dev/stdout --domtblout /dev/stdout $search_pfam - ";    
#   my $pfam_cmd = "hmmscan -o /dev/null --cpu 1 --noali --cut_nc --acc --notextw --tblout $workdir/INFILE.pfam.out --domtblout $workdir/INFILE.pfam.domtbl.out $search_pfam - ";
# this is what we could use for hhblits need $hhblits_exec and $hhblits_db (PFAM)
# my $hhblits_cmd = "$hhblits_exec -i INFILE -d $hhblits_db -oa3m BASENAME.a3m -n 2 -mact 0.5 >/dev/null "

    if ($prepare_pfam_only) {
        my $parafly_cmd_file = &multithread($pep_file,$pfam_cmd);
	my $parafly_cmd = "ParaFly -CPU $CPU -c $parafly_cmd_file --failed $parafly_cmd_file.failed -v";
        print "We have prepared the $parafly_cmd_file command file for you to run hmmscan separately (e.g. on a cluster).\n";
        print "Example for a single node with $CPU CPUs :\n\t$parafly_cmd\n\n";
        print "After your PFAM searches are complete, then concatanate all the out.db files using this command into $pep_file.pfam.out.\n";
        print "cat $workdir/*.out.db" 
              . '|tr -d \'\000\' ' 
              . "|grep -v '^#' > $pep_file.pfam.out\n\n";
        print "Then in order to restart transdecoder use the following command (along with any other options you want):\n";
        print "\t$0 -t $transcripts_file --search_pfam $search_pfam --pfam_out $pep_file.pfam.out --reuse --workdir $workdir\n\n";
        exit(0);
    }
    
    print "Processing with PFAM HMM searches...\n";
    
    if ($reuse && -s $pfam_out){
    #    open (CHECK,$pfam_out);
    #    my $check = <CHECK>;
    #    close CHECK;
    #    die "You've asked to use the $pfam_out file as hhsearch --tblout output, but it doesn't look like the right file!\n" unless $check && $check=~/full sequence/;
    #    my $check_number = `grep -c '^# target name' $pfam_out`;
    #    chomp($check_number);
    #    die "The number of sequences in $pfam_out do not equal the number of sequences in the protein file $pep_file\nDelete it and use --rerun and --workdir $workdir to resume.\n" unless $check_number == $number_of_peps;
    }else {
        my $parafly_cmd_file = &multithread($pep_file,$pfam_cmd);
	my $parafly_cmd = "ParaFly -CPU $CPU -c $parafly_cmd_file --failed $parafly_cmd_file.failed -v";
	&process_cmd("$parafly_cmd");
        if (-s "$parafly_cmd_file.failed"){
           die "Some sequences failed to be searched against PFAM. Please resolve the situation (see $parafly_cmd_file), delete $parafly_cmd_file and use --rerun and --workdir $workdir to resume.\n";
	}
        &process_cmd("cat $workdir/*.out.db" . '|tr -d "\000" | grep -v "^#" > ' . $pfam_out);
        &process_cmd("find $workdir -name '*.out.db' -delete");
        &process_cmd("find $workdir -name '*.out.db.idx' -delete");
        &process_cmd("find $workdir -name 'parafly.sh*' -delete");
            unless (-s $pfam_out) {
                die "Error, pfam results were not properly concatenated into file: $pfam_out";
            }
            
    }

    # parse PFAM output and split into DOM and TBL files. DOM file is reused by trinnotate
    open (my $fh, $pfam_out) or die "Error, cannot open file: $pfam_out";
    open (DOMOUT,">$pfam_out.domtbl");
    open (TBLOUT,">$pfam_out.tbl");
    while (my $ln=<$fh>) {
        next if $ln=~/^#/;
        my @x = split(/\s+/,$ln);
	next unless $x[2];
	# if this is the output of the domtbl part, then skip it and write it into another file.
	if ($x[2]=~/^\d+$/ && $x[5]=~/^\d+$/){
		print DOMOUT $ln;
		next;
	}
        print TBLOUT $ln;
        my $orf_acc = $x[2];
        $has_pfam_hit{$orf_acc} = 1;
    }
    close $fh;
    close DOMOUT;
    close TBLOUT;
}

# get accs for best entries
my $acc_file = "$cds_file.scores.selected";
{
	open (my $ofh, ">$acc_file") or die "Error, cannot write to $acc_file";
	open (my $ifh, "$cds_file.scores") or die "Error, cannot open file $cds_file.scores";
	while (<$ifh>) {
		chomp;
		my ($acc, @scores) = split(/\t/);
		
		my $score_1 = shift @scores;
		my $max_score_other_frame = max(@scores);
		if ($has_pfam_hit{$acc} 
            || 
            $orf_lengths{$acc} >= $RETAIN_LONG_ORFS
            ||
            ($score_1 > 0 && $score_1 > $max_score_other_frame)
            ) { 
			print $ofh "$acc\n";
            
            if ($has_pfam_hit{$acc}) {
                print STDERR "-$acc flagged as having a pfam domain.\n" if $verbose;
            }
            
        }
	}
	close $ifh;
	close $ofh;
}

# index the current gff file:
$cmd = "$UTIL_DIR/index_gff3_files_by_isoform.pl $gff3_file";
&process_cmd($cmd);

# retrieve the best entries:
$cmd = "$UTIL_DIR/gene_list_to_gff.pl $acc_file $gff3_file.inx > $cds_file.best_candidates.gff3";
&process_cmd($cmd);

{
    my $final_output_prefix = basename($transcripts_file) . ".transdecoder";
    
    # exclude shadow orfs (smaller orfs in different reading frame that are eclipsed by longer orfs)
    $cmd = "$UTIL_DIR/remove_eclipsed_ORFs.pl $cds_file.best_candidates.gff3 > $final_output_prefix.gff3";
    &process_cmd($cmd);
    


    ## write final outputs:
    
    ## make a BED file for viewing in IGV
    my $gff3_file = "$final_output_prefix.gff3";
    my $bed_file = $gff3_file;
    $bed_file =~ s/\.gff3$/\.bed/;
    $cmd = "$UTIL_DIR/gff3_file_to_bed.pl $gff3_file > $bed_file";
    &process_cmd($cmd);
    
    
    # make a peptide file:
    my $best_pep_file = $gff3_file;
    $best_pep_file =~ s/\.gff3$/\.pep/;
    $cmd = "$UTIL_DIR/gff3_file_to_proteins.pl $gff3_file $transcripts_file > $best_pep_file";
    &process_cmd($cmd);



    # make a CDS file:
    my $best_cds_file = $best_pep_file;
    $best_cds_file =~ s/\.pep$/\.cds/;
    $cmd = "$UTIL_DIR/gff3_file_to_proteins.pl $gff3_file $transcripts_file CDS > $best_cds_file";
    &process_cmd($cmd);
    
}

print STDERR "transdecoder is finished.\n";


exit(0);


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
   my $sequence_number = `wc -l < $fasta_file.db.idx`;chomp($sequence_number);
   die unless $sequence_number;
   return $sequence_number;
}


# MPIrun might require expert users here.... so just use parafly and be done with it (experts can run mpirun directly)
sub multithread(){
    my ($protein_file,$cmd) = @_;
    print "Preparing CMD for $number_of_peps files\n";
    my @fasta_files = &partition_transcript_db($protein_file);
    my $cmd_file = "$workdir/parafly.sh";
    unlink($cmd_file.'.completed');
    unlink($cmd_file.'.failed');
    open(OUT,">$cmd_file") || die; 
    foreach my $fasta_file (@fasta_files){
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
    close OUT;
    return "$cmd_file";
}


sub partition_transcript_db {
    my $transcript_db = shift;
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
    my @utils = qw(hmmscan ParaFly ffindex_apply_mpi ffindex_from_fasta);

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



