#!/usr/bin/perl -w

use strict;
use Getopt::Std;
use lib '.';
use GTF;

use vars qw($opt_c $opt_e $opt_f $opt_p $opt_s $opt_t $opt_k $opt_l $opt_m $opt_h);
getopts('ce:fpst:klmh');
my $usage = "usage: $0 [-fsmbB] [-t tx output filename] <gtf file> [sequence file]
Options:
 -t <file>: output transcript file
 -f: create a fixed gtf file (This may not be possible.  
     Always check the \"fixed\" file) 
 -e <count>: sets the maximum number of detailed error messages to return per 
     error to <count> (default is 5).
 -s: output list of inframe stop genes.
 -c: suppress warnings about missing start/stop
 -p: suppress warnings about bad splice site sequence
 -k: output a list of bad genes for \"super-clean\" training set
 -l: output a list of bad genes for training applications
 -m: output a list of bad genes for evaluation purposes
 -h: Display explanation on how to download GTF files from UCSC and clean
     them using this script. This is how to generate evaluation (RefSeq) sets!


";

if($opt_h){
    my$info=get_info();
    print $info;
    exit(0);
}
die $usage unless ((@ARGV == 1)||(@ARGV == 2));

my ($filename,$seqname) = @ARGV;
if($opt_t){
    unless(defined($seqname)){
	die "Must give a sequence file when you give the -t option.\n";
    }
    open(TX,">$opt_t");
}
else{
    open(TX,">/dev/null");
}
my $fix_file = $opt_f;
my $fix_filename;
my $suppress = [];
if($opt_c){
    $$suppress[17] = 1;
    $$suppress[19] = 1;
}
if($opt_p){
    $$suppress[26] = 1;
    $$suppress[27] = 1;
}
my $bad_genes = [];
if($opt_m){
	$$bad_genes[24] = 1;
	$$bad_genes[25] = 1;
	$$bad_genes[26] = 1;
	$$bad_genes[27] = 1;
	$$bad_genes[28] = 1;
	$$bad_genes[29] = 1;
	$$bad_genes[37] = 1;
	$$bad_genes[38] = 1;
	$$bad_genes[41] = 1;
	$$bad_genes[45] = 1;
	$$bad_genes[48] = 1;
}
if($opt_l){
	$$bad_genes[16] = 1;
	$$bad_genes[18] = 1;
	$$bad_genes[21] = 1;
	$$bad_genes[22] = 1;
	$$bad_genes[24] = 1;
	$$bad_genes[25] = 1;
	$$bad_genes[26] = 1;
	$$bad_genes[27] = 1;
	$$bad_genes[28] = 1;
	$$bad_genes[29] = 1;
	$$bad_genes[33] = 1;
	$$bad_genes[34] = 1;
	$$bad_genes[35] = 1;
	$$bad_genes[36] = 1;
	$$bad_genes[37] = 1;
	$$bad_genes[38] = 1;
	$$bad_genes[41] = 1;
	$$bad_genes[48] = 1;
}
if($opt_k){
	$$bad_genes[15] = 1;
	$$bad_genes[16] = 1;
	$$bad_genes[17] = 1;
	$$bad_genes[18] = 1;
	$$bad_genes[19] = 1;
	$$bad_genes[21] = 1;
	$$bad_genes[22] = 1;
	$$bad_genes[24] = 1;
	$$bad_genes[25] = 1;
	$$bad_genes[26] = 1;
	$$bad_genes[27] = 1;
	$$bad_genes[28] = 1;
	$$bad_genes[29] = 1;
	$$bad_genes[33] = 1;
	$$bad_genes[34] = 1;
	$$bad_genes[35] = 1;
	$$bad_genes[36] = 1;
	$$bad_genes[37] = 1;
	$$bad_genes[38] = 1;
	$$bad_genes[41] = 1;
	$$bad_genes[45] = 1;
	$$bad_genes[46] = 1;
	$$bad_genes[47] = 1;
	$$bad_genes[48] = 1;
	$$bad_genes[50] = 1;
}
my $info = {gtf_filename  => $filename,
	    warning_fh    => \*STDOUT,
	    fix_gtf       => $fix_file,
	    seq_filename  => $seqname,
	    inframe_stops => $opt_s,
	    tx_out_fh     => \*TX,
	    warning_skips => $suppress,
	    bad_list      => $bad_genes};
if($opt_e){
    if($opt_e =~ /^(\d+)$/){
	$info->{detailed_error_count} = $opt_e;
    }
    else{
	print STDERR "Bad value, $opt_e, for -e flag.  Should be an integer.\n";
    }
}

my $gtf = GTF::new($info);

close(TX);
if($fix_file){
    if($filename =~ /(\S*).g[t,f]f/){
	$fix_filename = "$1.fixed.gtf";
    }
    else{
	$fix_filename = "$filename.fixed.gtf";
    }
    if(open(FIX, ">$fix_filename")){
	$gtf->output_gtf_file(\*FIX);
    }
    else{
	$fix_file = 0;
	print STDERR "Could not open $fix_filename for output.\n";
	print STDERR "Will not create fixed gtf file.\n";
    }
}

exit;

  
sub get_info{
	my$info="
To generate an evaluation set, three steps are needed:
1. Download the annotation file from UCSC to the proper /bio/db/ directory
2. Convert the file to GTF and split per chromosome
3. Clean the set using this program

STEP 1. Downloading
Put the files in the correct /bio/db directory. The convention is:
/bio/db/<SPECIES>/assembly/<assembly id>/annotation/<downloaded set_version>
For Refseq, this would be:
/bio/db/Homo_sapiens/assembly/hg17/annotation/refseq_v1
Use ftp for downloading the annotation for the latest genome build:
>ftp hgdownload.cse.ucsc.edu 
>cd goldenPath/currentGenomes/<species>/database
>get <file>
Eg for getting human RefSeqs:
>cd goldenPath/currentGenomes/Homo_sapiens/database
>get refGene.gtf.gz

If you're unsure of the filename you can find it in UCSC's Table Browser
(http://genome.ucsc.edu/cgi-bin/hgTables), select the Track (RefSeq)
and see the name that pops up under Table.

STEP 2. Converting
Unzip the file, convert to gtf, split per chromosome
>gunzip <file>
>/bio/bin/ucsc2gtf.pl <file> > <file.gtf>
This command may generate some \"skipping <gene>\" error messages.
In those cases, no CDS can be found.
>/bio/bin/divide_gtfs_in_chrs.pl <file.gtf>
...patience...

STEP 3. Cleaning
First, create a badlist containing genes with inframe stops,
reading frame changes and other incorrectible problems. Then,
remove those genes from the files. After that, merge any 
overlapping transcripts into one gene model and remove all
transcripts that are identical to other transcripts:
The first step needs the chromosome sequence. Get this
from /bio/db/<SPECIES>/assembly/<assembly id>/chr_seq/chrN.fa
For every chromosome:
>/bio/bin/validate_gtf.pl -m chrN.gtf /bio/db/<SPECIES>/assembly/<assembly id>/chr_seq/chrN.fa > chrN.badlist.txt 
run this on the queue!
>/bio/bin/filter_badlist.pl chrN.gtf chrN.badlist.txt > chrN.filtered.gtf
>/bio/bin/merge_gtf_transcripts.py chrN.filtered.gtf > chrN.eval.gtf
The last two commands can be run locally

Cleanup: cat all *.badlist.txt into a file called Badlist.txt
Remove all intermediate gtf files, including <file.gtf>
Create a directory /info and move Badlist.txt, and the downloaded file into it


";

	return $info;
}

__END__
