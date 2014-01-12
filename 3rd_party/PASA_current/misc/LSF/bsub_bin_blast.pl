#!/usr/bin/env perl

use strict;
use FindBin;
use lib ("FindBin::Bin/../PerlLib", "$FindBin::Bin/PerlLib");
use Fasta_reader;
use Getopt::Std;
use strict;
use Carp;
use Cwd;
use Bsub;
use List::Util qw (shuffle);

our ($opt_d, $opt_q, $opt_s, $opt_Q, $opt_b, $opt_p, $opt_O, $opt_h, $opt_c, $opt_B, $opt_X, $opt_M);

&getopts ('dq:s:bp:O:hbc:B:Q:XM:');

my $usage =  <<_EOH_;

############################# Options ###############################
#
# -q query multiFastaFile (full or relative path)
# -s search multiFastaFile (full or relative path)
# -p program (e.g. blastn blastx tblastn)
# -O blast options  ie. "-filter xnu+seg"
# -b btab only
# -c cmds per node 
# 
# -B bin size  (input seqs per directory)  (default 5000)
# -Q bsub queue
# -M memory   (4000 = 4G is the default setting).  use nubmers only... don't say 4G!!!
# -X commands only!  don't launch.
#
###################### Process Args and Options #####################

_EOH_

    
    ;


if ($opt_h) {
    die $usage;
}

my $CMDS_ONLY = $opt_X;

my $bin_size = $opt_B || 5000;

my $queue = $opt_Q || "week";

our $DEBUG = $opt_d;

my $CMDS_PER_NODE = $opt_c;

unless ($opt_q && $opt_s && $opt_p) {
    die $usage;
}

my $queryFile = $opt_q;
unless ($queryFile =~ /^\//) {
    $queryFile = cwd() . "/$queryFile";
}
my $searchDB = $opt_s;
unless ($searchDB =~ /^\//) {
    $searchDB = cwd() . "/$searchDB";
}

unless (-s $searchDB || "-s $searchDB.pal") {
    die "Error, can't find $searchDB\n";
}

my $program = $opt_p;

my $progToken = $program;


my $progOptions = $opt_O;

my $memory = $opt_M || 4000;

## Create files to search

my $fastaReader = new Fasta_reader($queryFile);

my @searchFileList;

my $count = 0;
my $current_bin = 1;
my $bindir = "grp_" . sprintf ("%04d", $current_bin);
mkdir ($bindir) or die "Error, cannot mkdir $bindir";

while (my $fastaEntry = $fastaReader->next() ) {
    my $acc = $fastaEntry->{accession};
    $acc =~ s/\W/_/g;
    
    unless ($acc) {
        print "Error, fasta entry is missing an accession.\n";
        next;
    }
    
	$count++;
	
    my $filename = "$bindir/$acc";
                
    push (@searchFileList, $filename);
	
	my $sequence = $fastaEntry->{sequence};
    $sequence =~ s/(\w{60})/$1\n/g; #make fasta format.
    my $header = $fastaEntry->{header};
    open (TMP, ">$filename") or die "Can't create file ($filename)\n";
    print TMP ">$header\n$sequence\n";
    close TMP;
    chmod (0666, $filename);
    
	
	if ($count % $bin_size == 0) {
		# make a new bin:
		$current_bin++;
		$bindir = "grp_" . sprintf ("%04d", $current_bin);
		mkdir ($bindir) or die "Error, cannot mkdir $bindir";
	}
}

print "Sequences to search: @searchFileList\n";
my $numFiles = @searchFileList;
print "There are $numFiles sequences to be searched.\n";

my $curr_dir = cwd;

if  ($numFiles) {
    
    my @cmds;
    ## formulate blast commands:
    foreach my $searchFile (@searchFileList) {
        $searchFile = "$curr_dir/$searchFile";
        
		my $cmd = "blastall -p $program -d $searchDB -i $searchFile $progOptions > $searchFile.$progToken.result ";
		unless ($CMDS_ONLY) {
			$cmd .= "2>$searchFile.$progToken.stderr";
		}
        push (@cmds, $cmd);
    }
    
	
    @cmds = shuffle(@cmds);

    my $cmds_per_node;
    if ($CMDS_PER_NODE) {
        $cmds_per_node = $CMDS_PER_NODE;
    }
    else {
        $cmds_per_node = int ( scalar(@cmds) / 400); # split job across complete set of nodes available.
        if ($cmds_per_node < 1) {
            # use 10 as default.
            $cmds_per_node = 1;
        }
    }
    

	open (my $fh, ">cmds.list") or die $!;
	foreach my $cmd (@cmds) {
		print $fh "$cmd\n";
	}
	close $fh;
	
	unless ($CMDS_ONLY) {
		my $bsubber = new Bsub({cmds=>\@cmds,
								log_dir => $curr_dir,
								cmds_per_node => $cmds_per_node,
								queue => $queue,
								memory => $memory,
							}
			);
		
		$bsubber->bsub_jobs();
		
	}

} else {
    print STDERR "Sorry, no searches to perform.  Results already exist here\n";
}
## Cleanup

exit(0);
