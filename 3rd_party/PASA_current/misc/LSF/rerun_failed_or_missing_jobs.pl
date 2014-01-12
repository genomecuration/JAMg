#!/usr/bin/env perl

use strict;
use FindBin;
use lib ("$FindBin::Bin/PerlLib");
use Run_Bsub;
use Cwd;
use Getopt::Long qw(:config no_ignore_case bundling);
use List::Util qw (shuffle);

my $usage = <<_EOUSAGE_;

################################################################
# Required:
#
#  -B         bsub directory
#
# Optional:
#
#  -q         grid submissions queue (defaults to 'week')
#  -M         memory to reserve (eg. 4000, which means 4 G)
#
#  --cmds_per_node    commands for each grid node to process (recommend leaving this alone).
#  --mount_test       directory for grid nodes to check for existence of (and proper mounting)
#
#  --DEBUG            
#
####################################################################

_EOUSAGE_

	;

my $help_flag;
my ($cmds_per_node, $memory, $mount_test);

my $queue = 'week';
my $bsub_dir;
my $DEBUG = 0;


&GetOptions ( 'h' => \$help_flag,
			  'B=s' => \$bsub_dir,
			  'q=s' => \$queue,
			  'M=i' => \$memory,
			  'cmds_per_node=i' => \$cmds_per_node,
			  'mount_test=s' => \$mount_test,
			  'DEBUG' => \$DEBUG,
	);


if ($help_flag || ! $bsub_dir) {
	die $usage;
}



unless (-d $bsub_dir) {
    die "Error, $bsub_dir doesn't exist.";
}

my $cmds_file = "$bsub_dir/cmds_list.txt";
my @cmds;

my $index = 0;
open (CMDS, $cmds_file) or die $!;
while (<CMDS>) {
    chomp;
    my ($index_info, $cmd) = split (/\t/, $_, 2);
    if ($index_info ne "index($index)") {
        die "Error parsing cmds_file, index out of order: index:$index, $_\n";
    }
    $cmds[$index] = $cmd;
    $index++;
}
close CMDS;

## get failed or missing cmds:
my @missing_or_failed_entries;
for (my $i = 0; $i <= $#cmds; $i++) {
    my $retval_bin = int($i/1000);
	my $entry_ret_file = "$bsub_dir/retvals/$retval_bin/entry_$i.ret";

    print STDERR "\rexamining $entry_ret_file                ";
	
    if (-s $entry_ret_file) {
        my $retval = `cat $entry_ret_file`;
        if (int($retval) != 0) {
            push (@missing_or_failed_entries, $i);
        }
    }
    else {
        # never executed
        push (@missing_or_failed_entries, $i);
    }
}

my @cmds_to_rerun;
foreach my $entry (@missing_or_failed_entries) {
    print "Rerunning entry($entry): $cmds[$entry]\n";
    push (@cmds_to_rerun, $cmds[$entry]);
}

my $logdir = cwd();

print "LOGDIR: $logdir\n";

my $num_cds = scalar @cmds_to_rerun;
my $cmds_per_node = int ($num_cds / 250);
if ($cmds_per_node < 1) {
    $cmds_per_node = 1;
}




if ($DEBUG) {
    print "DEBUG mode: not launching jobs.\n";
    print "CMDS:\n" . join ("\n", @cmds_to_rerun);
    exit(1);

}
else {

	if ($queue) {
		&Run_Bsub::set_queue($queue);
	}
 
	if ($memory) {
		&Run_Bsub::set_memory($memory);
	}

	if ($mount_test) {
		&Run_Bsub::set_mount_test($mount_test);
	}

	if ($cmds_per_node) {
        $Run_Bsub::CMDS_PER_NODE = $cmds_per_node;
	}
	
	@cmds_to_rerun = shuffle(@cmds_to_rerun);
	my $ret = &Run_Bsub::run(@cmds_to_rerun);
	
    exit($ret);
	
}

exit(0);


