#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/PerlLib");
use Run_Bsub;
use List::Util qw (shuffle);

use Getopt::Long qw(:config no_ignore_case bundling);

my $usage = <<_EOUSAGE_;

################################################################
# Required:
#
#  -c         file containing list of commands
#  
# Optional:
#
#  -q         grid submissions queue (defaults to 'broad')
#  -M         memory to reserve (default: 4, which means 4 G)
#
#  --cmds_per_node    commands for each grid node to process (recommend leaving this alone).
#  --mount_test       directory for grid nodes to check for existence of (and proper mounting)
#
####################################################################

_EOUSAGE_

	;

my $help_flag;
my ($cmd_file, $cmds_per_node, $mount_test);

my $memory = 4;

my $queue = 'week';

if ($mount_test && ! -e $mount_test ) {
	die "Error, can't locate $mount_test ";
}

&GetOptions ( 'h' => \$help_flag,
			  'c=s' => \$cmd_file,
			  'q=s' => \$queue,
			  'M=i' => \$memory,
			  'cmds_per_node=i' => \$cmds_per_node,
			  'mount_test=s' => \$mount_test,
	
	);



unless ($cmd_file) { 
	die $usage;
}
if ($help_flag) {
	die $usage;
}


if ($cmds_per_node) {
	$Run_Bsub::CMDS_PER_NODE = $cmds_per_node;
}

main: {

	my $uname = `uname -n`;
	chomp $uname;

	print "SERVER: $uname, PID: $$\n";
	
    
    open (my $fh, $cmd_file) or die "Error, cannot open $cmd_file";
    my @cmds;

    while (<$fh>) {
        chomp;
        if (/\w/) {
            push (@cmds, $_);
        }
    }
    close $fh;

    @cmds = shuffle @cmds;  ## to even out load on grid nodes.  Some may topload their jobs!

	if ($queue) {
		&Run_Bsub::set_queue($queue);
	}
 
	if ($memory) {
		&Run_Bsub::set_memory($memory);
	}

	if ($mount_test) {
		&Run_Bsub::set_mount_test($mount_test);
	}
	
	my $ret = &Run_Bsub::run(@cmds);
	
    exit($ret);
}



    
    
