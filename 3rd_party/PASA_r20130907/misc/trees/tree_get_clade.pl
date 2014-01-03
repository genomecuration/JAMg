#!/usr/bin/env perl

use strict;
use warnings;

use Bio::TreeIO;
use Data::Dumper;

my $usage = "usage: $0 treeFile acc1 acc2 [acc3 ...]\n\n";

my $treeFile = $ARGV[0] or die $usage;
shift @ARGV;
my @accessions = @ARGV;

unless (scalar (@accessions) >= 2) { 
	die $usage;
}


main: {
	my $treeio = new Bio::TreeIO('-format' => 'newick',
								 '-file' => $treeFile);

	my $tree = $treeio->next_tree();

	my $root_node = $tree->get_root_node();
	
	my @nodes = $tree->get_nodes();
	
	## build a map of the nodes to accessions.
	my %acc_to_node;
	foreach my $node (@nodes) {
		if ($node->is_Leaf()) {
			my $id = $node->id();
			$acc_to_node{$id} = $node;
		}
	}

	## defensive prog.  Verify that each node exists in our parsed tree.
	my $not_found_flag = 0;
	foreach my $acc (@accessions) {
		unless (exists $acc_to_node{$acc}) {
			print STDERR "Error, cannot find node corresponding to accession [$acc]\n";
			$not_found_flag = 1;
		}
	}
	if ($not_found_flag) {
		die;
	}
		
	## start from one of the nodes, and climb to the ancestor that contains it and the others as descendants.
	my $node = $acc_to_node{$accessions[0]};
	$node = $node->ancestor();
	while (defined $node && $node ne $root_node) {
		my @descendents = $node->get_all_Descendents();
		my %accs;
		foreach my $descendent (@descendents) {
			if ($descendent->is_Leaf()) {
				my $id = $descendent->id();
				$accs{$id} = 1;
			}
		}
		
		my $found_all_flag = 1;
		foreach my $acc (@accessions) {
			if (! exists ($accs{$acc})) {
				$found_all_flag = 0; 
				last;
			}
		}

		if ($found_all_flag) {
			print join ("\n", sort keys %accs) . "\n";
			exit(0);
		}

		$node = $node->ancestor();
	}

	## if got here, then never found a clade that contained all of them.
	
	print STDERR "Error, no clade was found to contain all accessions: @accessions\n\n"
		. "Instead, all nodes in each accessions clade from the root are provided:\n\n";

	


	#################################################################
	## Retrieving all clades containing the acessions down from the root.
	##################################################################


	my %accs = ();

	foreach my $acc (@accessions) {
		
		my $node = $acc_to_node{$acc};
		
		## climb to node connected to the root:
		while ( (my $ancestor_node = $node->ancestor()) ne $root_node) {
			$node = $ancestor_node;
		}

		my @descendants = $node->get_all_Descendents();
		foreach my $othernode ($node, @descendants) {
			if ($othernode->is_Leaf()) {
				my $id = $othernode->id();
				$accs{$id} = 1;
			}
		}

	}

	print join ("\n", sort keys %accs) . "\n";
	exit(0);
		
}


