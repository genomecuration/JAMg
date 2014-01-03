#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
use CdbTools;
use Nuc_translator;

my $usage = "\nusage: $0 ordered_scaffold_data chained_coords_file genome\n\n";

my $ordered_scaffold_data = $ARGV[0] or die $usage;
my $coords_file = $ARGV[1] or die $usage;
my $genome = $ARGV[2];


my $MAX_ADJACENT_DISTANCE = 10;
my $SMALL_GAP_SIZE = 6; # short N-break, but still in-frame.
my %contig_to_match_list; 

main: {
	
	my @scaffold = &parse_scaffold($ordered_scaffold_data);
	
	%contig_to_match_list = &parse_match_list($coords_file);


	my @pseudo_map;
	## init pseudomolecule:
	my $first_contig = $scaffold[0];
	my ($first_acc, $orient) = ($first_contig->{acc}, $first_contig->{orient});
	my $pseudomolecule_sequence = cdbyank_linear($first_acc, $genome);
	if ($orient eq '-') {
		$pseudomolecule_sequence = &reverse_complement($pseudomolecule_sequence);
	}
			
	push (@pseudo_map, { acc => $first_acc,
						 orient => $orient,
						 contig_lend => 1,
						 contig_rend => length($pseudomolecule_sequence),
						 genome_lend => 1,
						 genome_rend => length($pseudomolecule_sequence),
					 } );
	
	

	## append the other contigs w or w/o gaps:
	for (my $i = 1; $i <= $#scaffold; $i++) {
		
		my $prev_scaffold = $scaffold[$i-1];
		my $curr_scaffold = $scaffold[$i];
		
		my $curr_acc = $curr_scaffold->{acc};
		my $curr_orient = $curr_scaffold->{orient};

		my $sequence = cdbyank_linear($curr_acc, $genome);

		if ($curr_orient eq '-') {
			$sequence = &reverse_complement($sequence);
		}

		if (! &directly_adjacent($prev_scaffold, $curr_scaffold)) {
			## add a 60 bp N-gap
			$pseudomolecule_sequence .= 'N' x 60;
		} 
		else {
			## add 6 bp spacer
			$pseudomolecule_sequence .= 'N' x $SMALL_GAP_SIZE;
		}
		
		push (@pseudo_map, { acc => $curr_acc,
							 orient => $curr_orient,
							 contig_lend => 1,
							 contig_rend => length($sequence),
							 genome_lend => length($pseudomolecule_sequence) + 1,
							 genome_rend => length($pseudomolecule_sequence) + length($sequence), 
						 } );
		
		$pseudomolecule_sequence .= $sequence;
	}
	
	## process results:
	
	$pseudomolecule_sequence =~ s/(\S{60})/$1\n/g;
	{
		open (my $fh, ">pseudomolecule.$$.sequence") or die $!;
		print $fh ">pseudomolecule\n$pseudomolecule_sequence\n";
		close $fh;
	}
	
	{
		# write the scaffold info:
		my $prev_contig_ref;
		open (my $fh, ">pseudomolecule.$$.contig_map") or die $!;
		foreach my $contig_ref (@pseudo_map) {
			
			if ((ref $prev_contig_ref) && $contig_ref->{genome_lend} != $prev_contig_ref->{genome_rend} + $SMALL_GAP_SIZE +1) {
				print $fh "# 60 N Gap\n";
			}
			
			print $fh $contig_ref->{acc} 
			. "\t" . $contig_ref->{contig_lend}
			. "\t" . $contig_ref->{contig_rend}
			. "\t" . $contig_ref->{orient}
			. "\tPseudomolecule"
				. "\t" . $contig_ref->{genome_lend}
			. "\t" . $contig_ref->{genome_rend}
			. "\n";
			
			$prev_contig_ref = $contig_ref;

		}
		close $fh;
	}
	
	
	exit(0);
}

####
sub directly_adjacent {
	my ($prev_contig, $next_contig) = @_;

	my $prev_acc = $prev_contig->{acc};
	my $next_acc = $next_contig->{acc};


	my @prev_matches = @{$contig_to_match_list{$prev_acc}};
	my @next_matches = @{$contig_to_match_list{$next_acc}};
	
	## capture prev genome end3's:
	my @prev_end3s;
	foreach my $prev_match (@prev_matches) {
		push (@prev_end3s, $prev_match->{ref_end3});
	}
	
	my @next_end5s;
	foreach my $next_match (@next_matches) {
		push (@next_end5s, $next_match->{ref_end5});
	}


	## now, see if any pairs of coordinates are properly arranged and adjacent:
	foreach my $prev_end3 (@prev_end3s) {

		foreach my $next_end5 (@next_end5s) {

			if (abs( $next_end5 - $prev_end3) <= $MAX_ADJACENT_DISTANCE) {
				return (1); # yes, considered directly adjacent
			}
		}
	}

	return(0);  # no, not considered adjacent
		
}
	
	


####
sub parse_scaffold {
	my ($scaffold_file) = @_;

	my @scaff;
	
	open (my $fh, $scaffold_file) or die "Error, cannot open file $scaffold_file";
	while (<$fh>) {
		chomp;
		my ($acc, $len, $orient) = split (/\t/);
		push (@scaff, { acc => $acc,
						orient => $orient,
					} );
	}
	close $fh;


	return (@scaff);
}

####
sub parse_match_list {
	my ($match_file) = @_;

	my %contig_to_match_list;

	open (my $fh, $match_file) or die "Error, cannot open file $match_file";
	while (<$fh>) {
		unless (/\w/) { next; }
		if (/^\#/) { next; }
		chomp;
		my @x = split (/\t/);
		unless (scalar @x == 8) { die "Error parsing chained coordinates.  Be sure to use the chained coordinate file and not the raw coords file";}
		my ($genome_acc, $ref_end5, $ref_end3, $contig_acc, $query_end5, $query_end3, $per_id) = @x;
		

		push (@{$contig_to_match_list{$contig_acc}}, { ref_end5 => $ref_end5,
													   ref_end3 => $ref_end3,
													   query_end5 => $query_end5,
													   query_end3 => $query_end3,
													   per_id => $per_id,
													   genome_acc => $genome_acc,
													   contig_acc => $contig_acc,
												   } );
	}
	close $fh;

	return (%contig_to_match_list);
}



	
	
	
	   
