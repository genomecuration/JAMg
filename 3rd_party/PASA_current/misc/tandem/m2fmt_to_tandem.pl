#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

use lib ($ENV{EUK_MODULES});
use Overlap_piler;
use DPchain;


my $usage = "usage: $0 selfmatches.m2fmt\n\n"
	. "recommended blast settings:\n"
	. " blastn seqFile seqFile -span M=5 N=-7 -mformat=2 -kap S=500\n\n";


my $match_file = $ARGV[0] or die $usage;

my $MAX_TR_SIZE = 5000;
my $MAX_DIST_BETWEEN_SEGS = 100;
my $MIN_PER_ID = 80;

main: {
	
	my @matches = &parse_matches($match_file);

	my @piles = &pile_matches(@matches);
	
	my $pile_count = 0;
	foreach my $pile (@piles) {
		$pile_count++;
		my @chains = &chain_matches(@$pile);
		
		my $chain_count = 0;
		foreach my $chain (@chains) {
			$chain_count++;
			
			my $contig = $chain->[0]->{coordsA}->{acc};
			
			my ($tandem_flag, $coordsA_range, $coordsB_range) = &extract_coord_ranges_from_chain($chain);
			
			unless ($tandem_flag) { next; }
			
			my ($period_length, $copies, $total_length, $range) = &get_tandem_info($coordsA_range, $coordsB_range);
			my $tandem_info = "Tandem_range: $range, Period: $period_length, Copies: $copies, Len: $total_length\n";
			
			print "#Tandem\t$contig " 
				. join ("-", @$coordsA_range) . " matches " 
				. join ("-", @$coordsB_range) . "\t"
				. $tandem_info . "\n";
			
			foreach my $match (@$chain) {
				print $match->{line};
			}
			print"\n";
		}
		print "\n\n";
	}
	
	
	exit(0);


}


####
sub chain_matches {
	my (@matches) = @_;
	
	## sort according to coordinates.
	
	@matches = sort sort_matches_by_coords @matches;
	
	
	## perform 'simple' DP chaining. (no gap penalties).
	
	my $get_base_score_sref = sub { 
		my ($ele) = @_;
		my $score = ($ele->{coordsA}->{rend} - $ele->{coordsA}->{lend} + 1) * $ele->{perID}/100 
			+ ($ele->{coordsB}->{rend} - $ele->{coordsB}->{lend} + 1) * $ele->{perID}/100;
		return ($score);
	};

	my $are_chainable_sref = sub {
		my ($eleA, $eleB) = @_;
		
		# B comes before A
		if ($eleB->{coordsA}->{lend} < $eleA->{coordsA}->{lend}
			||
			$eleB->{coordsB}->{lend} < $eleA->{coordsB}->{lend}
			) 
		{ # shouldn't happen anyway, but just in case. 
			return (0); 
		}
		# B is after A, but ends before A, must be encapsulated. 
		elsif ($eleB->{coordsA}->{rend} <= $eleA->{coordsA}->{rend}
			   &&
			   $eleB->{coordsB}->{rend} <= $eleA->{coordsB}->{rend}
			   ) { 
			return (1); 
		}
		
		# B is completely after A, no overlap
		elsif ( ($eleB->{coordsA}->{lend} > $eleA->{coordsA}->{rend}  && $eleB->{coordsA}->{lend} - $eleA->{coordsA}->{rend} > $MAX_DIST_BETWEEN_SEGS)
				||
				($eleB->{coordsB}->{lend} > $eleA->{coordsB}->{rend}  && $eleB->{coordsB}->{lend} - $eleA->{coordsB}->{rend} > $MAX_DIST_BETWEEN_SEGS)
				) {
			return (0); # out of range.
		}

		## if got here, then should be within range and an acceptable linkage.
		return (1);
	};


	## find diagonals that contain all matches.
	my %match_IDs;
	foreach my $match (@matches) {
		$match_IDs{$match->{ID}} = $match;
	}
	
	my @chains;
	
	while (%match_IDs) {
		my @dp_matches = values %match_IDs;
		@dp_matches = sort sort_matches_by_coords @dp_matches;
		
		my @chain = &DPchain::find_highest_scoring_chain(\@dp_matches, $get_base_score_sref, $are_chainable_sref);

		push (@chains, [@chain]);
		foreach my $ele (@chain) {
			delete $match_IDs{ $ele->{ID} };
		}
	}

	return (@chains);
}


####
sub pile_matches {
	my @matches = @_;
	
	
	## group into piles with overlapping coordinates.
	## given the max dist between linkable segs, this will shorten the search space for DP.
	my $overlap_piler = new Overlap_piler();
	
	my %match_ID_to_match_ref; #for convenience.

	foreach my $match (@matches) {
		my $match_id = $match->{ID};
		$match_ID_to_match_ref{$match_id} = $match;
		
		my $coordsA = $match->{coordsA};
		my ($lend, $rend) = ($coordsA->{lend}, $coordsA->{rend});
		$lend -= $MAX_DIST_BETWEEN_SEGS;
		$rend += $MAX_DIST_BETWEEN_SEGS;

		$overlap_piler->add_coordSet($match_id, $lend, $rend);
	}

	my @clusters = $overlap_piler->build_clusters();
	
	my @piles; # convert clusters above into piles of matches.
	foreach my $cluster (@clusters) {
		my @pile;
		foreach my $match_ID (@$cluster) {
			my $match_ref = $match_ID_to_match_ref{$match_ID};
			push (@pile, $match_ref);
		}
		push (@piles, [@pile]);
	}


	return (@piles);
}


####
sub parse_matches {
	my ($match_file) = @_;

	my @matches;

	my $match_counter = 0;
	my %seen; # track matches, avoid duplicate entries.

	open (my $fh, $match_file) or die "Error, cannot open file $match_file";
	while (<$fh>) {
		my $line = $_;
		chomp;
		my @x = split (/\t/);
		
		my ($accA, $accB, $per_ID, $lendA, $rendA, $lendB, $rendB) = ($x[0], $x[1], $x[10], $x[17], $x[18], $x[20], $x[21]);

		if ($per_ID < $MIN_PER_ID) { next; }
		
		my ($strandA, $strandB) = ($x[16], $x[19]);

		my $orient = ($strandA * $strandB > 0) ? '+' : '-';

		($lendA, $rendA) = sort {$a<=>$b} ($lendA, $rendA);
		($lendB, $rendB) = sort {$a<=>$b} ($lendB, $rendB);

		if ($accA ne $accB) { die "Error, should only have self alignments!"; }

		if ($lendA == $lendB && $rendB == $rendB) { next; }  #self identical alignment on main diagonal, meaningless.
		
		my $coordsA = { acc => $accA,
					   lend => $lendA,
					   rend => $rendA,
				   };
		my $coordsB = { acc => $accB,
						lend => $lendB,
						rend => $rendB,
					};

		if ($lendA > $rendA) {
			next; # wu-blast reports both alignments...
			($coordsA, $coordsB) = ($coordsB, $coordsA); # keep all data left->right
		}

		my $coords_key = $accA . ":" . join ("-", $coordsA->{lend}, $coordsA->{rend}). "," . join ("-", $coordsB->{lend}, $coordsB->{rend});
		if ($seen{$coords_key}) {
			next; # wu-blast reports different alignments for the same region.
		}
		$seen{$coords_key} = 1;

		my $delta = $coordsB->{lend} - $coordsA->{lend};  # approx size of a tandem repeat if these are adjacent periods
		if ($delta > $MAX_TR_SIZE) { next; }
		

		my $match = { coordsA => $coordsA,
					  coordsB => $coordsB,
					  perID => $per_ID,
					  
					  line => $line,
					  ID => ++$match_counter,
					  
				  };
		push (@matches, $match);
	}


	return (@matches);
}


#### 
sub sort_matches_by_coords {
	
	return ( $a->{coordsA}->{lend} <=> $b->{coordsA}->{lend}                                                                            
			 ||                                                                                                                         
				 $a->{coordsB}->{lend} <=> $b->{coordsB}->{lend});
}

####
sub extract_coord_ranges_from_chain {
	my ($matches_aref) = @_;

	my @coordsA;
	my @coordsB;

	foreach my $match (@$matches_aref) {
		my $coordsA_ref = $match->{coordsA};
		push (@coordsA, $coordsA_ref->{lend}, $coordsA_ref->{rend});
		
		my $coordsB_ref = $match->{coordsB};
		push (@coordsB, $coordsB_ref->{lend}, $coordsB_ref->{rend});
	}

	@coordsA = sort {$a<=>$b} @coordsA;
	@coordsB = sort {$a<=>$b} @coordsB;
	
	my $min_A = shift @coordsA;
	my $max_A = pop @coordsA;

	my $min_B = shift @coordsB;
	my $max_B = pop @coordsB;

	my $tandem_flag = ($min_A < $min_B && $max_A < $max_B  # staggered
					   &&
					   $min_A < $max_B && $min_B < $max_A # overlap
					   ) ? 1 : 0;
	
	return ($tandem_flag, [$min_A, $max_A], [$min_B, $max_B]);
}


####
sub get_tandem_info {
	my ($coordsA_range, $coordsB_range) = @_;

	my ($lend_A, $rend_A, $lend_B, $rend_B) = (@$coordsA_range, @$coordsB_range);

	if ($lend_A > $lend_B) { confess "Error, coords out of order:  $lend_A-$rend_A, $lend_B-$rend_B "; }
	
	my @coords = sort {$a<=>$b} ($lend_A, $rend_A, $lend_B, $rend_B);
	
	my $min = shift @coords;
	my $max = pop @coords;

	my $len = $max - $min + 1;

	my $period_length = $lend_B - $lend_A + 1;
	
	my $num_copies = sprintf ("%.1f", $len / $period_length);

	my $range = "$min-$max";
	
	return ($period_length, $num_copies, $len, $range);
}
