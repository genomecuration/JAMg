#!/usr/bin/env perl

use strict;
use warnings;

=DESCR 

What is the expectation maximization algorithm?
Do CB, Batzoglou S.
Nat Biotechnol. 2008 Aug;26(8):897-9.
PMID: 18688245

For example in the paper, run with params:

  script_name.pl 0.6 0.5 10

=cut

my $usage = "usage: $0 Qa Qb numRounds\n\n";

my $Qa = $ARGV[0] or die $usage;
my $Qb = $ARGV[1] or die $usage;
my $num_rounds = $ARGV[2] or die $usage;

my @experiments = ( [ qw (H T T T H H T H T H) ],
					[ qw (H H H H T H H H H H) ],
					[ qw (H T H H H H H T H H) ],
					[ qw (H T H T T T H H T T) ],
					[ qw (T H H H T H H H T H)  ] );


for my $round (1..$num_rounds) {
	
	my $sum_heads_A = 0;
	my $total_A = 0;

	my $sum_heads_B = 0;
	my $total_B = 0;

	foreach my $experiment (@experiments) {
		
		my $scoreA = &score_experiment($experiment, $Qa);
		my $scoreB = &score_experiment($experiment, $Qb);
		
		my $probA = exp(-1*$scoreA);
		my $probB = exp(-1*$scoreB);

		my $ratioA = sprintf("%.2f", $probA / ($probA + $probB));
		my $ratioB = 1-$ratioA;
		
		print "$round\t$probA\t$probB\t$ratioA\t$ratioB\n";
		
		my @heads = grep { $_ eq "H" } @$experiment;
		my $num_heads = scalar(@heads);
		my $num_total = scalar (@$experiment);
		
		$sum_heads_A += $ratioA * $num_heads;
		$total_A  += $ratioA * $num_total;

		$sum_heads_B += $ratioB * $num_heads;
		$total_B += $ratioB * $num_total;
		

	}
	
	print "\n";
	

	## compute new Qa and Qb

	$Qa = $sum_heads_A / $total_A;

	$Qb = $sum_heads_B / $total_B;

    printf("New Qa: %.2f, New Qb: %.2f\n\n", $Qa, $Qb);
	
	
}



####
sub score_experiment {
	my ($experiment, $prob_heads) = @_;
	
	# returns -1 * sum_log_score

	my $prob_tails = 1 - $prob_heads;

	my @tosses = @$experiment;

	my $sum_log_score = 0;

	foreach my $toss (@tosses) {
		
		my $p = ($toss eq "H") ? $prob_heads : $prob_tails;
		
		$sum_log_score += log($p);
	}

	return(-1 * $sum_log_score);
}

