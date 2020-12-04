#!/usr/bin/env perl

use strict;
use warnings;
use URI::Escape;

my $usage = "usage: $0 btab_file dbType=[P|N]\n\n";

my $btab_file = $ARGV[0] or die $usage;
my $db_type = $ARGV[1] or die $usage;

$db_type = ($db_type eq "N") ? "cDNA_match" : "translated_nucleotide_match";

my $match_counter = 0;

main: {
	my @matches;
	my $prev_match_ID = undef;
	my $prev_contig_id = undef;
	open (my $fh, $btab_file) or die "Error, cannot open file $btab_file";
	while (<$fh>) {
		unless (/\w/) { next; }
		if (/^\#/) { next; }
		my $line = $_;
		my @x = split (/\t/);
		my $match_id = $x[13];
		my $contig_id = $x[0];

				
		if ( (! defined($prev_match_ID)) || $prev_match_ID ne $match_id) {
			&process_matches(@matches) if @matches;
			@matches = ($line);
			$prev_match_ID = $match_id;
		}
		else {
			push (@matches, $line);
		}
		
		# print contig separator
		if (defined($prev_contig_id) && $prev_contig_id ne $contig_id) {
			print "####\n";
		}
		$prev_contig_id = $contig_id;

	}
	
	&process_matches(@matches); # get last ones.
	
	
	exit(0);
}

####
sub process_matches {
	my @matches = @_;
	

	my $first_match = $matches[0];
	my @genome_coords;
	my @hit_coords;
	
	my $orient;
	
	## get the genome span:
	my $sum_len = 0;
	my $sum_len_perid = 0;
	foreach my $match (@matches) {
		my @x = split (/\t/, $match);
		my ($genome_end5, $genome_end3) = ($x[6], $x[7]);
		my ($hit_lend, $hit_rend) = ($x[8], $x[9]);
		unless ($orient) {
			if ($genome_end5 < $genome_end3) {
				$orient = '+';
			}
			elsif ($genome_end5 > $genome_end3) {
				$orient = '-';
			}
		}
		push (@genome_coords, $genome_end5, $genome_end3);
		push (@hit_coords, $hit_lend, $hit_rend);
		
		my $len = abs ($genome_end5 - $genome_end3) + 1;
		my $per_id = $x[10];
		$sum_len_perid += $len * $per_id;
		$sum_len += $len;
		
	}
	
	@genome_coords = sort {$a<=>$b} @genome_coords;
	my $lend = shift @genome_coords;
	my $rend = pop @genome_coords;
	
	@hit_coords = sort {$a<=>$b} @hit_coords;
	my $hit_lend = shift @hit_coords;
	my $hit_rend = pop @hit_coords;
	
	my $avg_per_id = sprintf("%.2f", $sum_len_perid / $sum_len);
	
	
	## print the parent feature record:
	my @parent_info = split (/\t/, $first_match);
	my $genome_contig = $parent_info[0];
	my $hit_acc = $parent_info[5];
	my $chain_id = $parent_info[13];
	my $prot_name = $parent_info[15];
	$prot_name =~ s/=/_/g;
	$prot_name = uri_escape("$hit_acc $prot_name");
	
	#$prot_name = $hit_acc; ## DEBUG
	
	
	#print join ("\t", $genome_contig, "AAT", $db_type, $lend, $rend, $avg_per_id, $orient, ".", 
	#			"ID=chain_$chain_id;Name=$prot_name; Target=$hit_acc $hit_lend $hit_rend") . "\n";
	
	foreach my $match (@matches) {
		my @x = split (/\t/, $match);
		#$match_counter++;
		my ($genome_lend, $genome_rend) = sort {$a<=>$b} ($x[6], $x[7]);
		my ($hit_lend, $hit_rend) = ($x[8], $x[9]);
		my $per_id = $x[10];
		
		
		print join ("\t", $genome_contig, "AAT", $db_type, $genome_lend, $genome_rend, $per_id, $orient, ".", 
					"ID=chain_$chain_id;Target=$hit_acc $hit_lend $hit_rend") . "\n";

	}
	print "\n"; # spacer
	
	return;
}
