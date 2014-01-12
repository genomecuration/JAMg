#!/usr/local/bin/perl

use Data::Dumper;
use strict;
#use warnings;

my %asmbl_id_to_core_acc;
my %core_acc_to_list;

my $MAX_CLONE_SPAN = 20000000;

open (STDERR, ">&STDOUT");

while (<STDIN>) {
    chomp;
    my @x = split (/\t/);
    my ($asmbl_id, $acc, $end5, $end3) = @x;
    
    my ($core_acc, $read_type);

    if ($acc =~ /^(\w+)(T\w{1,2})$/) {
	$core_acc = $1;
	$read_type = $2;
    } else {
	warn "$acc is in an unexpected format!\n";
	next;
    }
    
    
    my $list_href = $asmbl_id_to_core_acc{$asmbl_id};
    unless (ref $list_href) {
	$list_href = $asmbl_id_to_core_acc{$asmbl_id} = {};
    }
    $list_href->{$core_acc} = 1;
    
    my $sequenced_orient = ($end5 < $end3) ? "+" : "-";
	
    my $struct = { 
	core_acc => $core_acc,
	read_type => $read_type,
	end5 => $end5,
	end3 => $end3,
	full_acc => $acc,
	sequenced_orient => $sequenced_orient,
	asmbl_id=> $asmbl_id
	};
    
    ## keep separate record for each alignment:
    my $other_list = $core_acc_to_list{$core_acc};
    unless (ref $other_list) {
	$other_list = $core_acc_to_list{$core_acc} = [];
    }
    push (@$other_list, $struct);

}

foreach my $asmbl_id (keys %asmbl_id_to_core_acc) {

    my @full_span_entries;
    
    my $list_href = $asmbl_id_to_core_acc{$asmbl_id};
    
    foreach my $core_acc (keys %$list_href) {
	
	my $min_lend_plus_orient = undef;
	my $max_rend_minus_orient = undef;
	
	my @struct_list = @{$core_acc_to_list{$core_acc}};
	
	foreach my $struct (@struct_list) {
	
	    my ($core_acc, $read_type, $end5, $end3, $seq_orient, $mapped_asmbl_id) = ($struct->{core_acc}, 
										$struct->{read_type},
										$struct->{end5}, 
										$struct->{end3},
										$struct->{sequenced_orient},
										$struct->{asmbl_id}
										);
	    if ($mapped_asmbl_id ne $asmbl_id) {
		next;
	    }
	    
	    print "CORE_ACC: $core_acc\t$end5, $end3, $seq_orient\n";
	    my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
	    
	    if ($seq_orient eq '+') {
		if (defined ($min_lend_plus_orient) && $lend < $min_lend_plus_orient) {
		    $min_lend_plus_orient = $lend;
		} elsif (!defined($min_lend_plus_orient)) {
		    $min_lend_plus_orient = $lend;
		}
		
	    } elsif ($seq_orient eq '-') {
		
		if (defined ($max_rend_minus_orient) && $rend > $max_rend_minus_orient) {
		    $max_rend_minus_orient = $rend;
		} 
		elsif (!defined($max_rend_minus_orient)) {
		    $max_rend_minus_orient = $rend;
		}
	    } else {
		die "Error!\n";
	    }
	    
	}
	
	if ($min_lend_plus_orient && $max_rend_minus_orient) {
	    
	    if ($min_lend_plus_orient < $max_rend_minus_orient) {
				
		my $length = $max_rend_minus_orient - $min_lend_plus_orient + 1;
		
		if ($length > $MAX_CLONE_SPAN) {
		    print "ERROR_SPAN_LENGTH:\t$core_acc\t$min_lend_plus_orient\t$max_rend_minus_orient (len: $length)\n";
		} else {
		    
		    print "CLONE_SPAN:\t$core_acc\t$min_lend_plus_orient\t$max_rend_minus_orient\n";    
		    # valid span
		    my $href = { core_acc => $core_acc };
		    
		    $href->{lend} = $min_lend_plus_orient;
		    $href->{rend} = $max_rend_minus_orient;
		    $href->{asmbl_id} = $asmbl_id;
		    push (@full_span_entries, $href);

		}
		
	    } else {
		print "INVALID_VALS, sorry $core_acc, min: $min_lend_plus_orient, max: $max_rend_minus_orient\n";
	    }
	} else {
	    print "LACK_END, sorry $core_acc doesn't have both orientations accounted for. ($min_lend_plus_orient, $max_rend_minus_orient)\n";
	}
    }
    
    &get_maximal_gene_spans(@full_span_entries);

}



####
sub get_maximal_gene_spans  {
    my @full_span_entries = @_;
    
    @full_span_entries = sort {$a->{lend}<=>$b->{lend}} @full_span_entries;

    my $i=0;
    foreach my $entry (@full_span_entries) {
	$entry->{index} = $i;
	my $length = $entry->{rend} - $entry->{lend} + 1;
	$entry->{score} = $length; ## init to ele length
	$entry->{link} = -1;
	$i++;
    }

    ## chain them:
    for (my $i=1; $i <= $#full_span_entries; $i++) {
	
	my $entry_i = $full_span_entries[$i];
	my $entry_i_lend = $entry_i->{lend};
	my $entry_i_rend = $entry_i->{rend};
	my $entry_i_length = $entry_i_rend - $entry_i_lend + 1;

	    
	my $max_scoring_index = -1;
	my $max_score = 0;
	

	for (my $j = $i-1; $j >= 0; $j--) {
	    my $entry_j = $full_span_entries[$j];
	    my $score_j = $entry_j->{score};
	    
	    my $entry_j_lend = $entry_j->{lend};
	    my $entry_j_rend = $entry_j->{rend};
	    
	    ## Entry j must come before entry i
	    if ($entry_j_rend < $entry_i_lend) {
		## ok
		my $current_score = $score_j + $entry_i_length;
		if ($current_score > $max_score) {
		    $max_score = $current_score;
		    $max_scoring_index = $j;
		}
	    }
	}
	if ($max_scoring_index >= 0) {
	    $entry_i->{link} = $max_scoring_index;
	    $entry_i->{score} = $max_score;
	}
    }

    # Get highest scoring index, do traceback:
    my $max_scoring_index = -1;
    my $max_score = 0;
    for (my $i = 0; $i <= $#full_span_entries; $i++) {
	my $entry = $full_span_entries[$i];
	my $score = $entry->{score};
	if ($score > $max_score) {
	    $max_score = $score;
	    $max_scoring_index = $i;
	}
    }

    ## traceback
    my $link_index = $max_scoring_index;
    while ($link_index != -1) {
	my $entry = $full_span_entries[$link_index];
	my ($asmbl_id, $acc, $lend, $rend) = ($entry->{asmbl_id}, $entry->{core_acc}, $entry->{lend}, $entry->{rend});
	print "Traceback: $asmbl_id, $acc, $lend, $rend\n";
	$link_index = $entry->{link};
    }

    print "Done.\n\n\n\n";
    
}

    
