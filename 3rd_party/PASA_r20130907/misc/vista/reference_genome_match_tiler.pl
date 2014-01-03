#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long qw(:config no_ignore_case bundling);
use Data::Dumper;

my ($mummer_coords_file, $mummer_delta_file, $help);

my $min_match_length = 0;
my $max_overlap_bases = 100; #default
my $max_dist_between_segs = 1000;
my $method = "";

my $DEBUG = 0;

&GetOptions ('mummer_coords_file=s' => \$mummer_coords_file,
             'mummer_delta_file=s' => \$mummer_delta_file,
             'min_match_length=i' => \$min_match_length,
             'max_overlap_bases=i' => \$max_overlap_bases,
             'max_dist_between_segs=i' => \$max_dist_between_segs,
             'method=s' => \$method,
             'h' => \$help,
			 'debug' => \$DEBUG,
			 );

my $usage = <<_EOUSAGE_;

###########################################################################################

Instructions:

$0 maps draft genome assemblies of related organisms to a single reference molecule.  
Requirements are:
-a single reference molecule in Fasta format (ie.  reference.fasta)
-a single fasta file containing all other assembled contigs (ie. assemblies.fasta)

First, use nucmer to align all assemblies to the reference:

%  nucmer reference.fasta assemblies.fasta

This will generate a nucmer alignment output file called 'out.delta'.  Obtain a tab-delimited
summary of the match coordinates by running the 'show-coords' utility on the out.delta file like so:

%  show-coords -H -T out.delta > out.coords

Now run this script to extract the best ordering of assemblies to the reference:

%  $0 --mummer_coords_file out.coords --method nucmer --mummer_delta_file out.delta > out.delta.refined

Running the above will generate the output file 'out.delta.refined', which will contain those relevant 
alignments between mapped assemblies and the reference.  It will also create a file called 
'ordered_query_list.txt', which creates the ordering and orientation of the assemblies as they match the reference.

Run mummerplot to verify the mappings and to visualize the relevant alignments like so:

%  mummerplot -Q ordered_query_list.txt out.delta.refined

Run show-coords on the out.delta.refined to obtain a more human readable alignment summary.



############################################################################################
#
# required:
#
# --mummer_coords_file     nucmer coords-formatted file 'use show-coords -T out.delta'
#
# --method                 nucmer | promer
#
# optional:
#  
# --min_match_length       default zero
# --mummer_delta_file      if specified, only the relevant subset will be reported.
# --max_overlap_bases      number of overlap allowed along the reference sequence and matching contigs (default 100)
# --max_dist_between_segs  maximum length between two chained segments along query or reference   (default 1000);
#
# --debug
#############################################################################################


_EOUSAGE_

    ;

if  ((! $mummer_coords_file) || $help || ! $method) {
    die $usage;
}

if ($mummer_delta_file && ! -s $mummer_delta_file) {
    die "Error, cannot locate $mummer_delta_file";
}

unless ($method =~ /^(nucmer|promer)$/i) {
    die "Error, don't understand method \"$method\"";
}

my %QUERY_LENGTHS;


open (my $log_fh, ">log.txt") or die $!;

main: {
    
    my @matches = &parse_match_file($mummer_coords_file);
    
    unless (@matches) {
        ## nothing to do.
        print STDERR "no matches reported in $mummer_coords_file\n";
        exit(0);
    }

    &build_DP_trellis(@matches);

    my $highest_scoring_node = &find_highest_scoring_node(@matches);

    ## trace back from highest scoring node. (node and match are synonymous)
    my $trace_node = $highest_scoring_node;
    
    my @traced_nodes;
    while ($trace_node) {
        push (@traced_nodes, $trace_node);
        $trace_node = $trace_node->{prev_best_node};
    }

	if ($DEBUG) {
		&write_match_list("debug.$$.afterDP", \@traced_nodes);
	}
	
    
    # report traced nodes:
    if ($mummer_delta_file) {
        
        ## extract matches from the mummer delta file.
        ## also, order matching contigs based on the midpt of the longest chain of matches.
        &retrieve_subset_delta_file(@traced_nodes);
    }
    else {
        ## chain together the query seqs based on same orientation and max gap
        &chain_matches_into_contig_regions(@traced_nodes);
    }
    
    exit(0);
}

####
sub chain_matches_into_contig_regions {
    my (@nodes) = @_;

    ## sort by reference position:
    @nodes = sort {$a->{ref_lend}<=>$b->{ref_lend}} @nodes;

    
    my $prev_node = shift @nodes;
    my @chains = [$prev_node];
    
    while (@nodes) {
        
        my $curr_chain = $chains[$#chains];
        
        my ($prev_query_acc, $prev_query_lend, $prev_query_rend, $prev_query_orient, $prev_ref_orient) = ($prev_node->{query_acc},
																										  $prev_node->{query_lend},
																										  $prev_node->{query_rend},
																										  $prev_node->{query_orient},
																										  $prev_node->{ref_orient});
        
        my $next_node = shift @nodes;
        my ($next_query_acc, $next_query_lend, $next_query_rend, $next_query_orient, $next_ref_orient) = ($next_node->{query_acc},
																						 $next_node->{query_lend},
                                                                                         $next_node->{query_rend},
                                                                                         $next_node->{query_orient},
                                                                                         $next_node->{ref_orient});
        
        if (  ($prev_query_acc eq $next_query_acc) 
			  &&
			  ($prev_query_orient eq $next_query_orient && $prev_ref_orient eq $next_ref_orient)
              && 
              (
               (abs ($prev_query_rend - $next_query_lend) < $max_dist_between_segs )
               || 
               (abs ($prev_query_lend - $next_query_rend) < $max_dist_between_segs) 
               ) 
              ) 
        {
                                                                                             
                                                                                             ## consistent segment, add to prev chain:
            push (@$curr_chain, $next_node);
        }
        else {
            ## create a new chain:
            push (@chains, [$next_node]);
        }
    
        $prev_node = $next_node;

    }

    foreach my $chain (@chains) {
        
        my $chain_text = "";
        
        my $ref_orientation;
        my $query_orientation;
        my @ref_coords;
        my @query_coords;

        my $ref_accession = $chain->[0]->{ref_acc};
        my $query_accession = $chain->[0]->{query_acc};
        
        foreach my $match (@$chain) {
            $chain_text .= &get_match_text($match);
            
            my ($ref_lend, $ref_rend, $ref_orient,
                $query_lend, $query_rend, $query_orient) = ($match->{ref_lend}, $match->{ref_rend}, $match->{ref_orient},
                                                            $match->{query_lend}, $match->{query_rend}, $match->{query_orient});
            
            $ref_orientation = $ref_orient;
            $query_orientation = $query_orient;
            push (@ref_coords, $ref_lend, $ref_rend);
            push (@query_coords, $query_lend, $query_rend);
        }

        @ref_coords = sort {$a<=>$b} @ref_coords;
        my $ref_lend = shift @ref_coords;
        my $ref_rend = pop @ref_coords;
        @query_coords = sort {$a<=>$b} @query_coords;
        my $query_lend = shift @query_coords;
        my $query_rend = pop @query_coords;

        print "# Alignment Chain:\t$ref_accession\t$ref_lend-$ref_rend\[$ref_orientation]\t$query_accession\t$query_lend-$query_rend\[$query_orientation]\n";
        print $chain_text . "\n";
        
    }
    
    return;
}


####
sub parse_match_file {
    my $match_file = shift;

    my @matches;

    open (my $fh, $match_file) or confess "Error, cannot open file $match_file\n";
    while (<$fh>) {
        unless (/^\d+/) { next; }
        chomp;
        my @x = split (/\t/);
        
        my ($ref_end5, $ref_end3, $query_end5, $query_end3, $ref_match_len, $query_match_len, $per_id, $ref_acc, $query_acc) = @x;

        if ($method =~ /promer/) {
            # output is a bit different...
            # set the ref and query accessions as the last two fields.
            
            $query_acc = pop @x;
            $ref_acc = pop @x;
        }
                

        my $ref_orient = ($ref_end5 < $ref_end3) ? '+' : '-';
        my $query_orient = ($query_end5 < $query_end3) ? '+' : '-';
        
        my ($ref_lend, $ref_rend) = sort {$a<=>$b} ($ref_end5, $ref_end3);
        my ($query_lend, $query_rend) = sort {$a<=>$b} ($query_end5, $query_end3);

        if ( (! exists $QUERY_LENGTHS{$query_acc}) || ($query_rend > $QUERY_LENGTHS{$query_acc}) ) {
            $QUERY_LENGTHS{$query_acc} = $query_rend;
        }
        
        unless ($ref_match_len >= $min_match_length) { next; }

        my $match_struct = {  
            ref_lend => $ref_lend,
            ref_rend => $ref_rend,
            ref_orient => $ref_orient,
            ref_acc => $ref_acc,
            
            query_lend => $query_lend,
            query_rend => $query_rend,
            query_orient => $query_orient,
            query_acc => $query_acc,
            
            per_id => $per_id,

            match_length => $ref_match_len,
            
            ## attributes for DP trace
            sum_novel_query_bp => $query_match_len,
            prev_best_node => undef, # link to previous best node in graph for DP tracing.
            query_matching_regions_thus_far => { $query_acc => [ 
                                                                 [$query_lend, $query_rend] 
                                                                 ] 
                                                             }, # list of all nonoverlapping regions of query matched up to this point in the graph, indexed by query accession.
                                                                 
                                                                 
           };
        
        push (@matches, $match_struct);
    }
    
    return (@matches);
    
}


####
sub build_DP_trellis {
    my @matches = @_;

    # sort matches by reference coordinates
    @matches = sort {$a->{ref_lend} <=> $b->{ref_lend}} @matches;
    
	
	if ($DEBUG) {
		&write_match_list("debug.$$.beforeDP", \@matches);
	}
	
	my $trellis_fh;
	if ($DEBUG) {
		open ($trellis_fh, ">debug.$$.trellis_build") or die "Error, cannot open trellis debug file";
	}
	
	# compare each node to the previous node
    for (my $i = 1; $i <= $#matches; $i++) {
        
        my $node_i = $matches[$i];
        my $node_i_novel_bp = $node_i->{sum_novel_query_bp};
        my $node_i_query_matching_regions_thus_far = $node_i->{query_matching_regions_thus_far};
		my $node_i_query_acc = $node_i->{query_acc};
		my ($node_i_ref_lend, $node_i_ref_rend) = ($node_i->{ref_lend}, $node_i->{ref_rend});
        
        # node j comes before node i in the graph.
        for (my $j = $i - 1; $j >= 0; $j--) {
                        
            my $node_j = $matches[$j];
            my $node_j_novel_bp = $node_j->{sum_novel_query_bp};
            my $node_j_query_matching_regions_thus_far = $node_j->{query_matching_regions_thus_far};
			my $node_j_query_acc = $node_j->{query_acc};
			my ($node_j_ref_lend, $node_j_ref_rend) = ($node_j->{ref_lend}, $node_j->{ref_rend});

			print $trellis_fh "$j $node_j_query_acc [$node_j_ref_lend-$node_j_ref_rend] vs. $i $node_i_query_acc [$node_i_ref_lend-$node_i_ref_rend] " if $DEBUG;

            unless ( (my ($ref_overlap) = &ref_coords_overlap($node_i, $node_j)) <= $max_overlap_bases) { 
                # comparison not allowed.  
				print $trellis_fh "not_allowed.  $ref_overlap along the refernece coordinates.\n" if $DEBUG;
                next; 
            }
            
            my %join_i_j_novel_coords = &join_i_j_novel($node_i_query_matching_regions_thus_far,
                                                        $node_j_query_matching_regions_thus_far);
            
            my $sum_bases = &sum_regions(%join_i_j_novel_coords);
	
	    
	    my $novel_bases_contributed_by_i = $sum_bases - $node_j_novel_bp;
	    
		
            if ($sum_bases > $node_i_novel_bp && $novel_bases_contributed_by_i > 0) {

                ## got a new higher score
				print $trellis_fh " adds " . ($sum_bases - $node_j_novel_bp) . " novel query seq to alignment.\n" if $DEBUG;
				$node_i_novel_bp = $node_i->{sum_novel_query_bp} = $sum_bases;
                $node_i->{prev_best_node} = $node_j;
                $node_i->{query_matching_regions_thus_far} = {%join_i_j_novel_coords};

            }
			else {
				print $trellis_fh " no novel addition.\n" if $DEBUG;
			}
        }
    }

	if ($DEBUG) {
		print $trellis_fh "Done.\n";
		close $trellis_fh;
	}
	
	return;
}

####
sub join_i_j_novel {
    my ($regions_B_href, $regions_A_href) = @_;

    # A may have many query_accs, B should have only one providing it's single component alignment
    
    my @query_accs = keys %$regions_B_href;
    
    if (scalar @query_accs != 1) {
        confess "Error, regions_B_href should have only one query alignment.";
    }
    
    my $query_acc = $query_accs[0];
    my $matches_B_aref = $regions_B_href->{$query_acc};

    my $matches_A_aref = $regions_A_href->{$query_acc};
    unless ($matches_A_aref) {
        # region_B contributes completely novel sequence here.
        my %combined = (%$regions_A_href, %$regions_B_href);
        return (%combined);
    }
    
    ## join them into non-redundant overlapping regions
    my @combined_coords = &Overlap_piler::simple_coordsets_collapser(@$matches_A_aref, @$matches_B_aref);
    
    my %combined = %$regions_A_href;
    $combined{$query_acc} = [@combined_coords];
    

    return (%combined);
}

####
sub sum_regions {
    my %regions = @_;
    
    my $sum_length = 0;
    
    foreach my $query_acc (keys %regions) {
        foreach my $match (@{$regions{$query_acc}}) {
            my ($lend, $rend) = @$match;
            my $len = abs ($rend - $lend) + 1;
            $sum_length += $len;
        }
    }
    return ($sum_length);
    
}

####
sub find_highest_scoring_node {
    my @matches = @_;

    my $best_scoring_node = shift @matches;
    foreach my $match (@matches) {
        if ($match->{sum_novel_query_bp} > $best_scoring_node->{sum_novel_query_bp}) {
            $best_scoring_node = $match;
        }
    }

    return ($best_scoring_node);
}

####
sub get_match_text {
    my ($match) = @_;

    my $ref_lend = $match->{ref_lend};
    my $ref_rend = $match->{ref_rend};
    my $ref_orient = $match->{ref_orient};
    my $ref_acc = $match->{ref_acc};

    my $query_lend = $match->{query_lend};
    my $query_rend = $match->{query_rend};
    my $query_orient = $match->{query_orient};
    my $query_acc = $match->{query_acc};

    my $per_id = $match->{per_id};
    my $match_length = $match->{match_length};

    my ($ref_end5, $ref_end3) = ($ref_orient eq '+') ? ($ref_lend, $ref_rend) : ($ref_rend, $ref_lend);
    my ($query_end5, $query_end3) = ($query_orient eq '+') ? ($query_lend, $query_rend) : ($query_rend, $query_lend);

    return ("$ref_acc\t$ref_end5\t$ref_end3\t$query_acc\t$query_end5\t$query_end3\t$per_id\t$match_length\n");
    
}

####
sub write_match_list {
	my ($file, $match_list_aref) = @_;

	open (my $fh, ">$file") or die "Error, cannot open file $file";
	foreach my $match (@$match_list_aref) {
		print $fh &get_match_text($match);
	}
	close $fh;
	
	return;
}


####
sub retrieve_subset_delta_file {
    my @traced_nodes = @_;

    my %alignment_tokens; #just those alignments to print
    foreach my $traced_node (@traced_nodes) {
        my $match_text = &get_match_text($traced_node);
        chomp $match_text;
        my ($ref_acc, $ref_end5, $ref_end3, $query_acc, $query_end5, $query_end3, $per_id, $match_len) = split(/\t/, $match_text);

        my $token = join ("__", $ref_acc, $ref_end5, $ref_end3, $query_acc, $query_end5, $query_end3);
        $alignment_tokens{$token} = 1;
    }

    ## parse the delta file:
    my $curr_header;
    my $curr_text = "";
    my ($ref_acc, $query_acc);

    open (my $fh, $mummer_delta_file) or die "error, cannot open $mummer_delta_file";
    # report the mummer header info (top 2 lines)
    my $line = <$fh>;
    print $line;
    $line = <$fh>;
    print $line;
    
    my @entries; ## store delta alignment text so we can sort it according to reference coordinates:
           
    my $ref_end5_stored;
    my $query_acc_stored;

    while (<$fh>) {
        if (/^>/) {
            if ($curr_text) {
                # print $curr_text;
                push (@entries, { text => $curr_text,
                                  coord => $ref_end5_stored,
                                  query_acc => $query_acc_stored } );
            }
            $curr_text = "";
            $curr_header = $_;
            my @rest;
            ($ref_acc, $query_acc, @rest) = split (/\s+/);
            $ref_acc =~ s/>//;
        }
        else {
            my @x = split (/\s+/);
            if (scalar (@x) > 1) {
                my ($ref_end5, $ref_end3, $query_end5, $query_end3, @rest) = @x;
                
                my $token = join ("__", $ref_acc, $ref_end5, $ref_end3, $query_acc, $query_end5, $query_end3);
                if ($alignment_tokens{$token}) {
                    # got one. record it.
                    unless ($curr_text) {
                        # add the header
                        $curr_text = $curr_header;
                        $ref_end5_stored = $ref_end5;
                        $query_acc_stored = $query_acc;
                    }
                    $curr_text .= $_;
                    while (my $line = <$fh>) {
                        $curr_text .= $line;
                        my $value = $line;
                        chomp $value;
                        if ($value == 0) {
                            last;
                        }
                    }
                }
            }
        }
    }
    if ($curr_text) {
        push (@entries, { text => $curr_text,
                          coord => $ref_end5_stored,
                          query_acc => $query_acc_stored } );
    }
    
    
    # print ordered list of query accessions for prettier nucmer plot.
    open ($fh, ">ordered_query_list.txt") or die "Error, cannot open ordered_query_list.txt";
   
    ## sort according to coordinate and print delta text:
    
    # sort based on an average coordinate weighting scheme.
    
    foreach my $entry (@entries) {
        my $text = $entry->{text};
        
        my $sum_lengths = 0;
        my $sum_match_pos = 0;
        my @lines = split (/\n/, $text);
        shift @lines; # rid the header;
        
        my %orientation_voter = ( '+' => 0, '-' => 0);
        
        my @component_matches;

        foreach my $line (@lines) {
            chomp $line;
            my @x = split (/\s+/, $line);
            if (scalar @x > 1) {
                my ($ref_lend, $ref_rend, $query_end5, $query_end3) = ($x[0], $x[1], $x[2], $x[3]);
                my $ref_length = abs ($ref_rend - $ref_lend) + 1;
                my $midpt = ($ref_lend + $ref_rend) / 2;
                $sum_lengths += $ref_length;
                $sum_match_pos += $ref_length * $midpt;
                
                my $match_orient = ($query_end5 < $query_end3) ? '+' : '-';
                $orientation_voter{$match_orient} += $ref_length;
                
                my ($query_lend, $query_rend) = sort {$a<=>$b} ($query_end5, $query_end3);
                
                push (@component_matches, {  ref_lend => $ref_lend,
                                             ref_rend => $ref_rend,
                                             query_lend => $query_lend,
                                             query_rend => $query_rend, 
                                             orient => $match_orient,
                                         } );
                
            }
        }
        
        if ($sum_lengths == 0) { 
            die "$text\nError, no sum lengths.\n";
        }
        
        # use weighted average for coordinate positioning:
        my $coord = $sum_match_pos / $sum_lengths;
        $entry->{coord} = $coord;

        ## set the orientation based on majority vote:
        my $contig_orient = ($orientation_voter{'+'} > $orientation_voter{'-'}) ? '+' : '-';
        
        $entry->{orient} = $contig_orient;
        
        $entry->{component_matches} = [@component_matches];
        
        my $CHAIN_CONTIG_SEGS = 1;
        if ($CHAIN_CONTIG_SEGS) {
            &assign_entry_coordinate_by_longest_chain_midpt($entry);
        }
        

    }

	if (1) {
		@entries = &DP_minimal_ref_overlap(@entries);
	}
    
	open (my $extended_fh, ">ordered_query_list.txt.w_coords") or die $!;
    ## sort by coordinate
    # @entries = sort {$a->{coord}<=>$b->{coord}} @entries;
    foreach my $entry (@entries) {
        
        print $entry->{text};
        
        my $query_acc = $entry->{query_acc};
        my $query_orient = $entry->{orient};
		
		print $fh "$query_acc\t" . $QUERY_LENGTHS{$query_acc} . "\t$query_orient\n";
        print $extended_fh "$query_acc\t" . $QUERY_LENGTHS{$query_acc} . "\t$query_orient\t" . $entry->{left_ref_span} . "-" . $entry->{right_ref_span} . "\n";
    
    }
    close $fh;
    close $extended_fh;

    return;
}

####
sub ref_coords_overlap {
    my ($match_A, $match_B) = @_;

    my ($match_A_lend, $match_A_rend) = ($match_A->{ref_lend}, $match_A->{ref_rend});
    my ($match_B_lend, $match_B_rend) = ($match_B->{ref_lend}, $match_B->{ref_rend});

    # if they do not overlap, then overlap = 0
    
    unless ($match_A_lend <= $match_B_rend && $match_A_rend >= $match_B_lend) {
        return (0);
    }

    return (&nucs_in_common($match_A_lend, $match_A_rend, $match_B_lend, $match_B_rend));
}

    
####
sub nucs_in_common {
    my ($e5, $e3, $g5, $g3) = @_;
    ($e5, $e3) = sort {$a<=>$b} ($e5, $e3);
    ($g5, $g3) = sort {$a<=>$b} ($g5, $g3);
    my $length = abs ($e3 - $e5) + 1;
    my $diff1 = ($e3 - $g3);
    $diff1 = ($diff1 > 0) ? $diff1 : 0;
    my $diff2 = ($g5 - $e5);
    $diff2 = ($diff2 > 0) ? $diff2 : 0;
    my $overlap_length = $length - $diff1 - $diff2;
    return ($overlap_length);
}


####
sub assign_entry_coordinate_by_longest_chain_midpt {
    my ($entry) = @_;

    my $query_acc = $entry->{query_acc};
    
    my @component_matches = @{$entry->{component_matches}};
    
    ## sort by reference coordinate position.
    @component_matches = sort {$a->{ref_lend}<=>$b->{ref_lend}} @component_matches;
    
    # init score and prev for DP
    foreach my $match (@component_matches) {
        my $length = $match->{query_rend} - $match->{query_lend} + 1;
        $match->{score} = $length;
        $match->{sum_score} = $length;
        $match->{prev} = undef;
    }
     
	my $verbose = 0;
	
    for (my $i = 1; $i <= $#component_matches; $i++) {
        
        my $match_i = $component_matches[$i];
        my $orient_i = $match_i->{orient};
        my ($i_lend, $i_rend) = ($match_i->{query_lend}, $match_i->{query_rend});
        my ($i_ref_lend, $i_ref_rend) = ($match_i->{ref_lend}, $match_i->{ref_rend});
        
        # score values.
        my $i_score = $match_i->{score};
        my $i_sum_score = $match_i->{sum_score};
        
        
        for (my $j = $i - 1; $j >= 0; $j--) {
            my $match_j = $component_matches[$j];
            my $orient_j = $match_j->{orient};
            my ($j_lend, $j_rend) = ($match_j->{query_lend}, $match_j->{query_rend});
            my ($j_ref_lend, $j_ref_rend) = ($match_j->{ref_lend}, $match_j->{ref_rend});
            
			print "comparing ($i) $i_lend-$i_rend $i_ref_lend-$i_ref_rend vs. ($j) $j_lend-$j_rend $j_ref_lend-$j_ref_rend\n" if $verbose;

            unless ($orient_i eq $orient_j) { 
				print "\topposite orients.\n" if $verbose;
				next; 
			} 
            
            if ($orient_j eq '-') {
                ## swap i and j query coordinates  
                ($i_lend, $i_rend, $j_lend, $j_rend) = ($j_lend, $j_rend, $i_lend, $i_rend);
            }
            
            unless ($j_rend < $i_lend + $max_overlap_bases) { 
				print "\ttoo much query overlap\n" if $verbose;
				next; 
			} # don't allow more than that overlap along query
			
            unless ($j_ref_rend < $i_ref_lend + $max_overlap_bases) { 
				print "\ttoo much reference overlap\n" if $verbose;
				next; 
			} # ditto along the reference.

            my $dist_between_segs = $i_lend - $j_rend;
            if ($dist_between_segs > $max_dist_between_segs) { 
				print "\ttoo much query gap ($dist_between_segs vs. $max_dist_between_segs)\n" if $verbose;
				next; 
			} # intermatch length too high!
            
            my $ref_dist_between_segs = $i_ref_lend - $j_ref_rend;  # i after j
            if ($ref_dist_between_segs > $max_dist_between_segs) { 
				print "\ttoo much reference gap ($ref_dist_between_segs vs $max_dist_between_segs)\n" if $verbose;
				next; 
			}
            
            my $j_sum_score = $match_j->{sum_score};

            if ($i_score + $j_sum_score > $i_sum_score) {
                ## got new best score
                $i_sum_score = $match_i->{sum_score} = $i_score + $j_sum_score;
                $match_i->{prev} = $match_j;
				print "\t** linking $i to $j\n" if $verbose;
			}
            
        }
    }
    
    ## get longest chain of matches;
    my $highest_scoring_match = shift @component_matches;
    foreach my $other_match (@component_matches) {
        if ($other_match->{sum_score} > $highest_scoring_match->{sum_score}) {
            $highest_scoring_match = $other_match;
        }
    }

    my @best_chain;
    while ($highest_scoring_match) {
        push (@best_chain, $highest_scoring_match);
        $highest_scoring_match = $highest_scoring_match->{prev};
    }

    ## determine midpt of chain:
    my @coords;
    
	my $sum_segment_lengths = 0;
    ### Note: Only the single best alignment chain is logged per contig!!!  
    print $log_fh "// Alignment Chain\n";
    foreach my $match (@best_chain) {
        
        print $log_fh "$query_acc\t" . $match->{ref_lend} . "-" . $match->{ref_rend} 
        . "\t" . $match->{query_lend} . "-" . $match->{query_rend} 
        . "\t" . $match->{orient} . "\n";
        
        push (@coords, $match->{ref_lend}, $match->{ref_rend});
    
		my $match_len = $match->{ref_rend} - $match->{ref_lend} + 1;
		$sum_segment_lengths += $match_len;
	}
    @coords = sort {$a<=>$b} @coords;

    my $lend = shift @coords;
    my $rend = pop @coords;

	## store midpt of longest chain and the span of the longest chain:
    my $midpt = ($rend + $lend)/2;

    $entry->{coord} = $midpt;
    
	$entry->{left_ref_span} = $lend;
	$entry->{right_ref_span} = $rend;
	$entry->{_DP_ref_score} = $sum_segment_lengths;
	
    print $log_fh "$query_acc\tMID: $midpt\n\n";
    
    # print STDERR Dumper (\@best_chain);
    
    return;
}


####
sub DP_minimal_ref_overlap {
	my @entries = @_;

	
	@entries = sort {$a->{left_ref_span}<=>$b->{left_ref_span}} @entries;

	## init:
	foreach my $entry (@entries) {
		$entry->{_DP_sum_score} = $entry->{_DP_ref_score};
	}
	
	for (my $i = 1; $i <= $#entries; $i++) {
		
		my $entry_i = $entries[$i];
		my $i_score = $entry_i->{_DP_ref_score};
		my ($i_ref_lend, $i_ref_rend) = ($entry_i->{left_ref_span}, $entry_i->{right_ref_span});
		my $i_sum_score = $entry_i->{_DP_sum_score};

		for (my $j = $i-1; $j >= 0; $j--) {

			my $entry_j = $entries[$j];
			my $j_score = $entry_j->{_DP_ref_score};
			my ($j_ref_lend, $j_ref_rend) = ($entry_j->{left_ref_span}, $entry_j->{right_ref_span});
			my $j_sum_score = $entry_j->{_DP_sum_score};
			
			if ($i_ref_lend > $j_ref_rend - $max_overlap_bases 
				&&
				$i_ref_lend > $j_ref_lend # no ties
				&& ($i_score + $j_sum_score > $i_sum_score)
				) {
				$i_sum_score = $entry_i->{_DP_sum_score} = ($i_score + $j_sum_score);
				$entry_i->{_DP_prev} = $entry_j;
			}
		}
	}


	my $highest_score = 0;
	my $node = undef;
	foreach my $entry (@entries) {
		if ($entry->{_DP_sum_score} > $highest_score) {
			$highest_score = $entry->{_DP_sum_score};
			$node = $entry;
		}
	}

	my @ret_nodes;
	while (defined $node) {
		push (@ret_nodes, $node);
		$node = $node->{_DP_prev};
	}
	
	return (reverse @ret_nodes);
}

our $SEE = 0;

package Overlap_piler;


use strict;
use warnings;
use Carp;

sub new {
    my $packagename = shift;
    my $self = {
        node_list => []
        };
    bless ($self, $packagename);
    return ($self);
}


#### Static method!!!
sub simple_coordsets_collapser {
    my @coordsets = @_; ## list of coordinates [a, b], [c, d], ...
    
    my $counter = 0;
    
    my $piler = new Overlap_piler();

    my %coords_mapping;
    foreach my $coordset (@coordsets) {
        $counter++;
        $coords_mapping{$counter} = [@$coordset];
        
        my ($lend, $rend) = @$coordset;
        if ($lend !~ /\d/ || $rend !~ /\d/) {
            confess "Error, coordinates [ $lend, $rend ] include a non-number";
        }

        $piler->add_coordSet($counter, @$coordset);
    }

    my @clusters = $piler->build_clusters();
    
    my @coord_spans;
    foreach my $cluster (@clusters) {
        my @eles = @$cluster;
        my @coords;
        foreach my $ele (@eles) {
            push (@coords, @{$coords_mapping{$ele}});
        }

        @coords = sort {$a<=>$b} @coords;
        
        my $min_coord = shift @coords;
        my $max_coord = pop @coords;

        push (@coord_spans, [$min_coord, $max_coord]);
    }

    return (@coord_spans);
    
}



####
sub add_coordSet {
    my $self = shift;
    my ($acc, $end5, $end3) = @_;
    my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
    my $node = CoordSet_node->new($acc, $lend, $rend);
    push (@{$self->{node_list}}, $node);
}




####
sub build_clusters {
    my $self = shift;
    my $node_list_aref = $self->{node_list};
    @{$node_list_aref} = sort {$a->{lend}<=>$b->{lend}} @{$node_list_aref}; #sort by lend coord.
    ## set indices
    for (my $i = 0; $i <= $#{$node_list_aref}; $i++) {
        $node_list_aref->[$i]->{myIndex} = $i;
    }
    
    my @clusters;
    my $first_node = $node_list_aref->[0];
    my $start_pos = 0;
    my ($exp_left, $exp_right) = ($first_node->{lend}, $first_node->{rend});
    print $first_node->{acc} . " ($exp_left, $exp_right)\n" if $SEE;
    for (my $i = 1; $i <= $#{$node_list_aref}; $i++) {
        my $curr_node = $node_list_aref->[$i];
        my ($lend, $rend) = ($curr_node->{lend}, $curr_node->{rend});
        print $curr_node->{acc} . " ($lend, $rend)\n" if $SEE;
        if ($exp_left <= $rend && $exp_right >= $lend) { #overlap
            $exp_left = &min($exp_left, $lend);
            $exp_right = &max($exp_right, $rend);
            print "overlap. New expanded coords: ($exp_left, $exp_right)\n" if $SEE;
        } else {
            print "No overlap; Creating cluster: " if $SEE;
            my @cluster;
            for (my $j=$start_pos; $j < $i; $j++) {
                my $acc = $node_list_aref->[$j]->{acc};
                push (@cluster, $acc);
                print "$acc, " if $SEE;
            }
            push (@clusters, [@cluster]);
            $start_pos = $i;
            ($exp_left, $exp_right) = ($lend, $rend);
            print "\nResetting expanded coords: ($lend, $rend)\n" if $SEE;
        }
    }
    
    print "# Adding final cluster.\n" if $SEE;
    if ($start_pos != $#{$node_list_aref}) {
        print "final cluster: " if $SEE;
        my @cluster;
        for (my $j = $start_pos; $j <= $#{$node_list_aref}; $j++) {
            my $acc = $node_list_aref->[$j]->{acc};
            print "$acc, " if $SEE;
            push (@cluster, $acc);
        }
        push (@clusters, [@cluster]);
        print "\n" if $SEE;
    } else {
        my $acc = $node_list_aref->[$start_pos]->{acc};
        push (@clusters, [$acc]);
        print "adding final $acc.\n" if $SEE;
    }
    return (@clusters);
}

sub min {
    my (@x) = @_;
    @x = sort {$a<=>$b} @x;
    my $min = shift @x;
    return ($min);
}

sub max {
    my @x = @_;
    @x = sort {$a<=>$b} @x;
    my $max = pop @x;
    return ($max);
}

#################################################################
package CoordSet_node;
use strict;

sub new {
    my $packagename = shift;
    my ($acc, $lend, $rend) = @_;
    my $self = { acc=>$acc,
                 lend=>$lend,
                 rend=>$rend,
                 myIndex=>undef(),
                 overlapping_indices=>[]
                 };
    bless ($self, $packagename);
    return ($self);
}

1; #EOM




