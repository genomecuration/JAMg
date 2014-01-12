#!/usr/bin/env perl

use strict;
use warnings;
use GD;
use Bio::TreeIO;
use lib ($ENV{EUK_MODULES});
use ColorGradient;
use BHStats;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case bundling);

my $DEBUG = 0;

my $gff3_file = undef;
my $seq_lengths_file = undef;
my $tree_file = undef;
my $min_genes_per_scaffold = 1;
my $canvas_width = 500;
my $canvas_height = 750;
my $suggested_number_of_tiers = 4;  # actual number of tiers depends on sequence lengths.


my $usage = <<_EOUSAGE_;

################################################################################
#   Required:
#   --gff3            gff3 file containing all genes to plot
#   --seqLengths      file containing lengths of genome sequences (format:      )

#   Optional:
#   --treeFile        tree file in newhampshire format        
#   --min_genes_per_scaffold     minimum number of genes required per scaffold (crank up to focus on tandem gene dups). (default 1)
#   --canvas_width                default 500
#   --canvas_height               default 750
#   --num_vertical_tiers          default 4   (a suggestion)
#
################################################################################

_EOUSAGE_
	;

&GetOptions (
			 'gff3_file=s' => \$gff3_file,
			 'seqLengths=s' => \$seq_lengths_file,
			 'treeFile=s' => \$tree_file,
			 'min_genes_per_scaffold=i' => \$min_genes_per_scaffold,
			 'canvas_width=i' => \$canvas_width,
			 'canvas_height=i' => \$canvas_height,
			 'num_vertical_tiers=i' => \$suggested_number_of_tiers,
			 );


unless (defined($gff3_file) && defined($seq_lengths_file) && -f $gff3_file && -f $seq_lengths_file) {
	die $usage;
}

##################################

##################################

## Constants:
# fixed sizes:

my $top_margin = 100;
my $bottom_margin = 50;
my $left_margin = 50;
my $right_margin = 50;

my $tree_depth_pixels = 10; #default, will be adjusted based on num tiers, etc.

my $neighboring_scaffold_spacer = 20;

# ratios:

my $chromo_box_panel_width_ratio = 0.1;  # rest for tree.
my $gene_tick_chromo_box_ratio = 0.5;


# derived values:
my $image_height = $canvas_height + $top_margin + $bottom_margin;
my $image_width = $canvas_width + $left_margin + $right_margin;

## color scale positioning (note this is written in the upper right margin).
my $color_map_width = 100;
my $color_map_height = 30;
my $color_map_right_margin = 20;
my $color_map_top_margin = 20;

# color map scale params:
my $color_scale_lowest_branch_length = 0;
my $color_scale_highest_branch_length = 0.005;


# misc Globals:
my %colors;
my $scaffold_color = 'black';
my $gene_color = 'black';




my @color_gradient;

my %gene_to_scaffold;
my $max_branch_length;
my $min_branch_length;
my $avg_branch_length;

my $color_using_branch_length = 0; # alternatively, colors based on tree position.
my %node_ID_to_color;


my %mRNA_to_gene_ID; # track child/parent relationships in gff3 file

main: {
	my %scaffold_lengths = &parse_scaffold_lengths($seq_lengths_file);

	my %scaffold_to_genes = &parse_genes_from_gff3_file($gff3_file);


	## remove those scaffolds that have less than two genes:
	foreach my $scaffold (keys %scaffold_to_genes) {
		my $num_genes = scalar (@{$scaffold_to_genes{$scaffold}});
		if ($num_genes < $min_genes_per_scaffold) {
			delete $scaffold_to_genes{$scaffold};;
		}
		else {
			my @genes = @{$scaffold_to_genes{$scaffold}};
			foreach my $gene (@genes) {
				my $acc = $gene->{acc};
				$gene_to_scaffold{$acc} = $scaffold;
			}
		}
	}
	
	my @scaffolds = reverse sort {$#{$scaffold_to_genes{$a}}<=>$#{$scaffold_to_genes{$b}}} keys %scaffold_to_genes;

	my $tree;
	if ($tree_file) {
		my $treeio = new Bio::TreeIO(-format => 'newick',
									 -file => $tree_file);                                                                                                                 
		$tree = $treeio->next_tree();  
	}
	
	
	my @tiers = &tier_scaffolds(\%scaffold_lengths, \@scaffolds, \%scaffold_to_genes, $tree);
	
	my $pixels_per_bp = &compute_smallest_pixels_per_bp(\%scaffold_lengths, \@tiers);
	
	
	my $image = new GD::Image($image_width, $image_height);
	&init_color_palette($image);

	my %gene_acc_to_XY_position;
	&draw_scaffolds_n_genes( { scaffold_lengths => \%scaffold_lengths,
							   scaffold_to_genes => \%scaffold_to_genes,
							   tiers => \@tiers,
							   image => $image,
							   pixels_per_bp => $pixels_per_bp,
							   gene_acc_to_XY_position => \%gene_acc_to_XY_position,
							   tree => $tree,
						   } );
	

	if ($tree_file) {
		&assign_tree_nodes($tree);
		unless ($color_using_branch_length) {
			&assign_colors_to_nodes_via_branching_pattern_from_root($tree);
		}
		

		# print STDERR Dumper (\%gene_acc_to_XY_position);
		
		&add_tree_to_image($tree, $image, \%gene_acc_to_XY_position);

		# &add_colormap_scale($image);

		## write the leaf node colors:
		open (my $fh, ">$gff3_file.leaf_colors.$$.txt") or die $!;
		foreach my $node ($tree->get_nodes()) {
			if ($node->is_Leaf()) {
				my $node_id = $node->id();
				my $color_rgb = $node_ID_to_color{$node_id};
				my ($r, $g, $b) = split (/,/, $color_rgb);
				my $hex_color = sprintf "#%02x%02x%02x", $r, $g, $b;
				print $fh "$node_id,$hex_color\n";
			}
		}
		close $fh;

	}
	
	print $image->png() if !$DEBUG;


	

	exit(0);

}


####
sub order_scaffolds_by_tree_position {
	my ($tree, $scaffolds_to_genes_href) = @_;

	my %gene_to_scaffold;
	foreach my $scaffold (keys %$scaffolds_to_genes_href) {
		my @gene_hrefs = @{$scaffolds_to_genes_href->{$scaffold}};
		
		foreach my $gene_href (@gene_hrefs) {
			my $acc = $gene_href->{acc};
			$gene_to_scaffold{$acc} = $scaffold;
		}
	}

	## turn nodes into list of contigs.
	
	my @nodes = $tree->get_leaf_nodes();
	my @contig_list;
	foreach my $node (@nodes) {
		my $node_ID = $node->id();
		my $scaffold = $gene_to_scaffold{$node_ID} or die "Error, cannot find scaffold based on node ID: $node_ID ";
		push (@contig_list, $scaffold);
	}

	## group into adjacent counts.
	my @group_counts;
	foreach my $contig (@contig_list) {
		unless (@group_counts) {
			push (@group_counts, [$contig, 1]);
			next;
		}
		
		my $prev_group = $group_counts[$#group_counts];
		if ($contig eq $prev_group->[0]) {
			$prev_group->[1]++;
		}
		else {
			push (@group_counts, [$contig, 1]);
		}
	}

	## determine max counts for each contig.
	my %max_counts;
	foreach my $group (@group_counts) {
		my ($contig, $count) = @$group;
		if ((!exists $max_counts{$contig}) || $count > $max_counts{$contig}) {
			$max_counts{$contig} = $count;
		}
	}
	
	## return in order of max count and non-redundant list.
	my %seen;
	my @ret_contigs;
	
	foreach my $group (@group_counts) {
		my ($contig, $count) = @$group;
		if ($seen{$contig}) { next; }
		if ($count == $max_counts{$contig}) {
			push (@ret_contigs, $contig);
			$seen{$contig} = 1;
		}
	}

	return (@ret_contigs);
	
}
	



####
sub draw_scaffolds_n_genes {
	my ($inputs_href) = @_;

	# unwrap inputs:
	my $scaffold_lengths_href = $inputs_href->{scaffold_lengths};
	my $scaffold_to_genes_href = $inputs_href->{scaffold_to_genes};
	my $tiers_aref = $inputs_href->{tiers};
	my $image = $inputs_href->{image};
	my $pixels_per_bp = $inputs_href->{pixels_per_bp};
	my $gene_acc_to_XY_position = $inputs_href->{gene_acc_to_XY_position};
	my $tree = $inputs_href->{tree};  

	## draw each scaffold and gene set one scaffold at a time:
	my $x_pos = $left_margin;
	my $y_pos = $image_height - $bottom_margin;
	
	my $num_tiers = scalar (@{$tiers_aref});
	my $tier_panel_width = int($canvas_width / $num_tiers);
	
	my $chromo_panel_width = $chromo_box_panel_width_ratio * $tier_panel_width;
	my $gene_width = int($chromo_panel_width * $gene_tick_chromo_box_ratio);
	my $chromo_width = $chromo_panel_width - $gene_width;
	
	my $remainder_width = $tier_panel_width - $chromo_panel_width;


	my  $tier_count = 0;
	foreach my $tier (@$tiers_aref) {
		$tier_count++;
		
		
		my $tier_x1_pos = $x_pos;
		my $tier_x2_pos = $x_pos + $tier_panel_width;
		my $tier_y_pos = $y_pos;
		# $image->rectangle($tier_x1_pos, $top_margin, $tier_x2_pos, $top_margin+$canvas_height, $colors{black});
		
				
		my $chromo_x1_pos = $tier_x1_pos;
		my $chromo_x2_pos = $tier_x1_pos + $chromo_width;
		
		my $gene_x1_pos = $chromo_x2_pos;
		my $gene_x2_pos = $gene_x1_pos + $gene_width;
		
		foreach my $scaffold (@$tier) {
			print STDERR "Drawing: tier($tier_count), scaffold: $scaffold\n";
			my $scaffold_length = $scaffold_lengths_href->{$scaffold};
			
			my $scaffold_height_pixels = int($scaffold_length * $pixels_per_bp);
			my $y2_pos = $tier_y_pos - $scaffold_height_pixels;
			$image->filledRectangle($chromo_x1_pos, $y2_pos, $chromo_x2_pos, $tier_y_pos, $colors{$scaffold_color});
			#$image->rectangle($chromo_x1_pos, $y2_pos, $chromo_x2_pos, $tier_y_pos, $colors{black});

			## draw the genes:
			foreach my $gene (@{$scaffold_to_genes_href->{$scaffold}}) {
				my ($acc, $lend, $rend) = ($gene->{acc}, $gene->{lend}, $gene->{rend});
				my $midpt = ($lend + $rend) / 2;
				my $gene_y_pos = $tier_y_pos - int($midpt * $pixels_per_bp);
				$gene_acc_to_XY_position->{$acc} = [$gene_x2_pos, $gene_y_pos];
				
				$image->line($gene_x1_pos, $gene_y_pos, $gene_x2_pos, $gene_y_pos, $colors{$gene_color});
				
			}
			
			## set up for next scaffold draw:
			$tier_y_pos -= $scaffold_height_pixels + $neighboring_scaffold_spacer;
		}
		
		$x_pos += $tier_panel_width;
	}

	
	if ($tree) {
		my $height = $tree->height();
		$tree_depth_pixels = int($remainder_width / $height);
	}
	
	return;
}


####
sub init_color_palette {
	my ($image) = shift;
	$colors{white} = $colors{"255,25,255"} =  $image->colorAllocate(255,255,255);
	$colors{black} = $colors{"0,0,0"} = $image->colorAllocate(0,0,0);
	$colors{red} = $colors{"255,0,0"} = $image->colorAllocate(255,0,0);
	$colors{blue} = $colors{"0,0,255"} = $image->colorAllocate(0,0,255);
	$colors{green} = $colors{"0,255,0"} = $image->colorAllocate(0,255,0);


	## build the color gradient w/ 200 colors:
	my @rgb_combos = &ColorGradient::get_RGB_gradient(200);
	foreach my $rgb_combo (@rgb_combos) {
		my ($r, $g, $b) = @$rgb_combo;
		my $color_token = join (",", $r, $g, $b);
		unless ($colors{$color_token}) { 
			$colors{$color_token} = $image->colorAllocate($r, $g, $b);
		}
		push (@color_gradient, $color_token);
	}
	
	return;
}


####
sub parse_scaffold_lengths {
	my ($seq_lengths_file) = @_;

	my %scaff_lengths;

	open (my $fh, $seq_lengths_file) or die "Error, cannot open file $seq_lengths_file";
	while (<$fh>) {
		chomp;
		my ($length, $acc) = split (/\s+/);
		$scaff_lengths{$acc} = $length;
	}
	close $fh;
	
	return (%scaff_lengths);
}

####
sub parse_genes_from_gff3_file {
	my ($gff3_file) = @_;

	my %scaffold_to_genes;

	open (my $fh, $gff3_file) or die "Error, cannot open file $gff3_file";
	while (<$fh>) {
		chomp;
		if (/^\#/) { next; }
		unless (/\w/) { next; }

		my @x = split (/\t/);
		my ($contig, $type, $lend, $rend, $gene_info) = ($x[0], $x[2], $x[3], $x[4], $x[8]);
		unless ($type eq 'gene' || $type eq 'mRNA') { next; }

		$gene_info =~ /ID=([^;]+)/ or die "Error, cannot parse gene ID from gene_info";
		my $gene_ID = $1;

		if ($type eq 'mRNA') {
			$gene_info =~ /Parent=([^;]+)/;
			my $parent_id = $1 or die "Error, cannot find parent ID from mRNA entry: $_";
			$mRNA_to_gene_ID{$gene_ID} = $parent_id; # note gene_ID here is actually the mRNA ID. 
		}
		
		unless ($type eq 'mRNA') { 
			next;
		}

		push (@{$scaffold_to_genes{$contig}}, { acc => $gene_ID,
												lend => $lend,
												rend => $rend, } );
		
	}

	close $fh;
	return (%scaffold_to_genes);
}


####
sub tier_scaffolds {
	my ($scaffold_lengths_href, $scaffolds_aref, $scaffolds_to_genes_href, $tree) = @_;
	my $num_scaffolds = scalar @$scaffolds_aref;                                                                                                                                                                 
    
	my @ordered_scaffolds = @$scaffolds_aref;
	
	if ($tree) {
		@ordered_scaffolds = &order_scaffolds_by_tree_position($tree, $scaffolds_to_genes_href);
	}
	

	print STDERR "Order of tiered scaffolds, bottom-up, left to right:\n"
		. join ("\n", @ordered_scaffolds) . "\n\n";
	

	my $sum_scaffold_length = 0;
    foreach my $scaffold (@ordered_scaffolds) {                                                                                                                                                                    
        my $len = $scaffold_lengths_href->{$scaffold} or die "Error, no scaff length for $scaffold";                                                                                                              
        $sum_scaffold_length += $len;                                                                                                                                                                      
    }       
	
	my $suggested_max_tier_bp = $sum_scaffold_length / $suggested_number_of_tiers;
	
	# print STDERR "$sum_scaffold_length / $suggested_number_of_tiers = suggested_max_tier_bp = $suggested_max_tier_bp\n";
	
	my @tiers = ([]);
	my $sum_tier_length = 0;

	foreach my $scaffold (@ordered_scaffolds) {
		my $scaff_len = $scaffold_lengths_href->{$scaffold};
		
		# print STDERR "scaff_len($scaffold) = $scaff_len\n";
		my $tier = $tiers[$#tiers];
		push (@$tier, $scaffold);
		
		$sum_tier_length += $scaff_len;
		if ($sum_tier_length > $suggested_max_tier_bp) {
			# start a new tier for next one.
			push (@tiers, []);
			print STDERR "-making new tier. (sum= $sum_tier_length)\n";
			$sum_tier_length = 0;
		}
	}

	## rid the last one if it's empty:
	unless (@{$tiers[$#tiers]}) {
		pop @tiers;
	}
	
	return (@tiers);
}

####
sub compute_smallest_pixels_per_bp {
	my ($scaffold_lengths_href, $tiers_aref) = @_;
	
	my @pixels_per_bp;
	
	foreach my $tier (@$tiers_aref) {
		my @scaffs = @$tier;
		my $num_scaffs = scalar (@scaffs);
		my $sum_length = 0;
		foreach my $scaff (@scaffs) {
			my $len = $scaffold_lengths_href->{$scaff};
			$sum_length += $len;
		}
		my $pixels_per_bp = ($canvas_height - ($num_scaffs -1)*$neighboring_scaffold_spacer) / $sum_length;
		push (@pixels_per_bp, $pixels_per_bp);
	}

	@pixels_per_bp = sort {$a<=>$b} @pixels_per_bp;
	my $smallest_pixels_per_bp = shift @pixels_per_bp;
	return ($smallest_pixels_per_bp);
}


####
sub assign_tree_nodes {
	my ($tree) = @_;
	
	my $root_node = $tree->get_root_node();
	$root_node->id("myRooot");
	
	## assign identifiers to nodes that lack them (ie. internal nodes).
	my $x=0;
		
	my @branch_lengths;
	
	my @nodes = $tree->get_nodes();
	foreach my $node (@nodes) {
		my $id = $node->id();
		unless ($id =~ /\w/) {
			$x++;
			$node->id("$x");
		}
		my $num_descendents = scalar ($node->each_Descendent());
		#print "node " . $node->id() . " has $num_descendents descendents\n";
	
		my $branch_length = $node->branch_length();
		
		if (defined($branch_length)) {
			push (@branch_lengths, $branch_length);
		}
		
	}
		
	## gather branch length stats:
	$avg_branch_length = &BHStats::avg(@branch_lengths);
	$max_branch_length = &BHStats::max(@branch_lengths);
	$min_branch_length = &BHStats::min(@branch_lengths);
	
	
	print STDERR "Branch length stats: min($min_branch_length), max($max_branch_length), avg($avg_branch_length)\n";

	return;
}
####
sub add_tree_to_image {
	my ($tree, $image, $gene_acc_to_XY_position_href) = @_;
		
	my $root_node = $tree->get_root_node();
	my @root_descendents = $root_node->each_Descendent(); # tree is unrooted, so three descendents.
	
	if (scalar @root_descendents == 2) {
		## rooted tree
		&recursively_connect_node_pairs($root_node, $image, $gene_acc_to_XY_position_href);
	}
	else {
		## unrooted tree
		foreach my $descendent (@root_descendents) {
			unless ($descendent->is_Leaf()) {
				&recursively_connect_node_pairs($descendent, $image, $gene_acc_to_XY_position_href);
			}
		}
	}
	return;
}


####
sub recursively_connect_node_pairs {
	my ($node, $image, $gene_acc_to_XY_position_href) = @_;

	my $node_ID = $node->id();
	
	if ($node->is_Leaf()) {
		die "Error, connect node pairs called on a leaf node";
	}

	my @descendents = $node->each_Descendent();
	# print STDERR "node: $node_ID, @descendents\n";
	if (scalar @descendents != 2) {
		die "Error, didn't get exactly two descendent nodes from a non-leaf node";
	}
	
	my ($nodeA, $nodeB) = @descendents;

	my $nodeA_ID = $nodeA->id();
	my $nodeB_ID = $nodeB->id();

	my $branch_length_A = $nodeA->branch_length();
	my $branch_length_B = $nodeB->branch_length();
	
	if ((! $nodeA->is_Leaf()) && !defined $gene_acc_to_XY_position_href->{$nodeA_ID}) {
		&recursively_connect_node_pairs($nodeA, $image, $gene_acc_to_XY_position_href);
	}

	if ( (! $nodeB->is_Leaf()) &&  !defined $gene_acc_to_XY_position_href->{$nodeB_ID}) {
		&recursively_connect_node_pairs($nodeB, $image, $gene_acc_to_XY_position_href);
	}

	my $coordsetA = $gene_acc_to_XY_position_href->{$nodeA_ID};
	my $coordsetB =  $gene_acc_to_XY_position_href->{$nodeB_ID};

	# print STDERR "Trying to connect [$nodeA_ID] to [$nodeB_ID] using ($coordsetA, $coordsetB)\n";
	
	if (ref $coordsetA && ref $coordsetB) {
		my ($AX, $AY) = @$coordsetA;
		my ($BX, $BY) = @$coordsetB;
		
		my $midX = int( ($AX+$BX)/2);
		my $midY = int( ($AY+$BY)/2);
		$midX += $tree_depth_pixels;

		
		my ($a_branchlength_colortoken, $b_branchlength_colortoken);
		if ($color_using_branch_length) {
			$a_branchlength_colortoken = &get_branchlength_colortoken($branch_length_A);
			$b_branchlength_colortoken = &get_branchlength_colortoken($branch_length_B);
		}
		else {
			## color based on tree branching pattern:
			$a_branchlength_colortoken = $node_ID_to_color{ $nodeA_ID };
			$b_branchlength_colortoken = $node_ID_to_color{ $nodeB_ID };
		}
		
		## draw line:
		$image->line($AX, $AY, $midX, $midY, $colors{$a_branchlength_colortoken});
		$image->line($BX, $BY, $midX, $midY, $colors{$b_branchlength_colortoken});
		# print STDERR "drawing tree lines.\n";
		$gene_acc_to_XY_position_href->{$node_ID} = [$midX, $midY];
	}
	elsif ( (ref $coordsetA) && (!ref $coordsetB) && $nodeB->is_Leaf()) {
		## assign internal node the coordset of A
		$gene_acc_to_XY_position_href->{$node_ID} = $coordsetA;
	}
	elsif ( (ref $coordsetB) && (! ref $coordsetA) && $nodeA->is_Leaf()) {
		$gene_acc_to_XY_position_href->{$node_ID} = $coordsetB;
	}
	else {
		$gene_acc_to_XY_position_href->{$node_ID} = -1; # not an array ref, also not !defined. :)
	}

	return;
}

####
sub get_branchlength_colortoken {
	my ($branchlength) = @_;

	my $max_color_index = $#color_gradient;

	my $index;
	
	if ($branchlength <= $color_scale_lowest_branch_length) {
		$index = 0;
	}
	elsif ($branchlength >= $color_scale_highest_branch_length) {
		$index = $max_color_index;
	}
	else {
		# branch length is somewhere in between.
		my $low_high_delta = $color_scale_highest_branch_length - $color_scale_lowest_branch_length;
		my $delta_low = $branchlength - $color_scale_lowest_branch_length;
		my $ratio = $delta_low / $low_high_delta;
		$index = int($ratio * $max_color_index);
	}
	
	
	return ($color_gradient[$index]);
}


####
sub add_colormap_scale {
	my ($image) = @_;

	## establish color map box boundaries:
	my $x1 = $image_width - $color_map_right_margin - $color_map_width;
	my $x2 = $image_width - $color_map_right_margin;
	my $y1 = $color_map_top_margin;
	my $y2 = $color_map_top_margin + $color_map_height;

	## go through every color:
	my $num_colors = scalar (@color_gradient);
	my $color_pixel_increments = $color_map_width / $num_colors;
	
	## rainbow color gradient:
	#for (my $i = 0; $i < $num_colors; $i++) {
	#	my $color_token = $color_gradient[$i];
	##	my $left_x = int($x1 + $i * $color_pixel_increments);
	#	my $right_x = $left_x + $color_pixel_increments;
	#	$image->filledRectangle($left_x, $y1, $right_x, $y2, $colors{$color_token});
	#}
		
	my $max_min_delta = $color_scale_highest_branch_length - $color_scale_lowest_branch_length;
	my $length_increment = $max_min_delta / $num_colors;
	my $init_len = $min_branch_length;
   
	 for (my $i = 0; $i < $num_colors; $i++) {                                                                                                                         
		 my $branch_length = $init_len + $i * $length_increment;
		 
		 my $color_token = &get_branchlength_colortoken($branch_length);
		 my $left_x = int($x1 + $i * $color_pixel_increments);                                                        
		 my $right_x = $left_x + $color_pixel_increments;                                                                                                               
		 $image->filledRectangle($left_x, $y1, $right_x, $y2, $colors{$color_token});                                                                                   
	 }               		
	
	
	return;
}

####
sub assign_colors_to_nodes_via_branching_pattern_from_root 	{
	my ($tree) = @_;

	my $root_node = $tree->get_root_node();
	
	&recursively_color_clades($root_node, [@color_gradient]);
	
	return;
}


####
sub recursively_color_clades {
	my ($node, $colors_avail_aref) = @_;

	my @descendents = $node->each_Descendent();
	
	
	my $sum_branch_lengths = &sum_branch_lengths($node);
	
	my $num_colors = scalar @$colors_avail_aref;

	my $prev_lower_index = 0;
	for (my $i = 0; $i <= $#descendents; $i++) {
		my $descendent = $descendents[$i];
		
		my $sum_descendent_branch_lengths = &sum_branch_lengths($descendent);

		# print STDERR "descendent($i): $sum_descendent_branch_lengths / $sum_branch_lengths\n";

		my $ratio_total;
		if ($sum_branch_lengths > 0) {
			$ratio_total = $sum_descendent_branch_lengths / $sum_branch_lengths;
		}
		else {
			$ratio_total = 1;
		}
		
		my $relative_index = int($ratio_total * $num_colors);
		my $high_index = $prev_lower_index + $relative_index;
		if ($high_index > $num_colors -1) {
			$high_index = $num_colors - 1;
		}
		
		my @colors_to_consider = @$colors_avail_aref[$prev_lower_index..$high_index];
		my $midpt_color = $colors_to_consider[ int($#colors_to_consider/2 + 4/9) ];
		
		# assign color:
		$node_ID_to_color{$descendent->id()} = $midpt_color;
		
		unless ($descendent->is_Leaf()) {
			&recursively_color_clades($descendent, [@colors_to_consider]);
		}
		
		$prev_lower_index = $high_index;
	}
		
	
	return;
}
		
		
####
sub sum_branch_lengths {
	my ($node) = @_;

	my $sum = 0;
	
	foreach my $descendent ($node, $node->get_all_Descendents()) {
		$sum += $descendent->branch_length() || 0;
	}


	return ($sum);
}

