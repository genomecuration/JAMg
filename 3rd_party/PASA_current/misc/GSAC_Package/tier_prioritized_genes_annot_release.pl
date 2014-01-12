#!/usr/bin/env perl

use strict;
use warnings;

use lib ($ENV{EUK_MODULES});
require "overlapping_nucs.ph";



my $usage = "usage: $0 gff3_filenames_list_file BOTH_STRANDS|SEPARATE_STRANDS [MAX_PERCENT_OVERLAP=20] [ids_to_exclude_file]\n\n";

my $file_list_file = $ARGV[0] or die $usage;
my $strand_type = $ARGV[1];
my $MAX_PERCENT_LENGTH_OVERLAP = $ARGV[2] || 20;
my $ids_to_exclude_file = $ARGV[3];

unless ($strand_type eq 'BOTH_STRANDS' || $strand_type eq 'SEPARATE_STRANDS') {
	die $usage;
}


my @annot_files =  `cat $file_list_file`;
foreach my $annot_file (@annot_files) {
	$annot_file =~ s/\s//g; # rid whitespace in filenames.
}

my %exclude_from_consideration;
if (defined($ids_to_exclude_file) && -s $ids_to_exclude_file) {
	open (my $fh, $ids_to_exclude_file) or die $!;
	while (<$fh>) {
		while (/(\S+)/g) {
			$exclude_from_consideration{$1} = 1;
		}
	}
	close $fh;
}



my %contig_id_to_selected_genes;

foreach my $file (@annot_files) {
	if ($file =~ /^\#/) { next; }
	
	my @gene_structs = &parse_gff3_file($file);
	
	foreach my $gene (@gene_structs) {
		unless ($exclude_from_consideration{ $gene->{ID} }) {
			&try_add_gene($gene);
		}
	}
}


## report results:
foreach my $contig (sort keys %contig_id_to_selected_genes) {
	
	my @genes = sort {$a->{lend}<=>$b->{lend}} @{$contig_id_to_selected_genes{$contig}};
	
	foreach my $gene (@genes) {
		
		my ($tier, $gene_ID, $lend, $rend) = ($gene->{tier}, $gene->{ID}, $gene->{lend}, $gene->{rend});
		
		my $gene_info = $gene->{gene_info};
		
		print "$contig\t$tier\t$gene_ID\t$lend\t$rend\t$gene_info\n";
	}

}

exit(0);






####
sub try_add_gene {
	my ($struct) = @_;
	
	my $contig = $struct->{contig};
	my $gene_ID = $struct->{ID};
	my $lend = $struct->{lend};
	my $rend = $struct->{rend};
	my $orient = $struct->{orient};

	my @segments = @{$struct->{segments_aref}};
	my $segments_length = $struct->{segments_length};
	


	if ($strand_type eq 'SEPARATE_STRANDS') {
		$contig = "$contig$orient";
	}
	
	my $contig_genes_aref = $contig_id_to_selected_genes{$contig};
	
	unless (ref $contig_genes_aref) {
		$contig_id_to_selected_genes{$contig} = [ $struct ];
		return;
	}

	## check for overlap:
	
	foreach my $gene (@$contig_genes_aref) {
		my ($gene_lend, $gene_rend) = ($gene->{lend}, $gene->{rend});
		if (&overlap($gene_lend, $gene_rend, $lend, $rend)) {
			
			my $num_overlapping_nucs = 0;
			my @gene_segments = @{$gene->{segments_aref}};
			foreach my $segment (@segments) {
				my ($seg_lend, $seg_rend) = @$segment;
				foreach my $other_segment (@gene_segments) {
					my ($other_lend, $other_rend) = @$other_segment;
					if (&overlap($seg_lend, $seg_rend, $other_lend, $other_rend)) {
						my $overlapping_nucs = &nucs_in_common($seg_lend, $seg_rend, $other_lend, $other_rend);
						$num_overlapping_nucs += $overlapping_nucs;
					}
				}
			}
			
			my $percentA = $num_overlapping_nucs / $segments_length * 100;
			my $percentB = $num_overlapping_nucs / $gene->{segments_length} * 100;
			
			if ($percentA > $MAX_PERCENT_LENGTH_OVERLAP || $percentB > $MAX_PERCENT_LENGTH_OVERLAP) {
				
				# Sorry, too much overlap, not adding it.
				#print "Sorry, $gene_ID overlaps " . $gene->{ID} . "\n";
				return;
			}
			
		}
	}
	
	## if got here, didn't overlap any genes on the current contig
	# add it
	push (@$contig_genes_aref, $struct);
	
	return;
}


####
sub parse_gff3_file {
	my ($gff3_file) = @_;

	my %ID_to_coords;
	my %parent_to_children_IDs;

	## parsing only gene, mRNA, and CDS.
	open (my $fh, $gff3_file) or die "Error, cannot open file $gff3_file";
	while (<$fh>) {
		if (/^\#/) { next; }
		unless (/\w/) { next; }
		chomp;
		my @x = split (/\t/);
		my ($contig, $feat_type, $lend, $rend, $orient, $gene_info) = ($x[0], $x[2], $x[3], $x[4], $x[6], $x[8]);

		my @contig_txt = split (/\s+/, $contig);
		$contig = $contig_txt[0];
		
		unless ($feat_type =~ /mRNA|CDS|gene|exon/) { next; }
		
		$gene_info =~ /ID=([^; ]+)/ or die "Error, cannot extract ID info from $_";
		
		my $feature_ID = $1;

		$ID_to_coords{$feature_ID} = { tier => $gff3_file,
									   type => $feat_type,
									   lend => $lend,
									   rend => $rend,
									   contig => $contig,
									   ID => $feature_ID,
									   gene_info => $gene_info,
									   orient => $orient,
								   };
		
		if ($gene_info =~ /Parent=([^; ]+)/) {
			my $parent_ID = $1;
			push (@{$parent_to_children_IDs{$parent_ID}}, $feature_ID);
		}
	}
	close $fh;


	## get the genes:
	my @gene_structs;
	foreach my $struct (values %ID_to_coords) {
		if ($struct->{type} eq 'gene') {
			push (@gene_structs, $struct);
		}
	}
	
	
	my @ret_genes;

	foreach my $gene (@gene_structs) {
		my $gene_ID = $gene->{ID};
		
		my $gene_info = $gene->{gene_info};

		my @mRNA_IDs = @{$parent_to_children_IDs{$gene_ID}};
		
		if (scalar (@mRNA_IDs) > 1) {
		    die "note, $gene_ID has multiple mRNAs: @mRNA_IDs\n";
		}
		elsif (! @mRNA_IDs) {
			die "Error, gene $gene_ID lacks an mRNA";
		}
		
		my $mRNA_ID = shift @mRNA_IDs;
		my @children = @{$parent_to_children_IDs{$mRNA_ID}};

		my @cds_coordsets;
		
		my @coords;
		foreach my $child (@children) {
			my $struct = $ID_to_coords{$child};
			if ($struct->{type} eq 'CDS' || ($gene_info =~ /pseudogene/ && $struct->{type} eq 'exon') ) {  # use exon coords for pseudogenes, since annotated cds is bogus.
				push (@coords, $struct->{lend}, $struct->{rend});
				push (@cds_coordsets, [$struct->{lend}, $struct->{rend}]);
			}
		}
		
		@coords = sort {$a<=>$b} @coords;
		my $gene_lend = shift @coords;
		my $gene_rend = pop @coords;
		
		## reset coords based on CDS span
		$gene->{lend} = $gene_lend;
		$gene->{rend} = $gene_rend;
		
		
		$gene->{segments_aref} = \@cds_coordsets;
		my $segment_lengths = 0;
		foreach my $segment (@cds_coordsets) {
			my ($lend, $rend) = @$segment;
			my $len = $rend - $lend + 1;
			$segment_lengths += $len;
		}
		$gene->{segments_length} = $segment_lengths;
		
		push (@ret_genes, $gene);
	}


	return (@ret_genes);
}


####
sub overlap {
	my ($lendA, $rendA, $lendB, $rendB) = @_;

	if ($lendA <= $rendB && $rendA >= $lendB) {
		return (1);
	}
	else {
		return(0);
	}
}


