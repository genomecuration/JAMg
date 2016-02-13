#!/usr/bin/env perl

use strict;
use warnings;
use FindBin;
use lib ("$FindBin::Bin/../..//PerlLib");
use File::Basename;

use Fasta_reader;
use GTF_to_geneobjs;
use GFF3_to_geneobjs;
use Getopt::Long qw(:config no_ignore_case bundling);
use Data::Dumper;
use BHStats;

umask(0000);

my $usage = <<_EOUSAGE_;

####################################################################################################
#  
#  --annot_file              genes file in corresponding format 
#  --format                  GTF or GFF3
#  -v                        verbose
#  -X                        eXport features
#####################################################################################################

_EOUSAGE_

    ;

my ($annot_file, $format, $help, $EXPORT_FLAG);

our $SEE;

&GetOptions ( 
			  "annot_file=s" => \$annot_file,
			  "format=s" => \$format,
			  "help|h" => \$help,
			  "v" => \$SEE,
			  "X" => \$EXPORT_FLAG,
              );

if ($help) { die $usage; }
unless ($annot_file && $format =~ /^(GTF|GFF3)$/) { 
    die $usage;
}

my $core_annot_filename = basename($annot_file);
my ($exonfh, $intronfh, $genicfh, $intergenicfh, $utrfh);
if ($EXPORT_FLAG) {
	open ($exonfh, ">$core_annot_filename.exons") or die $!;
	open ($intronfh, ">$core_annot_filename.introns") or die $!;
	open ($genicfh, ">$core_annot_filename.genes") or die $!;
	open ($intergenicfh, ">$core_annot_filename.intergenic") or die $!;
	open ($utrfh, ">$core_annot_filename.utr") or die $!;
}

# stats interested in:
my $gene_count = int(0);
my $mRNA_count = int(0);
my $alt_spliced_gene_count = int(0);
my $intron_containing_gene_count = int(0);
my $got_5utr = int(0);
my $got_3utr = int(0);
my $unique_exon_count = int(0);
my $unique_cds_count = int(0);
my $unique_intron_count = int(0);
my $unique_5utr_count = int(0);
my $unique_3utr_count = int(0);
my $alt_splice_diff_CDSs_count = int(0);
my $diff_splice_CDS_count = int(0);

my $num_intergenic_regions = int(0);
my $sum_intergenic_lengths = int(0);
my $sum_gene_lengths = int(0);
my $sum_intron_lengths = int(0);
my $sum_exon_lengths = int(0);
my $sum_5utr_lengths = int(0);
my $sum_3utr_lengths = int(0);
my ($got_5utr_long, $got_3utr_long)=(int(0),int(0));


main: {
        
    ## get the genes:
    my %gene_id_to_gene;
    my $seqname_map_href;
	
    if ($format eq 'GTF') {
        print STDERR "-processing GTF files\n" if $SEE;
        $seqname_map_href = &GTF_to_geneobjs::parse_file($annot_file, \%gene_id_to_gene);
	}
    else {
        print STDERR "-processing GFF3 files\n" if $SEE;
        $seqname_map_href = &GFF3_to_geneobjs::parse_file($annot_file, \%gene_id_to_gene);
	}

	
	foreach my $contig (keys %$seqname_map_href) {
        
        print STDERR "// processing $contig\n" if $SEE;
		my $gene_ids_aref = $seqname_map_href->{$contig};
		
		my %exons;
		my %cdss;
		my %introns;
		my %utr5;
		my %utr3;

		my @gene_spans;

		foreach my $gene_id (@$gene_ids_aref) {
                
			my $gene_obj = $gene_id_to_gene{$gene_id} or die "Error, no gene retrieved from reference via $gene_id";
                
			my ($gene_lend, $gene_rend) = sort {$a<=>$b} $gene_obj->get_coords();
			push (@gene_spans, [$gene_lend, $gene_rend]);
			
			$gene_count++;
		
			
			my %complete_CDS_tokens; # tracking the full CDS structure; ignores alt isoforms w/ splice vars in UTRs only
	
			if ($gene_obj->get_additional_isoforms()) {
				$alt_spliced_gene_count++;
			}


			my ($got_intron_containing_gene_flag);
						
			foreach my $isoform ($gene_obj, $gene_obj->get_additional_isoforms()) {
				$mRNA_count++;

				my @all_cds_coords;
				my @exons = $isoform->get_exons();

				foreach my $exon (@exons) {
					my $exon_token = join ("_", $exon->get_coords());
					$exons{$exon_token} = 1;

					if (my $cds = $exon->get_CDS_exon_obj()) {
						my $cds_token = join ("_", $cds->get_coords());
						$cdss{$cds_token} = 1;
					
						push (@all_cds_coords, $cds->get_coords());
					}
				}

				if (my @introns = $isoform->get_intron_coordinates()) {
					$got_intron_containing_gene_flag = 1;
					foreach my $intron_coordset (@introns) {
						my $intron_token = join ("_", @$intron_coordset);
						$introns{$intron_token}++;
					}
				}

				my $complete_cds_token = join ("_", sort @all_cds_coords);
				$complete_CDS_tokens{$complete_cds_token} = 1;
				
				if ( my @utrs = $isoform->get_5prime_UTR_coords()){
					$got_5utr++;
					foreach my $set (@utrs){
						my $prime5_utr_token = join ("_", @$set);
						$utr5{$prime5_utr_token}++;
					}
				}
				if ( my @utrs = $isoform->get_3prime_UTR_coords()){
					$got_3utr++;
					foreach my $set (@utrs){
						my $prime3_utr_token = join ("_", @$set);
						$utr3{$prime3_utr_token}++;
					}
				}

			}

			if ($got_intron_containing_gene_flag) {
				$intron_containing_gene_count++;
			}

			my $num_cds_tokens = scalar (keys %complete_CDS_tokens);
			
			$diff_splice_CDS_count += $num_cds_tokens;
			if ($num_cds_tokens > 1) {
				$alt_splice_diff_CDSs_count++; # count alt-spliced genes ignoring UTR diffs only.
			}
			
		} # end of foreach gene


		&analyze_intergenics($contig, \@gene_spans);
		
		# get exon and CDS stats:
		my @unique_exons = keys %exons;
		my @unique_cdss = keys %cdss;
		my @unique_introns = keys %introns;
		my @unique_5utr = keys %utr5;
		my @unique_3utr = keys %utr3;

		$unique_exon_count += scalar(@unique_exons);
		$unique_cds_count += scalar(@unique_cdss);
		$unique_intron_count += scalar(@unique_introns);
		$unique_5utr_count += scalar(@unique_5utr);
		$unique_3utr_count += scalar(@unique_3utr);

		## get sum lengths:
		foreach my $exon (@unique_exons) {
			my ($exon_lend, $exon_rend) = sort {$a<=>$b} split (/_/, $exon);
			my $exon_len =  ($exon_rend - $exon_lend) + 1;
			$sum_exon_lengths += $exon_len;
			print $exonfh "exon\t$contig\t$exon_lend\t$exon_rend\t$exon_len\n" if $EXPORT_FLAG;
		}

		foreach my $intron (@unique_introns) {
			my ($intron_lend, $intron_rend) = sort {$a<=>$b} split (/_/, $intron);
			my $intron_len = $intron_rend - $intron_lend + 1;
			$sum_intron_lengths += $intron_len;
			print $intronfh "intron\t$contig\t$intron_lend\t$intron_rend\t$intron_len\n" if $EXPORT_FLAG;
		}

		foreach my $utr (@unique_5utr){
			my ($lend, $rend) = sort {$a<=>$b} split (/_/, $utr);
			my $len = $rend - $lend + 1;
			$sum_5utr_lengths += $len;
			print $utrfh "UTR5\t$contig\t$lend\t$rend\t$len\n" if $EXPORT_FLAG;			
			next unless $len >= 30;
			$got_5utr_long++;
		}
		foreach my $utr (@unique_3utr){
			my ($lend, $rend) = sort {$a<=>$b} split (/_/, $utr);
			my $len = $rend - $lend + 1;
			$sum_3utr_lengths += $len;
			print $utrfh "UTR3\t$contig\t$lend\t$rend\t$len\n" if $EXPORT_FLAG;			
			next unless $len >= 30;
			$got_3utr_long++;
		}
		
		
	} # end of foreach contig

	
	## summarize statistics:
	
	print "\n\n";
	print &thousands($gene_count)." genes\n";
	print &thousands($mRNA_count)." mRNAs\n";
	print &thousands($unique_exon_count)." unique exons\n";
	print &thousands($unique_cds_count)." unique CDSs\n";
	print &thousands($unique_intron_count)." unique introns\n";

	printf ("\n%.1f exons per gene\n", $unique_exon_count / $gene_count);
	printf ("%.1f CDSs per gene\n", $unique_cds_count / $gene_count);
	printf ("%.1f introns per gene\n", $unique_intron_count / $gene_count);
	
	print "\n";
	print &thousands($alt_spliced_gene_count)." genes alternatively spliced\n";
	printf ("%.1f%% genes alternatively spliced\n", $alt_spliced_gene_count / $gene_count * 100);
	print "\n";
	printf ("%.2f transcripts per gene\n", $mRNA_count / $gene_count);

	print "Across all genes:\n";
	print &thousands($got_5utr)." have 5'UTR\n";
	print &thousands($got_3utr)." have 3'UTR\n";
	print &thousands($got_5utr_long)." have 5'UTR >= 30 bp\n";
	print &thousands($got_3utr_long)." have 3'UTR >= 30 bp\n";
	print "\n";
	
	if ($alt_spliced_gene_count) {
		my $genes_not_alt_spliced = $gene_count - $alt_spliced_gene_count;
		printf ("%.2f transcripts per alt-spliced gene\n", ( $mRNA_count - $genes_not_alt_spliced ) / $alt_spliced_gene_count );
	}

	if ($alt_splice_diff_CDSs_count) {	
		print &thousands($alt_splice_diff_CDSs_count)." genes alt spliced w/ altsplicing in coding regions.\n";
		my $genes_not_alt_spliced = $gene_count - $alt_splice_diff_CDSs_count;
		printf("%.2f CDS-structures per alt-spliced gene\n\n\n", ($diff_splice_CDS_count - $genes_not_alt_spliced) / $alt_splice_diff_CDSs_count);
		
	}



	## Lengths:
	print "\n";
	print &thousands($sum_gene_lengths)." bp total gene length\n";	
	print &thousands($sum_exon_lengths)." bp total unique exon length\n";	
	print &thousands($sum_intron_lengths)." bp total unique intron length\n";
	print &thousands($sum_5utr_lengths)." bp total 5'UTR length\n";	
	print &thousands($sum_3utr_lengths)." bp total 3'UTR length\n";	
	print &thousands($sum_intergenic_lengths)." bp total intergenic length\n";
	print &thousands(sprintf("%.2f", ($sum_gene_lengths / $gene_count)))." bp average gene length\n";
	print &thousands(sprintf("%.2f", ($sum_exon_lengths / $unique_exon_count)))." bp average exon length\n";
	print &thousands(sprintf("%.2f", ($sum_intron_lengths / $unique_intron_count)))." bp average intron length\n";
	print &thousands(sprintf("%.2f", ($sum_5utr_lengths / $unique_5utr_count)))." bp average 5'UTR length\n" if $unique_5utr_count >0;
	print &thousands(sprintf("%.2f", ($sum_3utr_lengths / $unique_3utr_count)))." bp average 3'UTR length\n" if $unique_3utr_count >0;
	print &thousands(sprintf("%.2f", ($sum_intergenic_lengths / $num_intergenic_regions)))." bp average distance between genes\n";
	
}


if ($EXPORT_FLAG) {
	close $exonfh;
	close $intronfh;
	close $genicfh;
	close $intergenicfh;
	close $utrfh;
}



exit(0);


			
####
sub analyze_intergenics {
	my ($contig, $gene_spans_aref) = @_;
	
	my @gene_spans = sort {$a->[0]<=>$b->[0]} @$gene_spans_aref;
	
	my @intergenics;

	my $gene_span = shift @gene_spans;
	my $prev_rend = $gene_span->[1];
	
	while (@gene_spans) {
		$gene_span = shift @gene_spans;
		my ($gene_lend, $gene_rend) = sort {$a<=>$b} @$gene_span;

		my $gene_length = $gene_rend - $gene_lend + 1;
		$sum_gene_lengths += $gene_length;

		print $genicfh "GENE\t$contig\t$gene_lend\t$gene_rend\t$gene_length\n" if $EXPORT_FLAG;
		
		push (@intergenics, [$prev_rend + 1, $gene_lend - 1] );
		
		$prev_rend = $gene_rend;
	}
	
	foreach my $intergenic (@intergenics) {
		my ($lend, $rend) = @$intergenic;
		my $len = $rend - $lend + 1;
		
		if ($len > 0) {
			$sum_intergenic_lengths += $len;
			$num_intergenic_regions++;
			print $intergenicfh "INTERGENIC\t$contig\t$lend\t$rend\t$len\n" if $EXPORT_FLAG;
		}
	}

	return;
}



		
		
sub thousands($){
        my $val = shift;
        return int(0) if !$val;
        $val = sprintf("%.0f", $val);
        1 while $val =~ s/(.*\d)(\d\d\d)/$1,$2/;
        return $val;
}

