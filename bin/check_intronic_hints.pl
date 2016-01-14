#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util qw(sum);
use Pod::Usage;
use File::Basename;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";
use Fasta_reader;

	my $hint_file = shift;
	my $genome = shift ||die;
	print "Reading Genome $genome\n";
	my $fasta_obj = new Fasta_reader($genome);
	my %fasta_data = $fasta_obj->retrieve_all_seqs_hash();
	print "Parsing HINT file $hint_file\n";
	# use the introns to set orientation for both introns and exons
	# (this would require parsing the file a couple of time)
	open (IN,$hint_file) || die;
	# # spit out splice sites while at it
	open (SPLICE1,">$hint_file.splice.augustus");
	open (SPLICE2,">$hint_file.splice.gmap");
	open (GFF,">$hint_file.intron.only");
	my $printed_counter=int(1);
	my (%intron_count,%intron_printed);
	while (my $ln=<IN>){
		my @data = split("\t",$ln);
		next unless $data[6] && $data[2] eq 'intron';
		my ($ref,$type,$start,$end,$strand) = ($data[0],$data[2],$data[3],$data[4],$data[6]);
		die "Cannot find sequence for $ref in genome FASTA file\n" unless $fasta_data{$ref};
		die "Unexpected GFF with start higher than the end\n$ln" if $start > $end;
		chomp($data[8]);
		# acceptor and donor sites (3' and 5' on + strand)
		my $site1 = uc(substr($fasta_data{$ref},($start-1),2));
		my $site2 = uc(substr($fasta_data{$ref},($end-1-1),2));
		my ($donor,$acceptor,$type1,$type2);
		if (
		(($site1 eq 'GT' || $site1 eq 'GC') && $site2 eq 'AG') ||
		($site1 eq 'AT' && $site2 eq 'AC')
		 ){
			$strand = '+';
			($donor,$acceptor,$type1,$type2) = ($site1,$site2,'dss','ass');
		}elsif (
		(($site2 eq 'AC' || $site2 eq 'GC') && $site1 eq 'CT') ||
		($site2 eq 'AT' && $site1 eq 'GT')
		 ){
			$strand = '-';
			($donor,$acceptor,$type1,$type2) = ($site2,$site1,'ass','dss');
		}else{
			# leave them unchanged
			$data[8] .='noncanonical=true;';
			print GFF join("\t",@data)."\n";
			next;
		}
		$intron_count{"$ref:$start..$end"}++;

		#only if supported by at least 30 introns, then save in database
		if (!$intron_printed{"$ref:$start..$end"} && $intron_count{"$ref:$start..$end"} && $intron_count{"$ref:$start..$end"} >= 30){
			# sequence 40 bp up and 40 bp upstream
			my $site1_seq = lc(substr($fasta_data{$ref},($start-1-40),82));
			my $site2_seq = lc(substr($fasta_data{$ref},($end-1-1-40),82));
			$site1_seq =~ tr/n/-/;$site2_seq =~ tr/n/-/;
			if ($strand eq '-'){
				$site1_seq = &revcomp($site1_seq);
				$site2_seq = &revcomp($site2_seq);
			}
			#augustus
			print SPLICE1 "$type1 $site1_seq\n" if $site1_seq;
			print SPLICE1 "$type2 $site2_seq\n" if $site2_seq;
			#gmap
			my $intron_gsnap_txt = ($strand eq '-') ? ">intron.$printed_counter $ref:$end..$start\n" : ">intron.$printed_counter $ref:$start..$end\n";
			print SPLICE2 $intron_gsnap_txt;
			$intron_printed{"$ref:$start..$end"}++;
			$printed_counter++;
		}

		$data[6] = $strand;
		$data[8].="donor=$donor;acceptor=$acceptor;noncanonical=false;";
		print GFF join("\t",@data)."\n";
	}
	close IN;
	close GFF;
	close SPLICE1;
	close SPLICE2;

