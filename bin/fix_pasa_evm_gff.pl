#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use HTML::Entities;

# reads the whole thing into memory.

my $infile = shift||die;
my $do_alias = shift;


open (IN,$infile)||die;
my (@file_lines,%states);
while (my $ln = <IN>){
	last if $ln=~/^##FASTA/;
	next if $ln=~/^##[^#]/;
	if ($ln=~/status=([^;]+)/){$states{$1}++;}
	push (@file_lines,$ln);
}
close IN;


$states{'undefined'}=1;
$states{'weird'}=1;
my %file_handles;
foreach my $state (keys %states){
	my $outfile = $infile .".$state.gff3";
	open (my $fh,">$outfile");
	print $fh "##gff-version 3\n";
	$file_handles{$state} = $fh;
}
my (%unique_id,%unique_name);

&process_webapollo_gff(\@file_lines);

foreach my $state (keys %file_handles){
	my $fh = $file_handles{$state};
	close $fh;
	my $outfile = $infile .".$state.gff3";
	unlink($outfile) if -s $outfile < 20;
}


##########################################################3
sub process_webapollo_gff(){
	my $ref = shift;
	my (@gene_data); 

	foreach my $ln (@$ref){
		next if $ln=~/^#\s/; #PASA comments
		last if $ln=~/^##FASTA/;
	        if ($ln=~/^###/ || $ln=~/^\s*$/){
			&process_gene_data(@gene_data);
			@gene_data=();
		}else{
			push(@gene_data,$ln);
		}
	}
	&process_gene_data(@gene_data);
	@gene_data=();
}



sub process_gene_data(){
	my @gene_data = @_;

	return if (!$gene_data[0] || !$gene_data[1]);

	my ($gene_id,$gene_name);

	#checks
	if ($gene_data[0] && $gene_data[0] !~/\bgene\b/){
#		warn "Unexpected first line type. Saving it in 'weird' file\n".join('',$gene_data[0]);
		my $fh = $file_handles{'weird'};
		print $fh join('',(@gene_data,"###\n"));	
	}
	if (!$gene_data[1]){
		warn "No data for this gene. Skipping\n".join('',@gene_data);
		return;
	}
	# more checks
	if ($gene_data[0]=~/ID=([^;]+)/){$gene_id = $1;}
	die "FATAL: No Gene ID" unless $gene_id;

	if ($unique_id{$gene_id}){
		warn "DANGEROUS: SKIPPING gene not unique: $gene_id (".$unique_id{$gene_id}.")\n".$gene_data[0]."\n";
		return;
	}

	if ($gene_data[0]=~/Name=([^;]+)/){$gene_name = $1;}
	unless ($gene_name){ 
		warn "Cannot find a name for a gene. Skipping\n".join('',@gene_data);
		return;
	}

	# need to consider each transcript separately.
	for (my $i=1;$i<@gene_data;$i++){
		my @data = split("\t",$gene_data[$i]);
		
		next if $data[2] ne 'mRNA';
		my ($transcript_name,$transcript_id,$old_transcript_id);
		if ($data[8]=~/ID=([^;]+)/){
			$transcript_id = $1;
		}
		die "FATAL: No Transcript ID" unless $transcript_id;
			if ($data[8]=~/Name=([^;]+)/){
			my $transcript_name = $1;
		}
		my $state = "undefined";
		if ($data[8]=~/status=([^;]+)/){
			$state = $1;
		}
		if (!$transcript_name){
			$transcript_name = $transcript_id;
		}
		if ($unique_id{$transcript_id}){
			$old_transcript_id = $transcript_id;
			while (	$unique_id{$transcript_id}){
				$transcript_id.='x';
			};
			warn "WARNING: Transcript not unique: $old_transcript_id \n".$gene_data[0].$gene_data[$i]." Renamed to $transcript_id\n";
		}
		$unique_id{$transcript_id}++;

		# if this is a curated transcript_name and the gene name is wrong, edit it.
		if ($transcript_name=~/^(\S+)[\.\-_]\d+$/ && $transcript_name!~/^scaffold/){
			my $new_gene_name = $1;
			if ($gene_name=~/^\w{32}-\d+$/){
				$gene_data[0]=~s/Name=[^;]+/Name=$new_gene_name/;
			}else{
			   $new_gene_name = $gene_name;	
			}
		}
		# now get all the other lines
		my @mrna_data;
		# we store the gene only once ($gene_data[0])
		push (@mrna_data,$gene_data[0]) unless $unique_id{$gene_id};

		if ($old_transcript_id){ # has been renamed
			$gene_data[$i]=~s/ID=$old_transcript_id([;\s])/ID=$transcript_id$1/;
		}

		push (@mrna_data,$gene_data[$i]);

		for (my $k=$i+1;$k<@gene_data;$k++){
			my @d = split("\t",$gene_data[$k]);
			next unless $d[2] eq 'exon' || $d[2] eq 'CDS';
			if ($old_transcript_id && $d[8]=~/=$old_transcript_id([\.;\s])/){ # has been renamed
				$d[8]=~s/=$old_transcript_id([\.;\s])/=$transcript_id$1/g;
				$gene_data[$k] = join("\t",@d);
			}
			push(@mrna_data,$gene_data[$k]) if $d[8]=~/Parent=$transcript_id[;\s]/;
		}
		# now store in GFF
		my $fh = $file_handles{$state};
		print $fh "###\n" unless $unique_id{$gene_id};
		print $fh join('',@mrna_data);
		#gene id goes within transcript so that we get it at least once
		$unique_id{$gene_id}++;
	}
}


