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
#	for (my $i=0;$i<@gene_data;$i++){
#		$gene_data[$i] = &pre_process_webapollo($gene_data[$i]);
#	}

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
	die "FATAL: Gene not unique: $gene_id (".$unique_id{$gene_id}.")\n".$gene_data[0]."\n" if $unique_id{$gene_id};

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


##########
sub pre_process_webapollo(){
	my (%seen_names,$death);
	my ($last_gene_name,$last_gene_id,$last_transcript_name);
	my $line = shift;
	next if ($line=~/^\s*$/ || $line=~/^#/);
	my @data=split("\t",$line);
	next unless $data[8];
			next if ($data[3] <= 0 || $data[4] <= 0);
			next if $data[2] eq 'non_canonical_five_prime_splice_site';
	        	$data[1] = 'WebApollo' if $data[1] eq '.';
	        	$data[2] = 'mRNA' if ($data[2] && $data[2] eq 'transcript');
        	        $data[8]=~s/Name=ID=/ID=/;
			if ($data[8]=~/Name=([^;]+)/){
				my $name = $1;
				if ($name=~s/\s+/_/g){
					$data[8]=~s/Name=[^;]+/Name=$name/;
				}
			}
		

			chomp($data[8]);
			if ($data[2] eq 'mRNA' || $data[2] eq 'gene'){
        	           if ($data[8]=~/Name=([^;]+)/){
				$data[8].= ';Alias='.$1 if $do_alias;
				$last_gene_name = $1 if $data[2] eq 'gene';
	                   }elsif ($data[8]=~/ID=([^;]+)/){
				$data[8].= ';Alias='.$1 if $do_alias;
				$last_gene_name = '' if $data[2] eq 'gene';
	                   }else{
				$last_gene_name = '' if $data[2] eq 'gene';
			   }
			}
			if ($data[2] eq 'gene'){
				$death.= "This user provided name already exists: $last_gene_name\n" if ($seen_names{$last_gene_name});
				$seen_names{$last_gene_name}++;
				$data[8] =~s/ID=[^;]+/ID=$last_gene_name/;
			}elsif ($data[2] eq 'mRNA'){
				$data[8] =~s/Parent=[^;]+/Parent=$last_gene_name/;
          	                $data[8]=~/Name=([^;]+)/;
                                $last_transcript_name=$1;
				#$death.= "This user provided name already exists: $last_transcript_name\n" if ($seen_names{$last_transcript_name});
				while ($seen_names{$last_transcript_name}){
					$last_transcript_name.='_';
	          	                $data[8]=~s/Name=[^;]+/Name=$last_transcript_name/;
				}
				$seen_names{$last_transcript_name}++;
				$data[8] =~s/ID=[^;]+/ID=$last_transcript_name/;				
			}else{
				$data[8] =~s/Parent=[^;]+/Parent=$last_transcript_name/;


			}
	
#wa bug fixed			$data[4]--; 			# && ($data[2] eq 'mRNA' || $data[2] eq 'gene' ||$data[2] eq 'exon' || $data[2] eq 'CDS'));
			$data[8] = encode_entities($data[8]);
	die "ERROR:\n".$death if $death;
	return @data;
}
