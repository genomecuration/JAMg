#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
my $file = shift||die;

open (IN,$file);
open (OUT,">$file.gtf");

my @gene_data;

while (my $ln=<IN>){
	next if $ln=~/^#/ || $ln!~/AUGUSTUS/;
	chomp($ln);
	my @data = split("\t",$ln);
	if ($data[2] eq "5'-UTR"){
		$data[2] = '5UTR';
	}elsif ($data[2] eq "3'-UTR"){
		$data[2] = '3UTR';
	}elsif($data[2] eq 'tss'){
next;
	}elsif($data[2] eq 'tts'){
next;
	}elsif($data[2] eq 'intron'){
next;
	}elsif($data[2] eq 'gene'){
  	 if (@gene_data){
		my @first_gene_data = split("\t",$gene_data[0]);
		my $strand = $first_gene_data[6];
		if ($strand eq '-'){
			for (my $g=0;$g<(@gene_data);$g++){
				my @data = split("\t",$gene_data[$g]);
				if ($data[2] eq 'CDS'){
					$data[3]--; # remove stop codon bug
					$gene_data[$g] = join("\t",@data);
					last;
				}
			}
		}else{
			for (my $g=(scalar(@gene_data)-1);$g>=0;$g--){
				my @data = split("\t",$gene_data[$g]);
				if ($data[2] eq 'CDS'){
					$data[4]-=3; # remove stop codon
					if ($data[4] > $data[3]){
						#sometimes augustus prints a 3 base CDS (which is the stop codon)
						$gene_data[$g] = join("\t",@data);
					}else{
						$gene_data[$g] = '';
					}
					last;
				}
			}
		}
		print OUT join("",@gene_data)."\n";
		undef(@gene_data);
	}
next;
		my $gene_id = $data[8];
		$data[8] = 'gene_id "'.$gene_id.'";';
		my $res = join("\t",@data)."\n";
		push(@gene_data,$res);
	}elsif($data[2] eq 'transcript'){
		my $tr_id = $data[8];
		my $gene_id = $tr_id;
		$gene_id=~s/\.\d+$//;
		$data[8] = 'transcript_id "'.$tr_id.'"; gene_id "'.$gene_id.'";';
		my $res = join("\t",@data)."\n";
		push(@gene_data,$res);
next;
	}
	$data[8]=~/gene_id\s(\S+)/;
	my $newdata ="gene_id $1";
	$data[8]=~/transcript_id\s(\S+)/;
	$newdata.=" transcript_id $1";
        $newdata.=" Name $1";
	$data[8] = $newdata;
	my $res = join("\t",@data)."\n";
	push(@gene_data,$res);
}
close IN;
	if (@gene_data){
		print OUT join("",@gene_data);
		undef(@gene_data);
	}
close OUT;
