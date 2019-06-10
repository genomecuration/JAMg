#!/usr/bin/env perl 

sub usage() {

 my @def = @_;

    die '
Usage:	-tabfile             	Input file
 	-score_cutoff        	Defaults to '.$def[0].'
	-align_size_cutoff	Defaults to '.$def[1].'
	-scaff_size_cutoff	Defaults to '.$def[2].'
	-reverse_query_target   If your map was the target, not the query
	-genome_fasta
';
}

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

my $score_cutoff = 1000;
my $align_size_cutoff = 3000;
my $scaff_size_cutoff = 5000;
my ($tabfile,$reverse_query_target,$genome_fasta);

&usage($score_cutoff,$align_size_cutoff,$scaff_size_cutoff) unless GetOptions(
	'score_cutoff:i' => \$score_cutoff,
	'align_size_cutoff:i' => \$align_size_cutoff,
	'scaff_size_cutoff:i'  => \$scaff_size_cutoff,
	'tabfile=s' =>  \$tabfile,
	'reverse_query_target' => \$reverse_query_target,
	'genome_fasta:s' => \$genome_fasta
);

$tabfile = shift if !$tabfile;
&usage($score_cutoff,$align_size_cutoff,$scaff_size_cutoff) unless $tabfile && -s $tabfile;

my (%hash,%size_hash);

open (TAB,$tabfile);

my $query_column = $reverse_query_target ? 6 : 1;
my $target_column = $reverse_query_target ? 1 : 6;

while (my $ln=<TAB>){
	next if $ln=~/^#/;
	chomp($ln);
	next unless $ln;
	my @data = split("\t",$ln);
	next unless $data[10];
	$size_hash{$data[$query_column]} = $data[5] if !$size_hash{$data[$query_column]};

	next if $data[5] < $scaff_size_cutoff || $data[10] < $scaff_size_cutoff;
	my $average_alignment_size = ( $data[3] + $data[8]) / 2;

	if ($data[0] >= $score_cutoff && $average_alignment_size >= $align_size_cutoff){
		my $chr = $data[$target_column];
		$chr =~s/:.+//;
		$hash{$data[$query_column]}{$chr}++;
	} 
}

close TAB;
my ($viol_counter,$support_counter,$viol_size,$support_size,$total_size) = (int(0),int(0),int(0),int(0),int(0));
my %scaffolds_found;

open (OUT1,">$tabfile.$score_cutoff.$align_size_cutoff.$scaff_size_cutoff.out");
foreach my $scaffold (sort {$size_hash{$b} <=> $size_hash{$a} } keys %size_hash){
	next unless $hash{$scaffold};
	my @chrs = sort keys %{$hash{$scaffold}};
	my $print = $scaffold."\t".$size_hash{$scaffold}."\t".scalar @chrs;
	foreach my $chr (@chrs){
		$print .= "\t".$chr.'('.$hash{$scaffold}{$chr}.')';
	}
	print OUT1 $print."\n";
	$total_size +=$size_hash{$scaffold};
	if ($size_hash{$scaffold} >= 1e5){
		if (scalar @chrs >1){
			$viol_counter++;
			$viol_size+=$size_hash{$scaffold};
			$scaffolds_found{$scaffold}++;
		}elsif(scalar @chrs == 1){
			$support_counter++;
			$support_size+=$size_hash{$scaffold};
			$scaffolds_found{$scaffold}++;
		}
	}
}

my $found_size = &thousands($viol_size + $support_size);
$total_size   = &thousands($total_size);
$viol_size    = &thousands($viol_size);
$support_size = &thousands($support_size);
close OUT1;

open (OUT2,">$tabfile.$score_cutoff.$align_size_cutoff.$scaff_size_cutoff.dump");
print OUT2 "#INFO: Searched ".scalar(keys %hash)." scaffolds of a total size $total_size b.p. (big & small). Found "
		.scalar(keys %scaffolds_found)." large (>100kb) scaffolds totaling $support_size bp\n".
	   "#BAD:  Large scaffolds with multiple chromosomes: $viol_counter occupying total of $viol_size bp\n".
           "#GOOD: Large scaffolds with a single chromosome:  $support_counter occupying total of $support_size bp\n";
print OUT2 Dumper \%hash;
close OUT2;

if ($genome_fasta && -s $genome_fasta){
	my $sizes_hash_ref = &process_fasta($genome_fasta);
	open (OUT3,">$tabfile.$score_cutoff.$align_size_cutoff.$scaff_size_cutoff.undetected");
	foreach my $id (keys %{$sizes_hash_ref}){
		next if $scaffolds_found{$id};
		print OUT3 "$id\t".$sizes_hash_ref->{$id}->{'length'}."\t".$sizes_hash_ref->{$id}->{'gaps'}
		 ."\t".sprintf("%.2f",$sizes_hash_ref->{$id}->{'gaps'}/$sizes_hash_ref->{$id}->{'length'})."\n";
	}
	close OUT3;
}

################################################################
sub thousands(){
        my $val = shift;
        $val = sprintf("%.0f", $val);
        return $val if length($val)<4;
        1 while $val =~ s/(.*\d)(\d\d\d)/$1,$2/;
        return $val;
}

sub process_fasta($){
	my $file=shift;
	my $orig_sep = $/;
	$/=">";
	my %genome_sizes_hash;
	open (IN,$file);
	while (my $rec=<IN>){
		chomp($rec);
		next unless $rec;
		my @lines=split("\n", $rec);
		my $id=shift(@lines);
		if ($id=~/^(\S+)\s/){$id=$1;}
		my $seq = join('',@lines);
		$seq=~s/\s+//g;
		next unless $id && $seq;
		$genome_sizes_hash{$id}{'length'} = length($seq);
		$genome_sizes_hash{$id}{'gaps'} = ( $seq =~ tr/N// );
	}
	close IN;
	$/=$orig_sep;
	return(\%genome_sizes_hash);
}
