#!/usr/bin/env perl

=pod

=head1 NAME

get_seqs_from_ffindex_index.pl

=head1 USAGE

 two files are needed: -id <File of IDs (see perldoc)> -db <hhblits database> [-invert] [-suffix suffix; defaults to a3m]
 
You can create specific subdatasets from a larger hhblits database, such as UniProt. Because it would have been too slow
otherwise, the entire database is stored in memory as a string (i.e. you will need ca 20-23gb of RAM for a3m, i.e. the size of the file). 
NB: The HHM database is much bigger

# e.g. hypothetical proteins:

 grep -a -E '^\W?#' uniprot20_2013_03_a3m_db| tr -d '\000' | grep -i -E 'hypothetical|putative' | sort -u > hypothetical.headers

 get_seqs_from_ffindex_index.pl -id hypothetical.headers -db uniprot20_2013_03_a3m_db
 get_seqs_from_ffindex_index.pl -id hypothetical.headers -db uniprot20_2013_03_hhm_db -suffix hhm
 get_seqs_from_ffindex_index.pl -id hypothetical.headers -db uniprot20_2013_03.cs219 -suffix cs219
 get_seqs_from_ffindex_index.pl -id hypothetical.headers -db uniprot20_2013_03.fsa.db -suffix fsa


# transposable:

 grep -a -E '^\W?#' uniprot20_2013_03_a3m_db | tr -d '\000' | grep -i -E 'RNA-?directed|transposon|transcriptase|transposable|Reverse|pol |pol-?like|gag |env ' | sort -u > retro.headers

 get_seqs_from_ffindex_index.pl -id retro.headers -db uniprot20_2013_03_a3m_db 
 get_seqs_from_ffindex_index.pl -id retro.headers -db uniprot20_2013_03_hhm_db -suffix hhm
 get_seqs_from_ffindex_index.pl -id retro.headers -db uniprot20_2013_03.cs219 -suffix cs219
 get_seqs_from_ffindex_index.pl -id retro.headers -db uniprot20_2013_03.fsa.db -suffix fsa

# everything but the above (time consuming)

 sort -u hypothetical.headers retro.headers > uniprot20_2013_03_just_useful.headers
 
 get_seqs_from_ffindex_index.pl -id uniprot20_2013_03_just_useful.headers -db uniprot20_2013_03_a3m_db -invert
 get_seqs_from_ffindex_index.pl -id uniprot20_2013_03_just_useful.headers -db uniprot20_2013_03_hhm_db -suffix hhm -invert
 get_seqs_from_ffindex_index.pl -id uniprot20_2013_03_just_useful.headers -db uniprot20_2013_03.cs219 -suffix cs219 -invert
 get_seqs_from_ffindex_index.pl -id uniprot20_2013_03_just_useful.headers -db uniprot20_2013_03.fsa.db -suffix fsa -invert

=cut

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Data::Dumper;

my ($id_file,$db,$invert);
my $suffix = 'a3m';

GetOptions(
	'suffix:s' => \$suffix,
	'id:s'     => \$id_file,
	'db:s'     => \$db,
	'invert'   => \$invert
);

pod2usage "No ID file\n" unless $id_file && -s $id_file;
pod2usage "No DB file\n" unless $db && -s $db;
my $db_index = $db.".index";
if (!-s $db_index){
	my $db_index = $db;
	$db_index =~s/data$/index/;
}

my (%ids,%huge_index);
print "Parsing $id_file...\n";
open (IN,$id_file);
while (my $ln=<IN>){
	my @data = split('\|',$ln);
	if ($data[1]){
		$data[1] =~s/\.[\.a-z0-9]+$//;
		next unless $data[1] && $data[1]!~/^\s*$/;
		$ids{$data[1]} = 1;
	}else{
		chomp($ln);
		$ln =~s/\.[\.a-z0-9]+$//;
		next unless $ln && $ln!~/^\s*$/;
		$ids{$ln} = 1;
	}
}	
close IN;
print "Found ".scalar(keys %ids)." entries\n";

if ($suffix eq 'cs219'){
	print "Loading huge file ($db) in memory...\n";
	my $orig_sep = $/;
	$/ = ">";
	open (HUGEFILE,$db);
	while (my $record=<HUGEFILE>){
		chomp($record);
		next if !$record || $record=~/^\s*$/;
		my @lines = split("\n",$record);
		my $id_line = shift(@lines);
		my $id;
		if ($id_line=~/^\w+\|(\w+)\|/){
			$id = $1;
		}elsif ($id_line=~/^(\w+)\.\w+$/){
			$id = $1;
		}elsif($id_line=~/^(\w+)$/){
			$id = $1;
		}else{
			warn "Didn't recognise this id line\n$id_line\n";
		}
		next unless $id;
		my $size;
		foreach my $ln (@lines){
			$size += length($ln);
		}
		unless ($size){
			warn "No data for $id_line";
			next;
		}
	        if ($invert && !$ids{$id}){
			$huge_index{$id}{'data'} = '>'.$record;
			$huge_index{$id}{'size'} = $size;
	        }elsif(!$invert && $ids{$id}){
			$huge_index{$id}{'data'} = '>'.$record;
			$huge_index{$id}{'size'} = $size;
	        }
	}
	close HUGEFILE;
	$/ = $orig_sep;
	my $total_ids = scalar(keys %huge_index);
	die "No data found!\n" if !$total_ids;
	print "Creating $id_file.$suffix...\n";
	my ($numcsfiles,$num_chars,$counter);
	my $outfile = $invert ? "$id_file.invert.$suffix" : "$id_file.$suffix";
	open (OUT,">$outfile");
	foreach my $id (keys %huge_index){
		$counter++;
		my $data = $huge_index{$id}{'data'};
		next unless $data;
		print OUT $data;
		print "$counter / $total_ids          \r" if $counter % 10000 == 0;
		$numcsfiles++;
		$num_chars += $huge_index{$id}{'size'};
	}
	print "$counter / $total_ids          \n";
	close OUT;
	open (OUT, ">$outfile.sizes");
	print OUT "$numcsfiles $num_chars\n";
	close OUT;
	
	print "Done\n";
	exit(0);
}


print "Loading index ($db_index) of ".scalar(keys %ids)." relevant IDs (or the inverse if -invert was given)...\n";
open (HUGEFILEIDX,$db_index);
while (my $ln=<HUGEFILEIDX>){
	chomp($ln);
	my @data = split("\t",$ln);
	my $id = $data[0];
	$id =~s/\.[\.a-z0-9]+$//;
	next unless $id && $data[2];
	if ($invert){
		@{$huge_index{$id}} = @data[1..2] if !$ids{$id};
	}else{
		@{$huge_index{$id}} = @data[1..2] if $ids{$id};
	}
}
close HUGEFILEIDX;
if (!$invert){
	foreach my $id (keys %ids){
		warn "Not found: $id\n" if !$huge_index{$id};
	}
}	

%ids = ();
my $total_ids = scalar(keys %huge_index);
print "Will create $id_file.$suffix with $total_ids IDs...\n";

print "Loading huge file ($db) in memory...\n";
my $hugestring = '';
open (HUGEFILE,$db);
while (my $ln=<HUGEFILE>){
	$hugestring .= $ln;
}
close HUGEFILE;
my $huge_length = length($hugestring);
print "Loaded $huge_length bytes.\n";

my $counter = int(0);
my $idx_c = int(0);
my $outfile = $invert ? "$id_file.invert.$suffix" : "$id_file.$suffix";
open (OUT,">$outfile");
open (OUTIDX,">$outfile.index");

foreach my $id (keys %huge_index){
	$counter++;
	my $index_data = $huge_index{$id};
	next unless $index_data && $index_data->[1];
	if ($huge_length < $index_data->[0] || $huge_length < ($index_data->[0] + $index_data->[1])){
		warn "ID $id has co-ordinates outside the score of the database file!\n";
		next; # die!
	}
	my $a3m  = substr($hugestring,$index_data->[0],$index_data->[1]);
	warn "Could not get $id if the error is 'substr outside of string' then this is an issue with your (old) version of perl itself. Nothing we can do...\n" unless $a3m;
	next unless $a3m;
	my $entry_length = length($a3m);
	print OUT $a3m;
	print OUTIDX "$id\t$idx_c\t$entry_length\n";
	$idx_c+= $entry_length;
	# in memory, it is really fast...
	print "$counter / $total_ids          \r" if $counter % 10000 == 0;
}
print "$counter / $total_ids             \n";
close OUT;
close OUTIDX;
print "Emptying memory...\n";
$hugestring = '';
%huge_index = ();
print "Done. Sorting...\n";

system("ffindex_build -as $id_file.$suffix $id_file.$suffix.index") if -s "$id_file.$suffix" && -s "$id_file.$suffix.index";

