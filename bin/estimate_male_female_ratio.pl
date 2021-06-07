#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use Time::localtime;
use Statistics::Descriptive;
use File::Basename;

# for JAMg
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";

=pod

=head1 NAME

 estimate_male_female_ratio.pl

=head1 SYNOPSIS

 estimate_male_female_ratio.pl -female [coverage file] -male [coverage file] -autosome CONTIGID -window 100

 Mandatory options
	-genome_fasta  s	=> Genome fasta file 
	-female_depth  s	=> Bigwig file with depth coverage of female data
	-male_depth    s	=> Bigwig file with depth coverage of male data
	-autosome      s	=> The ID of a contig/scaffold that is autosomal
	-window_size   i	=> The approximate size of the window to estimate the median coverage over (defaults to 100)
	-repeat_cutoff i        => Repeat cutoff: ignore coverage above this number (defaults to autocalculate of normalised autosome mean from both datasets * 2)

=cut

## -autosome for drosophila can use 2L: >NT_033779.5 Drosophila melanogaster chromosome 2L
## - Y contig: NW_001846231.1
## - X contig: NW_001845796.1



my $cwd = `pwd`;
chomp($cwd);
my ($debug, $verbose, $help);

my ($female_depth_file, $male_depth_file, $autosome_contig, $genome_fasta ); # mandatory options
my ($bigwigsummary_exec) = &check_program('bigWigSummary');

# optional:
my $window_size = 100;
my $repeat_cutoff;
my $outbasename;

&GetOptions(
	'female_depth:s' => \$female_depth_file,
	'male_depth:s'   => \$male_depth_file,
	'autosome:s'     => \$autosome_contig,
	'window_size:i'  => \$window_size,
	'verbose'	=> \$verbose,
	'debug'		=> \$debug,
	'help'		=> \$help,
	'genome_fasta:s' => \$genome_fasta,
	'repeat_cutoff:i' => \$repeat_cutoff,
	'outfile:s'	=> \$outbasename
);

pod2usage if !$female_depth_file || !$male_depth_file || !$autosome_contig || !$genome_fasta;
&check_files();


my $genome_data_hashref = &get_genome_data();
my $found_norm_factor = &estimate_norm_factor();
my ($autosome_mean,$autosome_median) = &estimate_coverage($autosome_contig);

print "Using $autosome_contig (size="
   . &thousands($genome_data_hashref->{$autosome_contig}->{'length'} )
   . ") as an autosomic reference for male / female coverage (using ".&estimate_bigwig_datapoints($autosome_contig)." bins from an approximate window of $window_size)\n";
#sleep(2);
	
print "Overall Male mean depth of coverage of $autosome_contig is: ". &estimate_whole_unnormalised_contig_mean($autosome_contig, $male_depth_file) ."\n";
print "Overall Female mean depth of coverage of $autosome_contig is: ". &estimate_whole_unnormalised_contig_mean($autosome_contig, $female_depth_file) ."\n";
print "Normalisation factor is ".sprintf("%.2f", $found_norm_factor)."\n";
print "Post normalisation and window-based ratios: Mean is $autosome_mean and Median is $autosome_median\n";
print "\n\n";

open (OUT1,">$outbasename.overall.tsv");
open (OUT2,">$outbasename.window.tsv");

print OUT1 "CONTIG_ID\tCONTIG_SIZE\tMedian\tMean\tCONTIG_DESCRIPTION\n";
print OUT2 "CONTIG_ID\tWINDOW_SIZE\tWindow_Mean\tWindow_number\n";


my $counter = 0;
my $number_of_contigs = scalar(keys %{$genome_data_hashref} ); 
foreach my $contig_id (sort keys %{$genome_data_hashref} ) {
	$counter++;
	my ($mean, $median, $ratio_ref, $true_window_size) =  &estimate_coverage($contig_id);
	print OUT1 "$contig_id\t"
		.$genome_data_hashref->{$contig_id}->{'length'}
		."\t$median\t$mean"
		."\t".$genome_data_hashref->{$contig_id}->{'description'}
		."\n";

#	foreach my $window_ratio (@$ratio_ref){
#		print OUT2 "$contig_id\t$window_ratio\n";
#	}	

	for (my $i=0; $i< scalar(@$ratio_ref); $i++){
		my $window_mean_ratio = $ratio_ref->[$i];
		print OUT2 "$contig_id\t$true_window_size\t$window_mean_ratio\t$i\n";
	}
	print "  Processed $counter/$number_of_contigs      \r" if $counter % 10 == 0;
}
close OUT1;
close OUT2;
print "  Processed $counter/$number_of_contigs      \n\n";


####################
sub estimate_coverage(){
	my $contig_id = shift;

	my $number_of_windows = &estimate_bigwig_datapoints($contig_id);
	my $true_window_size = int( $genome_data_hashref->{$contig_id}->{'length'} / $number_of_windows);

	# bigWigSummary -type=mean male_uniq.coverage.bw NW_001845718.1 0 499 10
	#121.043	328.18	555.98	745.6	973.56	1016.44	882.52	1108.46	1223.56	1184.04

	# get the means for each window for male and female datasets

	my $female_max_cmd = $bigwigsummary_exec 
		. " -type=max " . $female_depth_file . " " . $contig_id . " 0 " 
		. $genome_data_hashref->{$contig_id}->{'length'} . " " 
		. $number_of_windows;

	my $male_max_cmd = $bigwigsummary_exec 
		. " -type=max " . $male_depth_file . " " . $contig_id . " 0 " 
		. $genome_data_hashref->{$contig_id}->{'length'} . " " 
		. $number_of_windows;

	my $female_mean_cmd = $bigwigsummary_exec 
		. " -type=mean " . $female_depth_file . " " . $contig_id . " 0 " 
		. $genome_data_hashref->{$contig_id}->{'length'} . " " 
		. $number_of_windows;

	my $male_mean_cmd = $bigwigsummary_exec 
		. " -type=mean " . $male_depth_file . " " . $contig_id . " 0 " 
		. $genome_data_hashref->{$contig_id}->{'length'} . " " 
		. $number_of_windows;

	my @female_maxs = split("\t",`$female_max_cmd 2> /dev/null`);
	my @male_maxs = split("\t",`$male_max_cmd 2> /dev/null`);
	
	my @female_means = split("\t",`$female_mean_cmd 2>> $outbasename.err`);
	my @male_means = split("\t",`$male_mean_cmd  2>> $outbasename.err`);

	chomp($female_maxs[-1]) if $female_maxs[-1];
	chomp($male_maxs[-1]) if $male_maxs[-1];
	chomp($female_means[-1]) if $female_means[-1];
	chomp($male_means[-1]) if $male_means[-1];


	my @temp_ratios;

	for (my $i=0; $i < $number_of_windows; $i++){
	
		next if (
				( $female_maxs[$i] && $female_maxs[$i] ne 'n/a' && $female_maxs[$i] > $repeat_cutoff )
				|| 
				( $male_maxs[$i] && $male_maxs[$i] ne 'n/a' && $male_maxs[$i] > $repeat_cutoff)
			);

		if (!$male_means[$i] || $male_means[$i] eq 'n/a'){ $male_means[$i] = 1;}
		if (!$female_means[$i] ||   $female_means[$i] eq 'n/a'){
			$female_means[$i] = 1;
			$temp_ratios[$i] = $male_means[$i] / $female_means[$i];
			# equivalent to $ratios[$i] = $male_means[$i];
		}
		else { 
			$temp_ratios[$i] = $male_means[$i] / ( $female_means[$i] * $found_norm_factor );
		}
	}

	my @ratios;
	# clean up from empty values:
	foreach my $v (@temp_ratios){
		next if !$v;
		push(@ratios,$v);
	}


	# print Dumper \@ratios;

	my $stat = Statistics::Descriptive::Full->new();
	$stat->add_data(\@ratios);
	my $median = $stat->median();
	my $mean = $stat->mean();

	$median = $median ? sprintf("%.2f", $median) : 'NA';
	$mean = $mean ? sprintf("%.2f", $mean) : 'NA';


	#print "Mean is: $mean and median is $median\n";
	return ($mean,$median,\@ratios,$true_window_size);

}

sub estimate_norm_factor(){
	my $male_mean = &estimate_whole_unnormalised_contig_mean($autosome_contig,$male_depth_file); 
	my $female_mean = &estimate_whole_unnormalised_contig_mean($autosome_contig,$female_depth_file); 
	my $norm_factor = $male_mean / $female_mean;
	# the calc below is unnecessary but keep for making clear (could just be $male_mean * 2)
	$repeat_cutoff =  int ( ($male_mean + $female_mean * $norm_factor ) / 2) * 2 if !$repeat_cutoff;
	return($norm_factor);
}

sub estimate_whole_unnormalised_contig_mean(){
	my $contig_id = shift;
	my $dataset = shift;

	my $cmd = $bigwigsummary_exec 
		. " -type=mean " . $dataset . " " . $contig_id . " 0 " 
		. $genome_data_hashref->{$contig_id}->{'length'} . " " 
		. 1;

	my $result = `$cmd`;
	chomp($result);
	return $result;
}


sub estimate_bigwig_datapoints(){
	my $contig_id = shift;
	my $datapoints = $genome_data_hashref -> {$contig_id}->{'length'} / $window_size;
	return int( $datapoints);
}

sub get_genome_data(){
	my %hash;

	my $orig_sep = $/;
	$/ = '>';

	open (IN,$genome_fasta);
	while (my $record = <IN>){
		chomp($record);
		next if !$record;
	
		my @lines = split("\n",$record);
		next if !$lines[1];

		my $id_line = shift(@lines);
		my ($id, $description) = ('','');
		if ($id_line=~/^(\S+) (.+)/){
			$id = $1;
			$description = $2;
		}else{
			$id = $id_line;
		}

		my $sequence = uc(join('',@lines));
		my $seq_length = length($sequence);

		$hash{$id}{'length'} = $seq_length;
		$hash{$id}{'description'} = $description if $description;
		$hash{$id}{'description'} = '' if !$description;

		# $hash{$id}{'complexity'} = XXXXX
		# number of gaps

	}

	close IN;

	$/ = $orig_sep;

	return \%hash;
}


sub mytime() {
 my @mabbr =
   qw(January February March April May June July August September October November December);
 my @wabbr = qw(Sunday Monday Tuesday Wednesday Thursday Friday Saturday);
 my $sec   = localtime->sec() < 10 ? '0' . localtime->sec() : localtime->sec();
 my $min   = localtime->min() < 10 ? '0' . localtime->min() : localtime->min();
 my $hour =
   localtime->hour() < 10 ? '0' . localtime->hour() : localtime->hour();
 my $wday = $wabbr[ localtime->wday ];
 my $mday = localtime->mday;
 my $mon  = $mabbr[ localtime->mon ];
 my $year = localtime->year() + 1900;
 return "$wday, $mon $mday, $year: $hour:$min:$sec\t";
}


sub check_program() {
 my @progs = @_;
 my @paths;
 foreach my $prog (@progs) {
  my $path = `which $prog`;
  pod2usage "Error, path to a required program ($prog) cannot be found\n\n"
    unless $path =~ /^\//;
  chomp($path);
  push( @paths, $path );
 }
 return @paths;
}


sub process_cmd {
 my ( $cmd, $dir ) = @_;
 print &mytime . "CMD: $cmd\n" if $debug || $verbose;
 undef($dir) if $dir && $dir eq '.';
 chdir($dir) if $dir;
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
   chdir($cwd) if $dir;
   die "Error, cmd died with ret $ret\n";
 }
 return;
}


sub check_files(){
	&check_file_exists($female_depth_file);
	&check_file_exists($male_depth_file);
	&check_file_exists($genome_fasta);
	if (!$outbasename){
		$outbasename = $genome_fasta;
		$outbasename =~s/^[\.\/]+//;
		$outbasename =~s/\.\S+$//;
	}
	unlink("$outbasename.err");
	unlink("$outbasename.overall.tsv");
	unlink("$outbasename.window.tsv");
}

sub check_file_exists(){
	my $file = shift;
	die "File $file does not exist\n" if !-s $file;
}

sub thousands(){
        my $val = shift;
        $val = sprintf("%.0f", $val);
        return $val if length($val)<4;
        1 while $val =~ s/(.*\d)(\d\d\d)/$1,$2/;
        return $val;
}
