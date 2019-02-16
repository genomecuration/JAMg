#!/usr/bin/env perl

=pod

=head1 NAME

=head1 USAGE

	-in input file in FASTA/Q or posmap length file (e.g. posmap.scflen)
	-genome genome size in bp for estimating N lengths and indexes
	-single FASTA/Q has sequence in a single line (faster)
	-overwrite => Force overwrite
	-reads	=> Force processing as read data (no N50 statistics)
	-noreads => Force as not being read data. Good for cDNA assemblies with short contigs

=head1 AUTHORS

 Alexie Papanicolaou 1
	
	Ecosystem Sciences, CSIRO, Black Mountain Labs, Clunies Ross Str, Canberra, Australia
	alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

This software is released under the GNU General Public License version 3 (GPLv3).
It is provided "as is" without warranty of any kind.
You can find the terms and conditions at http://www.opensource.org/licenses/gpl-3.0.html.
Please note that incorporating the whole software or parts of its code in proprietary software
is prohibited under the current license.

=head1 BUGS & LIMITATIONS

None known so far.

=cut

use strict;
use warnings;
use Data::Dumper;
use Time::Progress;
use Getopt::Long;
use Pod::Usage;
use Statistics::Descriptive;
$|=1;

my (@infiles,$user_genome_size,$is_fasta,$is_fastq,$is_single,$overwrite,$is_reads,$isnot_reads);
GetOptions(
	'in=s{,}'    => \@infiles,
	'single' =>\$is_single,
	'genome:s' => \$user_genome_size,
	'overwrite' => \$overwrite,
	'reads'	=>\$is_reads,
	'noreads'	=>\$isnot_reads,
);
if (!@infiles){
	@infiles = @ARGV;
}
pod2usage "No input files!\n" if !@infiles;
die "Cannot ask for both reads and noreads options at the same time!\n" if $is_reads && $isnot_reads;
 
if ($is_reads && !$isnot_reads){
	print "Processing all data as reads\n";
}
$user_genome_size=~s/,//g if $user_genome_size;
if ($user_genome_size && $user_genome_size=~/\D$/){
	if ($user_genome_size=~/^(\d+)k/i){
		$user_genome_size=int($1.'000');
	}
	elsif ($user_genome_size=~/^(\d+)m/i){
		$user_genome_size=int($1.'000000');
	}
	elsif ($user_genome_size=~/^(\d+)g/i){
		$user_genome_size=int($1.'000000000');
	}
	print "Genome set to ".&thousands($user_genome_size)." b.p.\n";
}

foreach my $infile (@infiles){
	unless ($infile && -s $infile){warn("I need a posmap length file, e.g. .posmap.scflen for scaffolds\n");pod2usage;}
	my $outfile=$infile.'.n50';
	$outfile.='g' if ($user_genome_size);
	warn ("Outfile $outfile already exists\n") if -s $outfile && !$overwrite;
	next  if -s $outfile && !$overwrite;
	my $total=int(0);
	my ($seq_ref,$gap_distrib_ref,$mask_distrib_ref);
	my @head=`head $infile`;
	foreach (@head){
		if ($_=~/^>\S/){
			$is_fasta=1;
			print "FASTA file found!\n";
			last;
		}elsif($_=~/^@\S/){
			$is_fastq=1;
			print "FASTQ file found!\n";
			last;
		}
	}

	print "Parsing file $infile...\n";
	my $gaps;
	if ($is_fasta){
		($total,$seq_ref,$gap_distrib_ref,$mask_distrib_ref) = &process_fasta($infile);
	}
	elsif($is_fastq){
		($total,$gaps,$seq_ref) = &process_fastq($infile);
	}
	else {
		($total,$gaps,$seq_ref) = &process_csv($infile);
	}
	print "Preparing stats...\n";
	my ($mean,$n90,$n50,$n10,$n25,$n90_length,$n50_length,$n10_length,$n25_length,$scaffolds,$scaffolds_size,$smallest,$largest,$sequence_number,$sum,$genome_size) = &process_stats($seq_ref,$total);

	if ($mean){
		open (OUT,">".$outfile);
		my $stat = Statistics::Descriptive::Full->new();
		$stat->add_data($seq_ref);
		#my $skew='';sprintf("%.2f",$stat->skewness());
		my $mean = sprintf("%.2f",$mean);
		my $median =  sprintf("%.1f",$stat->median());
		my $var = sprintf("%.2f",$stat->variance());
		my $sd = sprintf("%.2f",$stat->standard_deviation());
		my $gap_mean = int(0);my $gap_median =int(0);my $gap_number=int(0);
		if ($gap_distrib_ref){
			my $gap_stat = Statistics::Descriptive::Full->new();
			$gap_stat ->add_data($gap_distrib_ref);
			$gaps = $gap_stat->sum();
			if ($gaps){
				$gap_number = $gap_stat->count();
				$gap_mean = sprintf("%.2f",$gap_stat->mean());
				$gap_median = sprintf("%.1f",$gap_stat->median());
			}
		}
		my $masks = int(0); my $mask_mean = int(0);my $mask_median =int(0);my $mask_number=int(0);
		if ($mask_distrib_ref){
			my $mask_stat = Statistics::Descriptive::Full->new();
			$mask_stat ->add_data($mask_distrib_ref);
			$masks = $mask_stat->sum();
			if ($masks){
				$mask_number = $mask_stat->count();
				$mask_mean = sprintf("%.2f",$mask_stat->mean());
				$mask_median = sprintf("%.1f",$mask_stat->median());
			}
		}

		if (!$scaffolds || $scaffolds == 0){
			$scaffolds=$sequence_number;
			$scaffolds_size=$smallest;
		}
		print OUT "File: $infile\n";
		print OUT "TOTAL: ".&thousands($total)." bp in ".&thousands($sequence_number)." sequences\n";
		print OUT "\tof total ".&thousands($gaps)." are Ns/gaps (occurences: $gap_number; mean: $gap_mean; median: $gap_median).\n";
		print OUT "\tof total ".&thousands($masks)." are masked (lc /X) (occurences: $mask_number; mean: $mask_mean; median: $mask_median).\n";
		print OUT "Mean: ".&thousands($mean)."\nStdev: ".&thousands($sd)."\n";
		print OUT "Median: ".&thousands($median)."\n";
		print OUT "Smallest: ".&thousands($smallest)."\nLargest: ".&thousands($largest)."\n";
		if (($mean >=1000 && !$is_reads) || $isnot_reads){
			print OUT "Assuming genome size of ".&thousands($user_genome_size)." for N50 Stats:\n" if $user_genome_size;
			print OUT "N10 length: ".&thousands($n10_length)."\nN10 Number: ".&thousands($n10)."\n";
			print OUT "N25 length: ".&thousands($n25_length)."\nN25 Number: ".&thousands($n25)."\n";
			print OUT "N50 length: ".&thousands($n50_length)."\nN50 Number: ".&thousands($n50)."\n";
			print OUT "N90 length: ".&thousands($n90_length)."\nN90 Number: ".&thousands($n90)."\n";
			print OUT "Assuming a genome size of "
				.&thousands($user_genome_size)
				." then the top ".&thousands($scaffolds)
				." account for it (min "
				.&thousands($scaffolds_size)
				." bp)\n" if $user_genome_size;
		}else{
			print OUT "Reads found! Read coverage estimated to ".sprintf("%.2f",$total/$user_genome_size)."x using user provided genome size of ".&thousands($user_genome_size)."\n" if $user_genome_size;
		}
		close (OUT);
		print "Done, see $outfile\n";
		system("cat $outfile");
	}else {
		open (OUT,">".$outfile);
		print OUT "File: $infile\n";
		print OUT "TOTAL: $total bp in $sequence_number sequences\n";
		close (OUT);
		warn "Non fatal warning: Something went wrong in estimating the statistics. Maybe the provided genome length is much larger than sequence length or maybe less than 3 sequences provided?\n";
	}
}
########################################################################
sub process_fasta(){
	print "Processing as FASTA\n";
	my $infile=shift;
	my @array ;
	my @gap_distrib ;
	my @mask_distrib ;
	my $total=int(0);
	my $timer   = new Time::Progress;
        $timer->attr( min => 0, max => -s $infile );
	my $counter = int(0);

	if ($is_single){
		open (IN,$infile)||die($!);
		while (my $seq_id=<IN>) {
			print $timer->report( "eta: %E min, %40b %p\r", $counter ) if ( $counter =~ /00000$/ );
	                my $seq=<IN>;
	                my $length=length($seq)-1; # newline
			$counter+=length($seq_id)+$length+1;
	                next unless $length;
			while ($seq=~/([nN\-]+)/g){
				push(@gap_distrib,length($1));
			}
			while ($seq=~/([xXatcg]+)/g){
                                push(@mask_distrib,length($1));
                        }
	                push(@array,$length);
        	        $total+=$length;
		}
		close IN;
	}else{
		my $orig_seq = $/;
		$/ = '>';
		open (IN,$infile)||die($!);
		while (my $record=<IN>) {
			$counter+=length($record);
			print $timer->report( "eta: %E min, %40b %p\r", $counter ) if ( $counter =~ /0000$/ );
			chomp($record);
			next if $record=~/^\s*$/;
			my @data = split("\n",$record);
			my $seq_id = shift(@data);
			my $seq = join('',@data);
			next unless $seq;
			$seq=~s/\s+//g;
			my $length=length($seq);
			next unless $length;
			while ($seq=~/([nN\-]+)/g){
				push(@gap_distrib,length($1));
			}
			while ($seq=~/([xXatcg]+)/g){
                                push(@mask_distrib,length($1));
                        }
			push(@array,$length);
	        	$total+=$length;
	        }
		close IN;
		$/ = $orig_seq;
	}
	print "\n";
	die "No data found or wrong format\n" unless $total;
	return ($total,\@array,\@gap_distrib,\@mask_distrib);
}
sub process_fastq(){
	print "Processing as FASTQ\n";
	my $infile=shift;
	my @array ;
	my $total=int(0);
	my $gaps=int(0);
	my $timer   = new Time::Progress;
        $timer->attr( min => 0, max => -s $infile );
	my $counter = int(0);
	open (IN,$infile);
	while (my $seq_id=<IN>) {
		print $timer->report( "eta: %E min, %40b %p\r", $counter ) if ( $counter =~ /00000$/ );
		my $seq=<IN>;
		my $scrap=<IN>.<IN>;
		my $length=length($seq)-1; #newline
		$counter+=length($seq_id)+$length+1;
		next unless $length;
		$gaps+=($seq=~tr/[xXNn\-]//);
                push(@array,$length);
                $total+=$length;
       }
	close IN;
	print "\n";
	die "No data found or wrong format\n" unless $total;
	return ($total,$gaps,\@array);
}
sub process_csv(){
	print "Processing as CSV\n";
	my $infile=shift;
	my @array ;
	my $total=int(0);
	my $gaps='N/A';
	my $timer   = new Time::Progress;
        $timer->attr( min => 0, max => -s $infile );
	my $counter = int(0);
	open (IN,$infile)||die ("Cannot open $infile\n");
	while (my $ln=<IN>){
		print $timer->report( "eta: %E min, %40b %p\r", $counter ) if ( $counter =~ /00000$/ );
		$counter+=length($ln);
		$ln=~/(\d+)$/;
		next unless $1;
		my $length= $1;
		push(@array,$1);
		$total+=$length;
	}
	close IN;
	die "No data found or wrong format\n" unless $total;
	return ($total,$gaps,\@array);
}

sub process_stats(){
	my $sequences_ref = shift;
	my $total = shift;
	my $genome_size = int(0);
	my $mean = $total / scalar(@$sequences_ref);
	my ($n90,$n50,$n10,$n25,$n90_length,$n50_length,$n10_length,$n25_length,$scaffolds,$scaffolds_size,$smallest,$largest,$sequence_number,$sum);
	print "Sorting...";
	my @sequences=sort{$b<=>$a} @$sequences_ref;
	$smallest=$sequences[-1];
	$largest=$sequences[0];
	print " done!\n";
	$|=0;
	if (($mean < 1000 && !$isnot_reads) || $is_reads ){
		print "Reads detected. Ignoring N* calculations.\n";
		$genome_size=$total;
		$sequence_number = scalar(@$sequences_ref);
         	return ($mean,$n90,$n50,$n10,$n25,$n90_length,$n50_length,$n10_length,$n25_length,$scaffolds,$scaffolds_size,$smallest,$largest,$sequence_number,$sum,$genome_size);
	}
	elsif (!$user_genome_size){
		print "Setting genome size for N* calculations to total consensus $total\n";
		$genome_size=$total;
	}elsif($user_genome_size){
		print "Setting genome size for N* calculations to user defined $user_genome_size\n";
		$genome_size = $user_genome_size;
	}
	
	foreach my $sequence_length ( @sequences){
		$sum+=$sequence_length;
		$sequence_number++;
		if($sum >= $genome_size*0.1 && !$n10){
			$n10=$sequence_number;
			$n10_length=$sequence_length;
		}
		elsif($sum >= $genome_size*0.25 && !$n25){
			$n25=$sequence_number;
			$n25_length=$sequence_length;
		}
		elsif($sum >= $genome_size*0.5 && !$n50){
			$n50 = $sequence_number;
			$n50_length=$sequence_length;
		}
		elsif($sum >= $genome_size*0.90 && !$n90){
			$n90=$sequence_number;
			$n90_length=$sequence_length;
		}
                elsif ($sum >= $genome_size && !$scaffolds){
			$scaffolds = $sequence_number;
			$scaffolds_size = $sequence_length;
		}

	}
	print "Processed $sequence_number sequences\n";
	return ($mean,$n90,$n50,$n10,$n25,$n90_length,$n50_length,$n10_length,$n25_length,$scaffolds,$scaffolds_size,$smallest,$largest,$sequence_number,$sum,$genome_size);
}

sub thousands($){
	my $val = shift;
	return int(0) if !$val;
	$val = sprintf("%.0f", $val);
	1 while $val =~ s/(.*\d)(\d\d\d)/$1,$2/;
	return $val;
}
