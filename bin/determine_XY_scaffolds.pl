#!/usr/bin/env perl

=head1 NAME

 determine_XY_scaffolds.pl

=head1 USAGE

Mandatory

    -hetero   s => hetero aligned depth data. Tab delimited samtools depth unless -bedgraph
    -homozy   s => homozy aligned depth data. Tab delimited samtools depth unless -bedgraph
    -fasta    s => Genome sequence fasta file

Optional
    
    -bedgraph        => Bedgraph format. Required.
    -checkseq  :s{,} => User provided fasta IDs to check
    -cutgaps   :i    => Maximum proportion of gaps allowed (0-100, def. 50)
    -cutsize   :i    => Minimum size of scaffolds to consider (def 1e5)
    -cutdepth  :i    => Depth cutoff to separate heterozygous from homozygous scaffolds. Autocalculated as half the genome median
    -cutmedian :f    => Multiplier added to autocalculation of -cutdepth median/2 (defaults to 0.20).
    -report_all
    -domean          => trimmed mean (30%)
    -noscale

=cut
 
use strict;
use warnings;
use Data::Dumper;
use Statistics::Descriptive;
use Getopt::Long;
use Pod::Usage;
$|=1;
my $median_multiplier = 0.20;
my $user_max_depth_cutoff;
my $size_cutoff = 1e5;
my $gap_cutoff = 50;
my $add_zeros;
my $use_bedgraph = 1;
my ($genome_fasta,$hetero_align_file,$homozy_align_file,@seqs_to_check,$scaffold_count,$do_all,$do_mean,$no_scale);

GetOptions(
    'hetero=s' => \$hetero_align_file,
    'homozy=s' => \$homozy_align_file,
    'fasta=s' => \$genome_fasta,
    'checkseq:s{,}' => \@seqs_to_check,
    'cutgaps:i' => \$gap_cutoff,
    'cutsize:i' => \$size_cutoff,
    'cutdepth:i' => \$user_max_depth_cutoff,
    'cutmedian:f' => \$median_multiplier,
    'bedgraph' => \$use_bedgraph,
    'report_all|all' => \$do_all,
    'domean|do_mean' => \$do_mean,
    'noscale|no_scale' => \$no_scale,
);
pod2usage "No hetero data\n" unless $hetero_align_file && -s $hetero_align_file;
pod2usage "No homozy data\n" unless $homozy_align_file && -s $homozy_align_file;

my $sizes_hash_ref = &process_fasta($genome_fasta);
my %user_seqs;
if (@seqs_to_check){
	foreach my $seq (@seqs_to_check){
		$user_seqs{$seq}++;
	}
}
$hetero_align_file=&process_depth($hetero_align_file);
$homozy_align_file=&process_depth($homozy_align_file);

&process_xy();

#&process_align($hetero_align_file);
#&process_align($homozy_align_file);
#&process_align($hetero_align_file,$user_max_depth_cutoff);


#################
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
		$scaffold_count++;
	}
	close IN;
	$/=$orig_sep;
	die unless $scaffold_count && $scaffold_count >0;
	print "Found $scaffold_count genome sequences in $file\n";
	return(\%genome_sizes_hash);
}

sub process_xy(){
        my (%hash_hetero,%hash_homozy,%scaffold_size);

	open (HETERO,$hetero_align_file) || die $!;
        while (my $ln=<HETERO>){
                chomp($ln);
		next unless $ln;
                my @data = split(',',$ln);
		next unless $data[1];
                if ($data[1] < 2){
#                        warn $data[0]." : skipped due to 0 coverage\n";
                }else{
                        $hash_hetero{'median'}{$data[0]} = $data[1];
                        $hash_hetero{'coverage'}{$data[0]} = $data[2];
			$scaffold_size{$data[0]} = $data[3];
                        $hash_homozy{'coverage'}{$data[0]} = $hash_homozy{'coverage'}{$data[0]} ? $hash_homozy{'coverage'}{$data[0]} : int(0);
                }
        }
        close (HETERO);
	open (HOMOZY,$homozy_align_file) || die $!;
        while (my $ln=<HOMOZY>){
                chomp($ln);
		next unless $ln;
                my @data = split(',',$ln);
		next unless $data[1];
                if ($data[1] < 2){
 #                       warn $data[0]." : skipped due to 0 coverage affecting median\n";
                }else{
                        $hash_homozy{'median'}{$data[0]} = $data[1];
                        $hash_homozy{'coverage'}{$data[0]} = $data[2];
			$scaffold_size{$data[0]} = $data[3] if !$scaffold_size{$data[0]};
                        $hash_hetero{'coverage'}{$data[0]} = $hash_hetero{'coverage'}{$data[0]} ? $hash_hetero{'coverage'}{$data[0]} : int(0);
                }
        }
        close (HOMOZY);
	my @median_array_hetero = sort{$a <=> $b} (values %{$hash_hetero{'median'}});
	my @median_array_homozy = sort{$a <=> $b} (values %{$hash_homozy{'median'}});

	my $stat_x = Statistics::Descriptive::Full->new();
	$stat_x->add_data( @median_array_hetero );
	my $genome_median_hetero = sprintf("%.2f", $stat_x->trimmed_mean( 0.3 ));

	$stat_x = Statistics::Descriptive::Full->new();
	$stat_x->add_data( @median_array_homozy );
	my $genome_median_homozy = sprintf("%.2f", $stat_x->trimmed_mean( 0.3 ));

#	my $genome_median_hetero = &median(\@median_array_hetero,scalar(@median_array_hetero),20);

#	my $genome_median_homozy = &median(\@median_array_homozy,scalar(@median_array_homozy),20);

	my ($scale_factor,$genome_median,$depth_cutoff,$breakpoint,$smaller_dataset) = (int(0),int(0),int(0),int(0));
	if ($genome_median_homozy < $genome_median_hetero){
		$scale_factor = $no_scale ? 1 : $genome_median_homozy / $genome_median_hetero;
		$genome_median = $genome_median_hetero * $scale_factor;
		$breakpoint = sprintf("%.2f", $genome_median_homozy / 0.69315);
		$smaller_dataset = 'homozy';
	}else{
		$scale_factor = $no_scale ? 1 : $genome_median_hetero / $genome_median_homozy;
		$genome_median = $genome_median_homozy * $scale_factor;
		$breakpoint = sprintf("%.2f", $genome_median_hetero / 0.69315);
		$smaller_dataset = 'hetero';
	}
	$depth_cutoff = $user_max_depth_cutoff ? $user_max_depth_cutoff : ( $genome_median * (0.5 + $median_multiplier ) );
	my $cutoff = $breakpoint < $depth_cutoff ?  $breakpoint : $depth_cutoff;	

	my $scale_factor_round = sprintf("%.2f", $scale_factor);
	print "Pre downscaling genome wide trimmed mean: Heterozygous:$genome_median_hetero & Homozygous:$genome_median_homozy. The $smaller_dataset dataset is smaller so other is scaled by about $scale_factor_round. Depth cutoff used is $depth_cutoff. Breakpoint would be around $breakpoint\n\n";
	
	print "Processing X chromosome\n";
	print "Scaffold\tSize\tGaps\tHeteroCov\tHomozyCov\tHeteroDepth\tHomozyDepth\tScaledHeteroDepth\tScaledHomozyDepth\n";
	my $total_size = int(0);

	foreach my $scaff (sort {$scaffold_size{$b} <=> $scaffold_size{$a}} keys %scaffold_size){
		my $scaff_size = $scaffold_size{$scaff};
		die "Can't find $scaff" if !$scaff_size;
		my $gaps_proportion = sprintf("%.2f",($sizes_hash_ref->{$scaff}->{'gaps'} / $scaff_size)*100);
		next if $gaps_proportion > $gap_cutoff || $scaff_size < $size_cutoff;;

		my $median_hetero = $hash_hetero{'median'}{$scaff} ? $hash_hetero{'median'}{$scaff} : int(0);
		my $median_homozy = $hash_homozy{'median'}{$scaff} ? $hash_homozy{'median'}{$scaff} : int(0);
		
		my $scaled_median_hetero = $smaller_dataset eq 'homozy' ? sprintf("%.1f", ( $median_hetero * $scale_factor )) : $median_hetero;
		my $scaled_median_homozy = $smaller_dataset eq 'hetero' ? sprintf("%.1f", ( $median_homozy * $scale_factor )) : $median_homozy;

		if ($do_all || $user_seqs{$scaff} || ($median_hetero < $median_homozy && ($median_hetero < $depth_cutoff))){ # || $median_hetero < $breakpoint)){
			print "$scaff\t$scaff_size\t$gaps_proportion\t"
			.$hash_hetero{'coverage'}{$scaff}."\t"
			.$hash_homozy{'coverage'}{$scaff}."\t"
			."$median_hetero\t$median_homozy\t"
			."$scaled_median_hetero\t$scaled_median_homozy\n"
			
		}
		$total_size+=$scaff_size;
	}
	print "Total: ".&thousands($total_size)."\n\n";

	print "Processing Y chromosome\n";
	print "Scaffold\tSize\tGaps\tHeteroCov\tHomozyCov\tHeteroDepth\tHomozyDepth\tScaledHeteroDepth\tScaledHomozyDepth\n";
	$total_size = int(0);

	foreach my $scaff (sort {$scaffold_size{$b} <=> $scaffold_size{$a}} keys %scaffold_size){
		my $scaff_size = $scaffold_size{$scaff};
		die "Can't find $scaff" if !$scaff_size;
		my $gaps_proportion = sprintf("%.2f",($sizes_hash_ref->{$scaff}->{'gaps'} / $scaff_size)*100);
		next if $gaps_proportion > $gap_cutoff || $scaff_size < $size_cutoff;;

		my $median_hetero = $hash_hetero{'median'}{$scaff} ? $hash_hetero{'median'}{$scaff} : int(0);
		my $median_homozy = $hash_homozy{'median'}{$scaff} ? $hash_homozy{'median'}{$scaff} : int(0);
		
		my $scaled_median_hetero = $smaller_dataset eq 'homozy' ? sprintf("%.1f", ( $median_hetero * $scale_factor )) : $median_hetero;
		my $scaled_median_homozy = $smaller_dataset eq 'hetero' ? sprintf("%.1f", ( $median_homozy * $scale_factor )) : $median_homozy;

		if ($do_all || $user_seqs{$scaff} || ($median_homozy < 10 && $median_hetero > $median_homozy)){
			print "$scaff\t$scaff_size\t$gaps_proportion\t"
			.$hash_hetero{'coverage'}{$scaff}."\t"
			.$hash_homozy{'coverage'}{$scaff}."\t"
			."$median_hetero\t$median_homozy\t"
			."$scaled_median_hetero\t$scaled_median_homozy\n"
		}
		$total_size+=$scaff_size;
	}
	print "Total: ".&thousands($total_size)."\n\n";
}


sub process_align($$){
	my $file = shift;
	my $max_depth_cutoff = shift;
	my %scaffold_size;
	next unless -s $file;
	open (IN,$file);
	my (%median_hash);
	while (my $ln=<IN>){
		chomp($ln);
		my @data = split(',',$ln);
		if ($data[1] < 1){
			warn $data[0]." : skipped due to 0 coverage\n";
		}else{
			$median_hash{$data[0]}=$data[1];
			$scaffold_size{$data[0]} = $data[3];
		}
	}
	close (IN);
	my @median_array = sort{$a <=> $b} (values %median_hash);
	my $genome_median = &median(\@median_array,scalar(@median_array));
	$max_depth_cutoff = $user_max_depth_cutoff ? $user_max_depth_cutoff : ( ( $genome_median/2 ) + ( $genome_median * $median_multiplier ) );;
	print "Processing $file (genome median of medians = $genome_median). Max cutoff set as $max_depth_cutoff\n";
	print "Scaffold ID\tSize (bp)\tProportion of Gaps (%)\tDepth (fold coverage)\n";
	my $total_size = int(0);
	foreach my $scaff (sort keys %median_hash){
		my $scaff_size = $scaffold_size{$scaff};
		my $gaps_proportion = sprintf("%.2f",($sizes_hash_ref->{$scaff}->{'gaps'} / $scaff_size)*100);
		next unless $median_hash{$scaff} < $max_depth_cutoff && $scaff_size >= $size_cutoff && $gaps_proportion < $gap_cutoff;
		next if $max_depth_cutoff != $user_max_depth_cutoff && $median_hash{$scaff} < $user_max_depth_cutoff;
		print "$scaff\t$scaff_size\t$gaps_proportion\t".$median_hash{$scaff}."\n";
		$total_size+=$scaff_size;
	}

	print "Total: ".&thousands($total_size)."\n\n";
	if (@seqs_to_check){
		my $total_size_user = int(0);
		print "Checking user provided sequences\n";
		foreach my $scaff (sort @seqs_to_check){
			next unless exists $sizes_hash_ref->{$scaff};
			my $scaff_size = $sizes_hash_ref->{$scaff}->{'length'};
			my $gaps_proportion = sprintf("%.2f",($sizes_hash_ref->{$scaff}->{'gaps'} / $scaff_size)*100);
			print "$scaff\t$scaff_size\t$gaps_proportion\t".$median_hash{$scaff}."\n";
			$total_size_user+=$scaff_size;
		}
		print "User provided sequences total: $total_size_user\n\n";
	}
}


sub process_depth(){
	my $file = shift;
	my $out = $do_mean ? "$file.trimmean.processed"  :  "$file.processed";
	$out .= '.zeros' if $add_zeros;
	return $out if (-s $out);
	pod2usage "No genome data\n" unless $genome_fasta && -s $genome_fasta;
	print "Processing $file as $out\n";
	open (DEPTH,$file);
	open (OUT,">$out");
	my @depth_array;
	my $previous_id;
	my $counter++;
	while (my $ln=<DEPTH>){
		chomp($ln);
		next unless $ln && $ln!~/^#/;
		my @data = split("\t",$ln);
		if ( $previous_id && $previous_id ne $data[0] ) {
			if ($add_zeros){
				my $size = $sizes_hash_ref->{$previous_id}->{'length'};
				my $extra = $size - scalar(@depth_array);
				if ($extra && $extra > 0){
					print STDERR "Adding extra $extra 0s for $previous_id until $size...\n";
					for (my $i=0;$i<$extra;$i++){
						push(@depth_array,int(0));
					}
				}
			}
			@depth_array = sort{$a <=> $b} @depth_array;
			my $median;
			if ($do_mean){
				my $stat_x = Statistics::Descriptive::Full->new();
				$stat_x->add_data( @depth_array );
				$median = sprintf("%.2f", $stat_x->trimmed_mean( 0.3 ));
			}else{
				$median = &median(\@depth_array,scalar(@depth_array));
			}
			my $coverage = sprintf("%.2f", (scalar(@depth_array)/ $sizes_hash_ref->{$previous_id}->{'length'}));
			print OUT "$previous_id,$median,$coverage,".$sizes_hash_ref->{$previous_id}->{'length'}."\n";
			$previous_id = $data[0];
			@depth_array = ();
			$counter++;
			print STDERR "Processed $counter / $scaffold_count    \r" if $counter=~/0$/;
		}
		if ($use_bedgraph){
			#Btry_100        2201    2202    29
			for (my $i=$data[1];$i<$data[2];$i++){
				push(@depth_array, int($data[3]));
			}
		}else{
			push(@depth_array, int($data[2]));
		}
		$previous_id = $data[0];
	}
	# last one
	if ($add_zeros){
		my $size = $sizes_hash_ref->{$previous_id}->{'length'};
		my $extra = $size - scalar(@depth_array);
		if ($extra && $extra > 0){
			print STDERR "Adding extra $extra 0s for $previous_id until $size...\n";
			for (my $i=0;$i<$extra;$i++){
				push(@depth_array,int(0));
			}
		}
	}
	@depth_array = sort{$a <=> $b} @depth_array;
	my $median;
	if ($do_mean){
		my $stat_x = Statistics::Descriptive::Full->new();
		$stat_x->add_data( @depth_array );
		$median = sprintf("%.2f", $stat_x->trimmed_mean( 0.3 ));
	}else{
		$median = &median(\@depth_array,scalar(@depth_array));
	}
	my $coverage = sprintf("%.2f", (scalar(@depth_array)/ $sizes_hash_ref->{$previous_id}->{'length'}));
	print OUT "$previous_id,$median,$coverage,".$sizes_hash_ref->{$previous_id}->{'length'}."\n";
	$counter++;
	print "Processed $counter / $scaffold_count    \n\n";
	close DEPTH;
	close OUT;
	return $out;
}

sub median(){
    my $ref = shift; #sorted
    my $len = shift;
    return -1 unless $ref && $len;

    if($len%2){
        return $ref->[int($len/2)];
    }
    else{
        return ($ref->[int($len/2)-1] + $ref->[int($len/2)])/2;
    }
}

sub thousands(){
        my $val = shift;
        $val = sprintf("%.0f", $val);
        return $val if length($val)<4;
        1 while $val =~ s/(.*\d)(\d\d\d)/$1,$2/;
        return $val;
}

