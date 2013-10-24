#!/usr/bin/perl
use strict; use warnings;
use IK;
use FAlite;
use DataBrowser;

die "usage: $0 <fasta> <k-mer>" unless @ARGV == 2;
my ($file, $K) = @ARGV;

# stage 1: get codon usage from longest ORF
my %codon;
open(IN, $file) or die;
my $fasta = new FAlite(\*IN);
while (my $entry = $fasta->nextEntry) {
	my $dna = uc $entry->seq;
	my @orf;
	for (my $i = 0; $i < 3; $i++) {
		my $tx1 = substr($dna, $i);
		my $tx2 = substr(IK::reverse_complement($dna), $i);
		my $aa1 = IK::translate($tx1);
		my $aa2 = IK::translate($tx2);
		push @orf, orfs($aa1, $tx1, '+');
		push @orf, orfs($aa2, $tx2, '-');
	}
	
	# sort by length only
	@orf = sort {length($b->{aa}) <=> length($a->{aa})} @orf;
	my $cds = $orf[0]{dna};
	for (my $i = 0; $i < length($cds); $i+=3) {
		$codon{substr($cds, $i, 3)}++;
	}
}
close IN;


# transform codon counts to log prob
my $total = 0;
foreach my $k (keys %codon) {$total += $codon{$k}}
foreach my $k (keys %codon) {$codon{$k} = log($codon{$k}/$total)}


# stage 2: keep longest ORF (ties broken by codon usage)
my @cds;
open(IN, $file) or die;
$fasta = new FAlite(\*IN);
while (my $entry = $fasta->nextEntry) {
	my $dna = uc $entry->seq;
	my @orf;
	for (my $i = 0; $i < 3; $i++) {
		my $tx1 = substr($dna, $i);
		my $tx2 = substr(IK::reverse_complement($dna), $i);
		my $aa1 = IK::translate($tx1);
		my $aa2 = IK::translate($tx2);
		push @orf, orfs($aa1, $tx1, '+');
		push @orf, orfs($aa2, $tx2, '-');
	}
	
	# find best ORF
	foreach my $orf (@orf) {$orf->{score} = score($orf->{dna})}
	@orf = sort {$b->{len} <=> $a->{len} or $b->{score} <=> $a->{score}} @orf;
	my $cds = $orf[0]{dna};
	push @cds, substr($cds, 6, length($cds) - 12);
}
close IN;


# stage 3: produce CDS model
my @model = (
	IK::blank_table($K, 1),
	IK::blank_table($K, 1),
	IK::blank_table($K, 1)
);
foreach my $cds (@cds) {
	for (my $frame = 0; $frame <= 2; $frame++) {
		for (my $i = $frame; $i < length($cds) -$K -$frame +1; $i+= 3) {
			my $kmer = substr($cds, $i, $K);
			$model[$frame]{$kmer}++;
			$total++;
		}
	}
}

#for (my $i = 0; $i < @model; $i++) {
#	print "frame$i\n";
#	my $count = 1;
#	foreach my $kmer (sort keys %{$model[$i]}) {
#		print "\t", $model[$i]{$kmer};
#		if ($count++ %4 == 0) {print "\n"}
#	}
#}
#die;

# output
my $table = IK::blank_table($K - 1);
my @alph = qw(A C G T);
print "Coding CDS 3 2 4 3 0.000\n";
for my $frame (0, 1, 2) {
	printf "\tframe%d LUT %d %d 4 0 0.000\n", $frame, $K, $K-1;

	foreach my $kmer (sort keys %$table) {
	
		my $total = 0;
		foreach my $nt (@alph) {
			$total += $model[$frame]{"$kmer$nt"};
		}
	
		foreach my $nt (@alph) {
			if ($total == 0) {
				print "\t\t.";
				next;
			}
			my $obs = $model[$frame]{"$kmer$nt"} / $total;
			if ($obs == 0) {
				print "\t\t.";
				next;
			}
			
			my $val = log($obs/0.25)/log(2);
			printf "\t\t%.3f", $val;
		}
		print "\n";
	}
}


##############################################################################

sub score {
	my ($dna) = @_;
	my $score = 0;
	for (my $i = 0; $i < length($dna); $i+=3) {
		my $k = substr($dna, $i, 3);
		if (exists $codon{$k}) {
			$score += $codon{$k};
		} else {
			$score += -10;
		}
	}
	return $score;
}


sub orfs {
	my ($pep, $dna, $s) = @_;
	my @orf;
	while ($pep =~ /(\w+)/g) {
		my $end = pos($pep);
		my $beg = $end - length($1);
		my $len = $end - $beg;
		my $orf = substr($dna, 3 * $beg, $len * 3);
		push @orf, {
			dna => $orf,
			aa => $1,
			strand => $s,
			len => length($1),
		};
	}
	return @orf;
}


__END__

In short transcripts, there may be more than 1 reasonable ORF. For this reason,
reads are processed in two steps. Step 1 selects the longest ORF and builds a
codon usage table. Step 2 also finds the longest ORF, but ties are broken by
similarity to the codon usage table.