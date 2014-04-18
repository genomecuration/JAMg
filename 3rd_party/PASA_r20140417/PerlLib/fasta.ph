#!/usr/local/bin/perl

use strict;

=head1 NAME

grasta.ph

=cut

=head1 DESCRIPTION

provides a subroutine for performing a grasta pairwise alignment between two sequences, protein or nucleotide.  It returns information for the top match as indicated in the subroutine specification below.

=cut




=item perform_fasta_align()

=over 4

B<Description:> Performs a pairwise alignment between two sequences (protein or nucleotide).

B<Parameters:> $seq1, $seq2

$seq1 and $seq2 are scalar values holding a raw line of text.  Do NOT provide a FASTA-formatted sequence.

B<Returns:> ($per_id, $overlap)

$per_id is the percent identity across the top match.
$overlap is the number of residues in this top alignment.

if the fasta program is not located in your PATH setting, set 
our \$FASTAPATH = /full/path/to/fasta

=back

=cut

    ;


our $FASTAPATH = "";
our $SEE;

####
sub perform_fasta_align {
    my $seq1 = shift;
    my $seq2 = shift;
    my $file1 = "/tmp/$$.seq1";
    open (SEQ1, ">$file1") or die;
    print SEQ1 ">seq1\n" . &FASTA_format($seq1);
    close SEQ1;
    my $file2 = "/tmp/$$.seq2";
    open (SEQ2, ">$file2") or die;
    print SEQ2 ">seq2\n" . &FASTA_format ($seq2);
    close SEQ2;
	
    my $result_file = "/tmp/$$.grasta_result";
    
    my $prog = `which fasta`;chomp($prog);
    if ($FASTAPATH) {
		$prog = $FASTAPATH;
    }
    die "Cannot find program fasta\n" unless $prog && -s $prog && -x $prog;
    
    my $cmd = "$prog $file1 $file2 > $result_file";
    my $ret = system ($cmd);
    if ($ret) {
		die "Error: couldn't perform alignment using prog $prog.\n";
    }
    
    ## Parse the output file
    open (OUTPUT, $result_file) or die "ERROR: Sorry, can't open $result_file";
    my $per_id = 0;
    my $overlap = 0;
    while (<OUTPUT>) {
		#print;
	
		if (/\s(\S+)% identity .* in (\d+) (nt|aa) overlap/) {
			$per_id = $1;
			$overlap = $2;
			last;
		}
		
	}
    close OUTPUT;
    print "FASTA: per_id: $per_id, overlap: $overlap\n" if $SEE;
    unlink ($file1, $file2, $result_file);
	
    return ($per_id, $overlap);
}


sub FASTA_format {
    my $seq = shift;
    $seq =~ s/(\w{60})/$1\n/g;
    return ($seq);
}


1; #end of ph.
