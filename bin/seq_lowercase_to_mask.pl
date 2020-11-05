#!/usr/bin/perl -w
package seq_uppercase;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::SeqIO;
my $infile;
my ( $outfile );
GetOptions( 
	'o|outfile:s' => \$outfile, 
	'in' => \$infile,
	 );
$infile=shift unless $infile;

unless ( $infile && -s $infile ) { pod2usage; }
if ( !$outfile ) { $outfile = $infile . ".hardmasked"; }

	my $filein  = new Bio::SeqIO( -file => $infile,     -format => 'fasta' );
	my $fileout = new Bio::SeqIO( -file => ">$outfile", -format => 'fasta' );
	while ( my $seq = $filein->next_seq() ) {
		my $sequence = $seq->seq();
		$sequence=~tr/[a-z]/N/;

		$seq->seq( $sequence );
		$fileout->write_seq($seq);
	}
