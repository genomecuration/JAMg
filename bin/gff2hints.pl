#!/usr/bin/env perl

use strict;
use warnings;

my $gff       = shift;
my $golden    = shift;
die "<GFF_filename> [is golden]\n" unless $gff && -s $gff;
my $src       = $golden ? 'GLD' : 'XNT';
my $type_suffix = $golden ? '' : 'part';
my $priority = $golden ? 7 : 5;
my $delimiter = "\n\n";
my $orig_sep  = $/;
open( GFF, $gff )          || die $!;
open( OUT, ">$gff.hints" ) || die $!;
$/ = $delimiter;

RECORD: while ( my $record = <GFF> ) {
    my @lines = split( "\n", $record );

    if ($lines[0]=~/\bcDNA_match\b/){
	for ( my $i = 0 ; $i < scalar(@lines) ; $i++ ) {
		my @data = split( "\t", $lines[$i] );
		$data[2] = 'exonpart';
		$data[1] = 'PASA_assembly';
		my $grp='';
		if ($data[8]=~/ID=(\S+);/){
			$grp = "grp=".$1.';';
		}
		print OUT join( "\t", @data[ 0 .. 7 ] )
	          . "\tsrc=PASA;pri=3;$grp\n";
	}
	next RECORD;
    }

    #my $gene_line = $lines[0];
    my $mRNA_line = $lines[1];
    $mRNA_line =~ /ID=([^;]+)/;
    my $mRNA_id = $1;

    # CDS intron
    for ( my $i = 2 ; $i < scalar(@lines) ; $i++ ) {
      my @data = split( "\t", $lines[$i] );
      next unless $data[8];
      

      if ( $data[2] eq 'mRNA' ) {
        $data[2] = 'genicpart';
        print OUT join( "\t", @data[ 0 .. 7 ] )
          . "\tsrc=$src;pri=$priority;grp=$mRNA_id\n";
      } elsif (    $data[2] eq 'exon'
                || $data[2] eq 'CDS'
                || $data[2] eq 'intron' )
      {
        $data[2].=$type_suffix;
        print OUT join( "\t", @data[ 0 .. 7 ] )
          . "\tsrc=$src;pri=$priority;grp=$mRNA_id\n";
      } elsif ( $data[2] eq 'splice_junction' ) {
        $data[2] = $data[8] =~ /splice3/ ? 'ass' : 'dss';
        print OUT join( "\t", @data[ 0 .. 7 ] )
          . "\tsrc=$src;pri=$priority;grp=$mRNA_id\n";
      } elsif ( $data[2] eq 'UTR' ) {
        $data[2] = 'UTRpart';
        print OUT join( "\t", @data[ 0 .. 7 ] )
          . "\tsrc=$src;pri=$priority;grp=$mRNA_id\n";
      }
    }
  }
  close GFF;
  close OUT;
  $/ = $orig_sep;

