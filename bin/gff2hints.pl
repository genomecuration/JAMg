#!/usr/bin/perl -w
use strict;

  my $gff       = shift;
  my $golden    = shift;
  my $src       = $golden ? 'GLD' : 'XNT';
  my $type_suffix = $golden ? '' : 'part';
  my $priority = $golden ? 7 : 5;
  my $delimiter = "\n\n";
  my $orig_sep  = $/;
  open( GFF, $gff )          || die;
  open( OUT, ">$gff.hints" ) || die;
  $/ = $delimiter;

  while ( my $record = <GFF> ) {
    my @lines = split( "\n", $record );

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

