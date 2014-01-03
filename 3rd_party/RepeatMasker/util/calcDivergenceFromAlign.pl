#!/usr/local/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) calcDivergenceFromAlign.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A utility script to calculate a new divergence measure on the
##      RM alignment files.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2009 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
#
# ChangeLog
#
#     $Log: calcDivergenceFromAlign.pl,v $
#     Revision 1.16  2013/06/20 17:49:19  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
#
# To Do:
#

=head1 NAME

calcDivergenceFromAlign.pl - Recalculate the divergences in an align file.

=head1 SYNOPSIS

  calcDivergenceFromAlign.pl [-version] [-noCpG] *.align > new.align

=head1 DESCRIPTION

  A utility script to calculate a new divergence measure on the
  RM alignment files.  Currently we only calculate the Kimura 2-Parameter
  divergence metric.  

The options are:

=over 4

=item -version

Displays the version of the program

=item -noCpG

Ignore "CG" sites in the consensus sequence while calculating
the divergence of the sequence.

=back

=head1 SEE ALSO

ReapeatMasker

=head1 COPYRIGHT

Copyright 2009 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use FindBin;
use lib $FindBin::Bin;
use lib "$FindBin::Bin/../";
use Getopt::Long;
use Data::Dumper;
use CrossmatchSearchEngine;
use File::Basename;

#
# Version
#
my $Version = 0.1;
my $DEBUG   = 0;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number parameters
#
my @getopt_args = (
                    '-version',    # print out the version and exit
                    '-noCpG',
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n";
  exec "pod2text $0";
  exit( 1 );
}

if ( $options{'version'} ) {
  print "$Version\n";
  exit;
}

usage() if ( !$ARGV[ 0 ] );

my $annotationFile = $ARGV[ 0 ];

my $noCpG = undef;
$noCpG = 1 if ( $options{'noCpG'} );

#
# Open up a search results object
#
my $searchResults =
    CrossmatchSearchEngine::parseOutput( searchOutput => $annotationFile );

for ( my $i = 0 ; $i < $searchResults->size() ; $i++ ) {
  my $result = $searchResults->get( $i );

  my $consensus = $result->getSubjString();
  my $subject   = $result->getQueryString();

  print ""
      . $result->toStringFormatted( SearchResult::AlignWithQuerySeq ) . "\n";

  print "  KIMURA 2-Parameter Correction = "
      . kimuraDivergence( $consensus, $subject, $noCpG ) . "\n\n";

}

## TODO: This is now available from the SearchResult object directly.
sub kimuraDivergence {
  my $consensusSeq = shift;
  my $subjectSeq   = shift;
  my $noCpGs       = shift;

  my $alignedBases = 0;
  my $transI       = 0;
  my $transV       = 0;

  for ( my $i = 0 ; $i < length( $consensusSeq ) ; $i++ ) {
    my $cBase   = uc( substr( $consensusSeq, $i, 1 ) );
    my $CpGPair = uc( substr( $consensusSeq, $i, 2 ) );
    my $sBase   = uc( substr( $subjectSeq,   $i, 1 ) );

    # Gaps don't count towards aligned bases
    next if ( $cBase eq '-' || $sBase eq '-' );

    # Ignore CpG sites by not counting them in the aligned bases
    # and not counting their mismatches
    if ( defined $noCpGs && $noCpGs > 0 && $CpGPair eq "CG" ) {
      $i++;
      next;
    }

    # What to do with consensi containing IUB codes?  The following
    # if statement assumes that it's a match, increments the aligned
    # bases and moves on.
    if ( $cBase !~ /[ACGT]/ ) {
      $alignedBases++;
      next;
    }

    $transI++
        if (    $cBase . $sBase eq "CT"
             || $cBase . $sBase eq "TC"
             || $cBase . $sBase eq "AG"
             || $cBase . $sBase eq "GA" );
    $transV++
        if (    $cBase . $sBase eq "GT"
             || $cBase . $sBase eq "TG"
             || $cBase . $sBase eq "GC"
             || $cBase . $sBase eq "CG"
             || $cBase . $sBase eq "CA"
             || $cBase . $sBase eq "AC"
             || $cBase . $sBase eq "AT"
             || $cBase . $sBase eq "TA" );
    $alignedBases++;
  }

  return ( 0 ) if ( $alignedBases == 0 );

  my $p          = $transI / $alignedBases;
  my $q          = $transV / $alignedBases;
  my $logOperand = ( ( 1 - ( 2 * $p ) - $q ) * ( 1 - ( 2 * $q ) )**0.5 );

  print "Transitions = $transI Transverions = $transV "
      . "AlignedBases = $alignedBases\n";

  if ( $logOperand <= 0 ) {
    return ( 1 );
  }
  else {
    return ( abs( ( -0.5 * log( $logOperand ) ) ) );
  }

}

1;
