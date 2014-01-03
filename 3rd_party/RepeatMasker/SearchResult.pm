#!/usr/local/bin/perl -w
##---------------------------------------------------------------------------##
##  File:
##      @(#) SearchResult.pm
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      An for holding a generic biological sequence similarity
##      search result.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2004 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
# Implementation Details:
#
# bless( {
#          'qryGCBackground' => '43',
#          'querySeq' => 'AGCAA..TGTAAA',
#          'sbjEnd' => '2760',
#          'qryEnd' => '319032',
#          'qryBegin' => '318751',
#          'percDel' => '4.26',
#          'sbjBegin' => '2484',
#          'qryName' => 'ctg12382',
#          'score' => '1005',
#          'sbjName' => 'Charlie1',
#          'percIns' => '6.03',
#          'sbjOrient' => 'C',
#          'qryLeft' => '4582070',
#          'sbjLeft' => '1',
#          'id' => '5',
#          'uniqId' => '1383820',
#          'percDiv' => '16.31',
#          'subjSeq' => 'AGCGGT...AA',
#          'matrix' => '35p40g.matrix',
#          'transV' => '3',
#          'transI' => '1',
#          'overlap' => '*',
#        }
#     }, 'SearchResult' );
#
###############################################################################
# ChangeLog
#
#     $Log: SearchResult.pm,v $
#     Revision 1.99  2013/06/20 17:49:19  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
# To Do:
#    - Expose overlap, transitions and transversions attributes!
#
#

=head1 NAME

SearchResult

=head1 SYNOPSIS

use SearchResult

Usage: 

    $SearchResultsCollection = SearchResult->new();

  or 

    $SearchResultsCollection = SearchResult->new( 
                                 queryName=>value, subjName=>value,
                                 pctInsert=>value, pctDelete=>value, 
                                 queryStart=>value, queryEnd=>value, 
                                 score=>value, pctDiverge=>value, 
                                 subjRemaining=>value, subjType=>value,
                                 queryRemaining=>value, id=>value,
                                 orientation=>value, queryString=>value,
                                 subjString=>value, matrix=>value,
                                 id=>value, lineageId=>value );

=head1 DESCRIPTION

A class for storing a result from the crossmatch search engine.

=head1 SEE ALSO

=over 4

SearchResultCollection

The sequences from a cross_match alignment are stored relative to
the forward oriented query sequence.  I.e

  212 13.48 0.00 0.00 chr10 1  5 ( 0 ) + AluSp 1 5 (308 )
  chr10  1 AGCCG 5
            i
  AluSp  1 AACCG 5

  queryString = AGCCG
  subjString =  AACCG

However the following is reverse/complemented before storing:

  212 13.48 0.00 0.00 chr10 1  5 ( 0 ) C AluSp (308) 5 1 
  C chr10  5 ATCCG 1
              i
    AluSp  1 AACCG 5

  queryString = CGGAT
  subjString =  CGGTT


=back

=head1 COPYRIGHT

Copyright 2004 Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=head1 INSTANCE METHODS

=cut 

package SearchResult;
use strict;
use Data::Dumper;
use Matrix;
use Carp;
use vars qw(@ISA @EXPORT @EXPORT_OK %EXPORT_TAGS $VERSION);

use constant NoAlign => 1;

# Indicates that alignments should always be reported with
# the query sequences all in the forward direction
use constant AlignWithQuerySeq => 2;

# Indicates that alignments should always be reported with
# the subj sequences all in the forward direction
use constant AlignWithSubjSeq          => 3;
use constant OutFileFormat             => 4;
use constant CompressedAlignCSV        => 5;
use constant PSL                       => 6;
use constant RangeHighlightedAlignment => 7;

# Names to specify query/subject inputs
use constant Query   => 1;
use constant Subject => 2;

require Exporter;

@ISA = qw(Exporter);

@EXPORT = qw();

@EXPORT_OK = qw();

%EXPORT_TAGS = ( all => [ @EXPORT_OK ] );

my $CLASS = "SearchResult";
my $DEBUG = 0;

##-------------------------------------------------------------------------##
## Constructor:
##-------------------------------------------------------------------------##
sub new {
  my $class          = shift;
  my %nameValuePairs = @_;

  # Create ourself as a hash
  my $this = {};

  # Bless this hash in the name of the father, the son...
  bless $this, $class;

  # Allow import of values
  if ( %nameValuePairs ) {
    while ( my ( $name, $value ) = each( %nameValuePairs ) ) {
      my $method = "set" . _ucFirst( $name );
      unless ( $this->can( $method ) ) {
        croak(
             "SearchResult::add: Instance variable $name doesn't exist." . "" );
      }
      $this->$method( $value );
    }
  }

  return $this;
}

##-------------------------------------------------------------------------##

=head2 clone()

  Use: my $newObj = $obj->clone();

  Clone a SearchResult *duplicating* all the values of the old
  object in the new one.

=cut

##-------------------------------------------------------------------------##
sub clone {
  my $this = shift;

  my %newHash = %{$this};
  my $newObj  = \%newHash;

  bless $newObj, ref( $this );

  return $newObj;
}

##-------------------------------------------------------------------------##
## Get and Set Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2 get_setMatrixName()

  Use: my $value    = getMatrixName( );
  Use: my $oldValue = setMatrixName( $value );

  Get/Set the name of the matrix.

=cut

##-------------------------------------------------------------------------##
sub getMatrixName {
  my $obj = shift;

  my $value = $obj->{'matrixName'};

  return $value;
}

sub setMatrixName {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'matrixName'};
  $obj->{'matrixName'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setLineageId()

  Use: my $value    = getLineageId( );
  Use: my $oldValue = setLineageId( $value );

  Get/Set the lineage id.  NOTE: This alternate ID is used as the 
  refinment id by RepeatMasker/ProcessRepeats.  

=cut

##-------------------------------------------------------------------------##
sub getLineageId {
  my $obj = shift;

  my $value = $obj->{'lineageId'};

  return $value;
}

sub setLineageId {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'lineageId'};
  $obj->{'lineageId'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setQueryName()

  Use: my $value    = getQueryName();
  Use: my $oldValue = setQueryName( $value );

  Get/Set the name of the query sequence.

=cut

##-------------------------------------------------------------------------##
sub getQueryName {
  my $obj = shift;

  my $value = $obj->{'qryName'};

  return $value;
}

sub setQueryName {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'qryName'};
  $obj->{'qryName'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setSubjName()

  Use: my $value    = getSubjName();
  Use: my $oldValue = setSubjName( $value );

  Get/Set the name of the subject sequence.

=cut

##-------------------------------------------------------------------------##
sub getSubjName {
  my $obj = shift;

  my $value = $obj->{'sbjName'};

  return $value;
}

sub setSubjName {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'sbjName'};
  $obj->{'sbjName'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setSubjType()

  Use: my $value    = getSubjType();
  Use: my $oldValue = setSubjType( $value );

  Get/Set the type of the subject sequence ( for out files only ).

=cut

##-------------------------------------------------------------------------##
sub getSubjType {
  my $obj = shift;

  my $value = $obj->{'sbjType'};

  return $value;
}

sub setSubjType {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'sbjType'};
  $obj->{'sbjType'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setPctInsert()

  Use: my $value    = getPctInsert();
  Use: my $oldValue = getPctInsert( $value );

  Get/Set the percent insert value.

=cut

##-------------------------------------------------------------------------##
sub getPctInsert {
  my $obj = shift;

  my $value = $obj->{'percIns'};

  return $value;
}

sub setPctInsert {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'percIns'};
  $obj->{'percIns'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setPctDelete()

  Use: my $value    = getPctDelete();
  Use: my $oldValue = setPctDelete( $value );

  Get/Set the percent deletion value.

=cut

##-------------------------------------------------------------------------##
sub getPctDelete {
  my $obj = shift;

  my $value = $obj->{'percDel'};

  return $value;
}

sub setPctDelete {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'percDel'};
  $obj->{'percDel'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setSubjStart()

  Use: my $value    = getSubjStart();
  Use: my $oldValue = setSubjStart( $value );

  Get/Set the subject start position (1 based).

=cut

##-------------------------------------------------------------------------##
sub getSubjStart {
  my $obj = shift;

  my $value = $obj->{'sbjBegin'};

  return $value;
}

sub setSubjStart {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'sbjBegin'};
  $obj->{'sbjBegin'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setQueryStart()

  Use: my $value    = getQueryStart();
  Use: my $oldValue = setQueryStart( $value );

  Get/Set the query start position (1 based).

=cut

##-------------------------------------------------------------------------##
sub getQueryStart {
  my $obj = shift;

  my $value = $obj->{'qryBegin'};

  return $value;
}

sub setQueryStart {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'qryBegin'};
  $obj->{'qryBegin'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setQueryEnd()

  Use: my $value    = getQueryEnd();
  Use: my $oldValue = setQueryEnd( $value );

  Get/Set the query end position (1 based).

=cut

##-------------------------------------------------------------------------##
sub getQueryEnd {
  my $obj = shift;

  my $value = $obj->{'qryEnd'};

  return $value;
}

sub setQueryEnd {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'qryEnd'};
  $obj->{'qryEnd'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setSubjEnd()

  Use: my $value    = getSubjEnd();
  Use: my $oldValue = setSubjEnd( $value );

  Get/Set the subject end position (1 based).

=cut

##-------------------------------------------------------------------------##
sub getSubjEnd {
  my $obj = shift;

  my $value = $obj->{'sbjEnd'};

  return $value;
}

sub setSubjEnd {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'sbjEnd'};
  $obj->{'sbjEnd'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setScore()

  Use: my $value    = getScore();
  Use: my $oldValue = setScore( $value );

  Get/Set the score value.

=cut

##-------------------------------------------------------------------------##
sub getScore {
  my $obj = shift;

  my $value = $obj->{'score'};

  return $value;
}

sub setScore {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'score'};
  $obj->{'score'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setBitScore()

  Use: my $value    = getBitScore();
  Use: my $oldValue = setBitScore( $value );

  Get/Set the bit score value.

=cut

##-------------------------------------------------------------------------##
sub getBitScore {
  my $obj = shift;

  my $value = $obj->{'bitScore'};

  return $value;
}

sub setBitScore {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'bitScore'};
  $obj->{'bitScore'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setEValue()

  Use: my $value    = getEvalue();
  Use: my $oldValue = setEvalue( $value );

  Get/Set the Evalue.

=cut

##-------------------------------------------------------------------------##
sub getEvalue {
  my $obj = shift;

  my $value = $obj->{'evalue'};

  return $value;
}

sub setEvalue {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'evalue'};
  $obj->{'evalue'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setPValue()

  Use: my $value    = getPValue();
  Use: my $oldValue = setPValue( $value );

  Get/Set the PValue.

=cut

##-------------------------------------------------------------------------##
sub getPValue {
  my $obj = shift;

  my $value = $obj->{'PValue'};

  return $value;
}

sub setPValue {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'PValue'};
  $obj->{'PValue'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setPctDiverge()

  Use: my $value    = getPctDiverge();
  Use: my $oldValue = setPctDiverge( $value );

  Get/Set the percent divergence value.

=cut

##-------------------------------------------------------------------------##
sub getPctDiverge {
  my $obj = shift;

  my $value = $obj->{'percDiv'};

  return $value;
}

sub setPctDiverge {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'percDiv'};
  $obj->{'percDiv'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setSubjRemaining()

  Use: my $value    = getSubjRemaining();
  Use: my $oldValue = setSubjRemaining( $value );

  Get/Set the subject remaining length value.

=cut

##-------------------------------------------------------------------------##
sub getSubjRemaining {
  my $obj = shift;

  my $value = $obj->{'sbjLeft'};

  return $value;
}

sub setSubjRemaining {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'sbjLeft'};
  $obj->{'sbjLeft'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setQueryRemaining()

  Use: my $value    = getQueryRemaining();
  Use: my $oldValue = setQueryRemaining( $value );

  Get/Set the query remaining length value.

=cut

##-------------------------------------------------------------------------##
sub getQueryRemaining {
  my $obj = shift;

  my $value = $obj->{'qryLeft'};

  return $value;
}

sub setQueryRemaining {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'qryLeft'};
  $obj->{'qryLeft'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setOverlap()

  Use: my $value    = getOverlap();
  Use: my $oldValue = setOverlap( $value );

  Get/Set the the overlap value.

=cut

##-------------------------------------------------------------------------##
sub getOverlap {
  my $obj = shift;

  my $value = $obj->{'overlap'};

  return $value;
}

sub setOverlap {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'overlap'};
  $obj->{'overlap'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setId()

  Use: my $value    = getId();
  Use: my $oldValue = setId( $value );

  Get/Set the the identifier value.

=cut

##-------------------------------------------------------------------------##
sub getId {
  my $obj = shift;

  my $value = $obj->{'id'};

  return $value;
}

sub setId {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'id'};
  $obj->{'id'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setOrientation()

  Use: my $value    = getOrientation();
  Use: my $oldValue = setOrientation( $value );

  Get/Set the orientation of the sequence.  The orientation is 
  interpreted as the orientation of the subject sequence.
  The query is always assumed to be in the forward direction.

=cut

##-------------------------------------------------------------------------##
sub getOrientation {
  my $obj = shift;

  my $value = $obj->{'sbjOrient'};

  return $value;
}

sub setOrientation {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'sbjOrient'};
  $obj->{'sbjOrient'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setQueryString()

  Use: my $value    = getQueryString();
  Use: my $oldValue = setQueryString( $value );

  Get/Set the query portion of the alignment string.

=cut

##-------------------------------------------------------------------------##
sub getQueryString {
  my $obj = shift;

  my $value = $obj->{'querySeq'};

  return $value;
}

sub setQueryString {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'querySeq'};
  $obj->{'querySeq'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setSubjString()

  Use: my $value    = getSubjString();
  Use: my $oldValue = setSubjString( $value );

  Get/Set the subject portion of the alignment string.

=cut

##-------------------------------------------------------------------------##
sub getSubjString {
  my $obj = shift;

  my $value = $obj->{'subjSeq'};

  return $value;
}

sub setSubjString {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'subjSeq'};
  $obj->{'subjSeq'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##

=head2 get_setPPString()

  Use: my $value    = getPPString();
  Use: my $oldValue = setPPString( $value );

  Get/Set the HMMER PP portion of the alignment string.

=cut

##-------------------------------------------------------------------------##
sub getPPString {
  my $obj = shift;

  my $value = $obj->{'ppSeq'};

  return $value;
}

sub setPPString {
  my $obj      = shift;
  my $value    = shift;
  my $oldValue = undef;

  $oldValue = $obj->{'ppSeq'};
  $obj->{'ppSeq'} = $value;

  return $oldValue;
}

##-------------------------------------------------------------------------##
## General Object Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##

=head2

  Use: my $trimmed = trimHMMERAlignment( threshold => # );

  threshold:  Percentage threshold

  The APP line of HMMER records a per-position metric 
  for how confident nhmmer is of the position being part of the
  core model.  

      Currently ( beta release ):
             0 : <70%, 
             1 :  70 - <73% ... 
             2 :  73 - <76
             3 :  76-79
             4 :  79-82% ... 
             7 :  88-91% ... 
                  91-94
                  94-97
                      94 8
                      95 8.33
                      97 9
             * : >97%.  

      Planned Final Release:
             L : < 80%
             M : 80-90%
             H : >90% 
         The ranges for L,M,H will be configurble so
         we need some way of knowing what they have been
         set to.

=cut

##-------------------------------------------------------------------------##
sub trimHMMERAlignment {
  my $this           = shift;
  my %nameValuePairs = @_;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  unless ( exists $nameValuePairs{'threshold'} ) {
    croak "$subroutine: \"threshold\" parameter missing!\n";
  }
  my $thresh      = $nameValuePairs{'threshold'};
  my $minAPPValue = "0";
  if ( $thresh >= 70 && $thresh < 97 ) {
    $minAPPValue = int( ( $thresh - 70 ) / 3 ) + 1;
  }
  elsif ( $thresh == 97 ) {
    $minAPPValue = "9";
  }
  else {
    $minAPPValue = "*";
  }

  my $leftTrim  = 0;
  my $rightTrim = 0;

  my $ppSeq  = $this->getAPPString();
  my $qrySeq = $this->getQueryString();
  my $sbjSeq = $this->getSubjString();
  if ( $ppSeq eq "" ) {
    croak "$subroutine: This search result doesn't contain a APP sequence.\n"
        . "Perhaps this is not a HMMER alignment?\n";
  }

  # leftTrim
  for ( my $i = 0 ; $i < length( $ppSeq ) ; $i++ ) {
    my $APPValue = substr( $ppSeq, $i, 1 );
    last if ( $APPValue eq "*" || $APPValue >= $minAPPValue );
    $leftTrim++;
  }

  my $qBase = substr( $qrySeq, $leftTrim, 1 );
  my $sBase = substr( $sbjSeq, $leftTrim, 1 );

  if ( $qBase ne $sBase ) {

    # Move right
    while ( $qBase ne $sBase && $leftTrim <= length( $ppSeq ) ) {
      $leftTrim++;
      $qBase = substr( $qrySeq, $leftTrim, 1 );
      $sBase = substr( $sbjSeq, $leftTrim, 1 );
    }
  }
  else {

    # Move left
    while ( $qBase eq $sBase && $leftTrim > 0 ) {
      $leftTrim--;
      $qBase = substr( $qrySeq, $leftTrim, 1 );
      $sBase = substr( $sbjSeq, $leftTrim, 1 );
    }
    $leftTrim++ if ( $qBase ne $sBase );
  }

  # rightTrim
  for ( my $i = length( $ppSeq ) ; $i >= 0 ; $i-- ) {
    my $APPValue = substr( $ppSeq, $i, 1 );
    last if ( $APPValue eq "*" || $APPValue >= $minAPPValue );
    $rightTrim++;
  }

  $qBase = substr( $qrySeq, length( $ppSeq ) - $rightTrim - 1, 1 );
  $sBase = substr( $sbjSeq, length( $ppSeq ) - $rightTrim - 1, 1 );

  if ( $qBase ne $sBase ) {

    # Move left
    while ( $qBase ne $sBase && $rightTrim <= length( $ppSeq ) ) {
      $rightTrim++;
      $qBase = substr( $qrySeq, length( $ppSeq ) - $rightTrim - 1, 1 );
      $sBase = substr( $sbjSeq, length( $ppSeq ) - $rightTrim - 1, 1 );
    }
  }
  else {

    # Move right
    while ( $qBase eq $sBase && $rightTrim > 0 ) {
      $rightTrim--;
      $qBase = substr( $qrySeq, length( $ppSeq ) - $rightTrim - 1, 1 );
      $sBase = substr( $sbjSeq, length( $ppSeq ) - $rightTrim - 1, 1 );
    }
    $rightTrim++ if ( $qBase ne $sBase );
  }

  if ( ( $leftTrim + $rightTrim ) >= length( $qrySeq ) ) {
    return ( -1 );
  }

  my $leftQuery = substr( $qrySeq, 0, $leftTrim );
  my $leftQueryGap = ( $leftQuery =~ tr/-/-/ );
  my $rightQuery = substr( $qrySeq, length( $qrySeq ) - $rightTrim - 1 );
  my $rightQueryGap = ( $rightQuery =~ tr/-/-/ );

  $qrySeq =
      substr( $qrySeq, $leftTrim, length( $qrySeq ) - $leftTrim - $rightTrim );

  my $leftSubj = substr( $sbjSeq, 0, $leftTrim );
  my $leftSubjGap = ( $leftSubj =~ tr/-/-/ );
  my $rightSubj = substr( $sbjSeq, length( $sbjSeq ) - $rightTrim - 1 );
  my $rightSubjGap = ( $rightSubj =~ tr/-/-/ );
  $sbjSeq =
      substr( $sbjSeq, $leftTrim, length( $sbjSeq ) - $leftTrim - $rightTrim );

#print "TRIM: $leftTrim - $rightTrim $leftQuery ( $leftQueryGap, $leftSubjGap ) $rightQuery ( $rightQueryGap, $rightSubjGap )\n";
  $this->setQueryString( $qrySeq );
  $this->setQueryStart( $this->getQueryStart() + $leftTrim - $leftQueryGap );
  $this->setQueryEnd( $this->getQueryEnd() - $rightTrim + $rightQueryGap );
  $this->setSubjString( $sbjSeq );
  if ( $this->getOrientation ne "C" ) {
    $this->setSubjStart( $this->getSubjStart() + $leftTrim - $leftSubjGap );
    $this->setSubjEnd( $this->getSubjEnd() - $rightTrim + $rightSubjGap );
  }
  else {
    $this->setSubjStart( $this->getSubjStart() + $rightTrim - $rightSubjGap );
    $this->setSubjEnd( $this->getSubjEnd() - $leftTrim + $leftSubjGap );
  }
  return ( $leftTrim + $rightTrim );
}

##-------------------------------------------------------------------------##

=head2

  Use: my $count = getSeedPatternCount( $bitPattern );

    $bitPattern    : The seed pattern to run over the alignment.
                     ie. "10101011110111" or "111111111"

  Run the seed pattern over the alignment and calculate the number of
  seed matches using given a seed bit pattern. This can be used to
  determine the sensitivity of a seed pattern on a given alignment.

  NOTE: The un-interrupted seed may not be visible in the final alignment
        due to gaps being inserted.  Use the ungapped sequences to 
        find the potential seed matches.
=cut

##-------------------------------------------------------------------------##
sub getSeedPatternCount {
  my $obj        = shift;
  my $bitPattern = shift;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  my $query = $obj->{querySeq};
  $query =~ s/-//g;
  my $subject = $obj->{subjSeq};
  $subject =~ s/-//g;

  my $patLen     = length( $bitPattern );
  my $matchCount = 0;

  # Blast uses the query to build it's word hash
  my %queryWordHash = ();
  for ( my $i = 0 ; $i < length( $query ) - $patLen ; $i++ ) {
    my $word = "";
    for ( my $j = 0 ; $j < $patLen ; $j++ ) {
      my $patBit = substr( $bitPattern, $j, 1 );
      $word .= substr( $query, $i + $j, 1 ) if ( $patBit eq "1" );
    }
    $queryWordHash{$word}++;
  }

  #print "$query\n";
  for ( my $i = 0 ; $i < length( $subject ) - $patLen ; $i++ ) {
    my $word = "";
    for ( my $j = 0 ; $j < $patLen ; $j++ ) {
      my $patBit = substr( $bitPattern, $j, 1 );
      $word .= substr( $subject, $i + $j, 1 ) if ( $patBit eq "1" );
    }
    $matchCount++ if ( defined $queryWordHash{$word} );

    #print " "x($i) . $bitPattern . "\n" if ( defined $queryWordHash{$word} );
  }

  return ( $matchCount );
}

##-------------------------------------------------------------------------##

=head2

  Use: toStringFormatted( $format );

    $format        : SearchResult::NoAlign
                     SearchResult::AlignWithQuerySeq
                     SearchResult::AlignWithSubjSeq
                     SearchResult::OutFileFormat
                     SearchResult::CompressedAlignCSV
                     SearchResult::RangeHighlightedAlignment

  Create an string representation of a single result.

=cut

##-------------------------------------------------------------------------##
sub toStringFormatted {
  my $obj           = shift;
  my $format        = shift;
  my $displayParams = shift;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  $format = SearchResult::NoAlign
      if ( !defined $format );

  if (    $format == SearchResult::NoAlign
       || $format == SearchResult::AlignWithQuerySeq
       || $format == SearchResult::AlignWithSubjSeq
       || $format == SearchResult::RangeHighlightedAlignment )
  {
    return $obj->_toCrossMatchFormat( $format, $displayParams );
  }
  elsif ( $format == SearchResult::OutFileFormat ) {
    return $obj->_toOUTFileFormat();
  }
  elsif ( $format == SearchResult::CompressedAlignCSV ) {
    return $obj->_toCSVFormat( $displayParams );
  }
  else {
    croak $CLASS . "::toStringFormatted: Unknown format " . "( $format )\n";
  }
}

##-------------------------------------------------------------------------##

=head2

  Use: parseFromCSVFormat( $csvString );

  Populate object with values stored in CSV format.

=cut

##-------------------------------------------------------------------------##
sub parseFromCSVFormat {
  my $record = shift;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  my @flds   = split( /\,/, $record );
  my $rawSeq = $flds[ 16 ];

  my $inSub  = 0;
  my $inDel  = 0;
  my $inIns  = 0;
  my $qrySeq = "";
  my $sbjSeq = "";
  while ( $rawSeq =~ s/^(\S)// ) {
    my $char = $1;
    if ( $char eq "/" ) {
      $inSub = 1;
    }
    elsif ( $char eq "-" ) {
      $inDel ^= 1;
    }
    elsif ( $char eq "+" ) {
      $inIns ^= 1;
    }
    else {
      if ( $inSub == 1 ) {
        substr( $sbjSeq, length( $sbjSeq ) - 1, 1 ) = $char;
        $inSub = 0;
      }
      elsif ( $inDel == 1 ) {
        $qrySeq .= $char;
        $sbjSeq .= "-";
      }
      elsif ( $inIns == 1 ) {
        $qrySeq .= "-";
        $sbjSeq .= $char;
      }
      else {
        $qrySeq .= $char;
        $sbjSeq .= $char;
      }
    }
  }

  my $orient = "";
  if ( $flds[ 13 ] == 1 ) {
    $orient = "C";
  }

  my $retVal = SearchResult->new(
                                  score          => $flds[ 0 ],
                                  pctDiverge     => $flds[ 1 ],
                                  pctDelete      => $flds[ 2 ],
                                  pctInsert      => $flds[ 3 ],
                                  queryName      => $flds[ 4 ],
                                  queryStart     => $flds[ 5 ],
                                  queryEnd       => $flds[ 6 ],
                                  queryRemaining => $flds[ 7 ],
                                  subjName       => $flds[ 8 ],
                                  subjType       => $flds[ 9 ],
                                  subjStart      => $flds[ 10 ],
                                  subjEnd        => $flds[ 11 ],
                                  subjRemaining  => $flds[ 12 ],
                                  overlap        => $flds[ 14 ],
                                  id             => $flds[ 15 ],
                                  orientation    => $orient,
                                  queryString    => $qrySeq,
                                  subjString     => $sbjSeq
  );
  return $retVal;
}

##-------------------------------------------------------------------------##

=head2

  Use: my ( \@newSearchResults ) = fragmentSearchResult( 
                     regionList => \@ranges );

=cut

##-------------------------------------------------------------------------##
sub fragmentSearchResult {
  my $this            = shift;
  my %nameValueParams = @_;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  my $regionList;
  unless ( defined( $regionList = $nameValueParams{'regionList'} ) ) {
    croak "$subroutine: Error missing regionList parameters!\n";
  }

  my @newResults = ();
  for ( my $i = 0 ; $i < $#{$regionList} ; $i += 2 ) {
    my $alignStart = $regionList->[ $i ];
    my $alignEnd   = $regionList->[ $i + 1 ];
    my $newResult  = $this->clone();
    my $seqLength = $newResult->getQueryRemaining() + $newResult->getQueryEnd();

    my $qryPos = $newResult->getQueryStart();
    my $sbjPos = $newResult->getSubjStart();
    my $qrySeq = $newResult->getQueryString();
    my $sbjSeq = $newResult->getSubjString();
    if (    $alignStart < 0
         || $alignEnd > length( $qrySeq )
         || $alignEnd > length( $sbjSeq ) )
    {
      croak "$subroutine: Error range is outside alignment boundaries"
          . " $alignStart-$alignEnd ( alignment seqlen = "
          . length( $qrySeq ) . " / "
          . length( $sbjSeq ) . "\n";
    }

    my $j = 0;
    while ( $j <= $alignEnd ) {
      if ( $j == $alignStart ) {
        $newResult->setQueryStart( $qryPos );
        $newResult->setSubjStart( $sbjPos );
      }
      if ( $j == $alignEnd ) {
        $newResult->setQueryEnd( $qryPos );
        $newResult->setSubjEnd( $sbjPos );
        last;
      }
      $qryPos++ if ( substr( $qrySeq, $j, 1 ) ne "-" );
      $sbjPos++ if ( substr( $sbjSeq, $j, 1 ) ne "-" );
      $j++;
    }
    $newResult->setQueryRemaining( $seqLength - $newResult->getQueryEnd() );
    $newResult->setQueryString(
                  substr( $qrySeq, $alignStart, $alignEnd - $alignStart + 1 ) );
    $newResult->setSubjString(
                  substr( $sbjSeq, $alignStart, $alignEnd - $alignStart + 1 ) );
    push @newResults, $newResult;
  }

  return ( \@newResults );

}

##-------------------------------------------------------------------------##

=head2

  Use: my $bitScore = &rawToBitScore( $lambda, $mu );

  Calculate the bitscore given using the rawscore stored in the
  searchresult and the scoring system's lambda, and mu parameters.

=cut

##-------------------------------------------------------------------------##
sub rawToBitScore {
  my $this   = shift;
  my $lambda = shift;
  my $mu     = shift;

  return ( ( ( $this->getScore() * $lambda ) - log( $mu ) ) / log( 2 ) );
}

##-------------------------------------------------------------------------##

=head2

  Use: my ( $score, $divergence, $cpgsites, $percIns, $percDel, 
            \@positionScores, \@xdrop_fragments )  = rescoreAlignment( 
                             scoreMatrix => $matrixFileName,
                             gapOpenPenalty => #,
                             gapExtPenalty => #,
                             [noCpG => 1],
                             [xDrop => 1],
                             [complexityAdjust => 1] );

  Use the provided scoring system to rescore the alignment data
  stored in the object.  Does not alter objects values.  This
  routine will report the score ( complexity_adjusted if specified )
  , the Kimura divergence ( note - not the same as CrossMatch reports ),
  the number of CpG sites in alignment, the percent insertions, and 
  percent deletions.  

=cut

##-------------------------------------------------------------------------##
sub rescoreAlignment {
  my $this            = shift;
  my %nameValueParams = @_;

  my $subroutine = ( caller( 0 ) )[ 0 ] . "::" . ( caller( 0 ) )[ 3 ];

  my $matrix;
  unless ( defined( $matrix = $nameValueParams{'scoreMatrix'} )
           && $matrix->isa( "Matrix" ) )
  {
    croak "$subroutine: Error rescore alignment requires a "
        . "Matrix object as the scoreMatrix parameter!\n";
  }

  my $gapOpen      = $nameValueParams{'gapOpenPenalty'};
  my $gapExtension = $nameValueParams{'gapExtensionPenalty'};
  my $noCpG        = $nameValueParams{'noCpG'};

  my $qrySeq = $this->getQueryString();
  my $sbjSeq = $this->getSubjString();

  if ( $qrySeq eq "" || $sbjSeq eq "" ) {
    croak "$subroutine: Missing alignment data:\n" . ""
        . $this->toStringFormatted( SearchResult::NoAlign ) . "\n";
  }

  my $DEBUG = 0;

  my $matScoresRef     = $matrix->getMatrixValuesRef();
  my $matAlphabet      = join( "", @{ $matrix->getAlphabetRef() } );
  my $score            = 0;
  my $ungappedRawScore = 0;
  my @matCounts        = ();
  my $inGap            = "";
  my $transitions      = 0;
  my $transversions    = 0;
  my $CpGSites         = 0;
  my $cIndex           = index( $matAlphabet, "C" );
  my $gIndex           = index( $matAlphabet, "G" );
  my $cgScore          = $$matScoresRef[ $cIndex ][ $cIndex ] +
      $$matScoresRef[ $gIndex ][ $gIndex ];

  my $wellCharacterizedBases = 0;

  # Deletions ( Assuming Subject is conensus )
  my $deletionInits = 0;
  my $deletionExtns = 0;

  # Insertions ( Assuming Subject is conensus )
  my $insertionInits = 0;
  my $insertionExtns = 0;

  # Rescore
  my $i                   = 0;
  my @positionScores      = ();
  my $deletionOpenCounted = 0;
  while ( $i < length( $qrySeq ) ) {
    my $sbjSeqTrimmed = substr( $sbjSeq, $i );
    my $qrySeqTrimmed = substr( $qrySeq, $i );

    my $qBase = substr( $qrySeq, $i, 1 );
    my $sBase = substr( $sbjSeq, $i, 1 );

    # Assumption: Subject is the consensus sequence
    # Assumption: Matrix rows = ancestral state, cols = derived state
    #                or matrix[$sIdx][$qIdx]
    my $qIdx = index( $matAlphabet, $qBase );
    my $sIdx = index( $matAlphabet, $sBase );

    # Tabulate CpGSites even if we aren't ignoring them
    if ( $sbjSeqTrimmed =~ /^C([-]*)G.*/ ) {
      $CpGSites++;
    }

    if ( $noCpG && $sbjSeqTrimmed =~ /^C([-]*)G.*/ ) {
      my $gapLenCpG = length( $1 );
      my $CpGSubst  =
          substr( $qrySeq, $i, 1 ) . substr( $qrySeq, $i + $gapLenCpG + 1, 1 );

      #print "Comparing $qBase:$sBase $matAlphabet $qIdx $sIdx\n" if ( $DEBUG );
      # First handle the insertion
      if ( $gapLenCpG ) {
        $insertionInits++;
        $score += $gapOpen;
        push @positionScores, $score;
        for ( my $j = 0 ; $j < $gapLenCpG - 1 ; $j++ ) {
          $insertionExtns++;
          $score += $gapExtension;
          push @positionScores, $score;
        }
      }

      # Now handle deletions starting at one of these two bases
      if ( $CpGSubst =~ /^(.-|-.|--)$/ ) {
        my $numDel = ( $CpGSubst =~ tr/-/-/ );
        $deletionInits++;
        $score += $gapOpen;
        push @positionScores, $score;
        for ( my $j = 0 ; $j < $numDel - 1 ; $j++ ) {
          $deletionExtns++;
          $score += $gapExtension;
          push @positionScores, $score;
        }
        $deletionOpenCounted = 1;
      }

      my $qIdx2 = index( $matAlphabet, substr( $CpGSubst, 1, 1 ) );

      # Subtract one transition out of the total at a
      # CpG site ( owing to the fact that transitions are
      # really really common at CpG sites but after the first
      # transition the second could be counted as if it
      # occured at a non-CpG site ).  Transversions at
      # CpG sites are always counted.

      my $posScore;

      #$posScore += $$matScoresRef[ $qIdx ][ $sIdx ]
      $posScore += $$matScoresRef[ $sIdx ][ $qIdx ]
          if ( $qBase ne "-" );

      #$posScore += $$matScoresRef[ $qIdx2 ][ $gIndex ]
      $posScore += $$matScoresRef[ $gIndex ][ $qIdx2 ]
          if ( substr( $CpGSubst, 1, 1 ) ne "-" );

      if ( $CpGSubst eq "TA" ) {

        # Two transitions count as one
        $transitions++;
      }
      elsif ( $CpGSubst =~ /T./ || $CpGSubst =~ /.A/ ) {

        # One transition counts as 1/10 of a normal transition at these sites
        $transitions += 1 / 10;
      }
      elsif ( $CpGSubst =~ /[AG]G/ ) {
        $transversions++;
      }
      elsif ( $CpGSubst =~ /C[TC]/ ) {
        $transversions++;
      }
      elsif ( $CpGSubst =~ /[AG][TC]/ ) {
        $transversions += 2;
      }

      $ungappedRawScore += $posScore;
      $score            += $posScore;
      push @positionScores, $score;

      if ( $CpGSubst =~ /([ACGT]*)/ ) {
        $wellCharacterizedBases += length( $1 );
      }

      $matCounts[ $qIdx ]++  if ( $qIdx >= 0 );
      $matCounts[ $qIdx2 ]++ if ( $qIdx2 >= 0 );

      # Messing with the loop counter here. We processed the
      # 2+ positions here.
      $i += $gapLenCpG + 2;

    }
    else {
      if ( $sbjSeqTrimmed =~ /^(-[-]*)[^-].*/ ) {

        # Insertion
        print "Insertion init = $gapOpen  ext = " if ( $DEBUG );
        $insertionInits += 1;
        $insertionExtns += ( length( $1 ) - 1 );
        $score          += $gapOpen;
        push @positionScores, $score;
        for ( my $j = 0 ; $j < length( $1 ) - 1 ; $j++ ) {
          print "" . $gapExtension . " " if ( $DEBUG );
          $score += $gapExtension;
          push @positionScores, $score;
        }
        print "\n" if ( $DEBUG );
        $i += length( $1 );
      }
      elsif ( $qrySeqTrimmed =~ /^(-[-]*)[^-].*/ ) {

        # Deletion
        print "Deletion " if ( $DEBUG );

        # Account for deletion started in CpG above
        if ( !$deletionOpenCounted ) {
          $deletionInits += 1;
          print "init = $gapOpen ext = " if ( $DEBUG );
          $score += $gapOpen;
        }
        else {
          $deletionOpenCounted = 0;
        }
        push @positionScores, $score;
        for ( my $j = 0 ; $j < length( $1 ) - 1 ; $j++ ) {
          $deletionExtns++;
          print "" . $gapExtension . " " if ( $DEBUG );
          $score += $gapExtension;
          push @positionScores, $score;
        }
        print "\n" if ( $DEBUG );
        $i += length( $1 );
      }
      else {

        # Match/Mismatch
        #my $posScore = $$matScoresRef[ $qIdx ][ $sIdx ];
        my $posScore = $$matScoresRef[ $sIdx ][ $qIdx ];
        print "qIdx/sIdx = $qIdx/$sIdx Posscore = $posScore\n" if ( $DEBUG );
        $ungappedRawScore += $posScore;
        $score            += $posScore;
        $wellCharacterizedBases++
            if (    $sBase =~ /[ACGT]/
                 && $qBase =~ /[ACGT]/ );
        $matCounts[ $qIdx ]++ if ( $qIdx >= 0 );

        $transitions++
            if (    $qBase . $sBase eq "CT"
                 || $qBase . $sBase eq "TC"
                 || $qBase . $sBase eq "AG"
                 || $qBase . $sBase eq "GA" );
        $transversions++
            if (    $qBase . $sBase eq "GT"
                 || $qBase . $sBase eq "TG"
                 || $qBase . $sBase eq "GC"
                 || $qBase . $sBase eq "CG"
                 || $qBase . $sBase eq "CA"
                 || $qBase . $sBase eq "AC"
                 || $qBase . $sBase eq "AT"
                 || $qBase . $sBase eq "TA" );
        push @positionScores, $score;
        $i++;
      }
    }
  }

  #
  # Perform post-alignment xDrop calculation and record
  # where HSPs ( sub-alignments ) would have been reported.
  #
  my @xdropFragments = ();
  if ( exists $nameValueParams{'xDrop'} ) {
    my $xDrop               = $nameValueParams{'xDrop'};
    my $lastHighestScore    = 0;
    my $lastHighestScorePos = 0;
    my $startPoint          = 0;
    my $subtractScore       = 0;
    my $i;
    for ( $i = 0 ; $i <= $#positionScores ; $i++ ) {
      my $adjScore = $positionScores[ $i ] - $subtractScore;

      if ( $adjScore < 0 ) {
        $adjScore            = 0;
        $startPoint          = $i;
        $subtractScore       = $positionScores[ $i ];
        $lastHighestScore    = 0;
        $lastHighestScorePos = $i + 1;
      }

     #print "[$i] AdjScore = $adjScore, starting = $startPoint " .
     #      "lastHighestScore = $lastHighestScore pos = $lastHighestScorePos\n";
      if ( $adjScore >= $lastHighestScore ) {
        $lastHighestScore    = $positionScores[ $i ] - $subtractScore;
        $lastHighestScorePos = $i;
        next;
      }
      if ( ( $lastHighestScore - $adjScore ) > $xDrop ) {

        #print "Resetting score ( $startPoint, $lastHighestScorePos )\n";
        push @xdropFragments, ( $startPoint, $lastHighestScorePos );
        $subtractScore       = $positionScores[ $i ];
        $startPoint          = $i + 1;
        $lastHighestScore    = 0;
        $lastHighestScorePos = $i + 1;
      }
    }

    # Don't save trailing fragments
    if ( $lastHighestScorePos - $startPoint > 5 ) {
      push @xdropFragments, ( $startPoint, $lastHighestScorePos );
    }
  }

  my $qryBases = $this->getQueryEnd() - $this->getQueryStart() + 1;
  my $sbjBases = $this->getSubjEnd() - $this->getSubjStart() + 1;
  if ( $wellCharacterizedBases < 1 || $sbjBases < 1 || $qryBases < 1 ) {
    croak "$subroutine: Error corrupt search result!:\n"
        . $this->toStringFormatted( SearchResult::AlignWithQuerySeq ) . "\n";

  }
  my $p          = $transitions / $wellCharacterizedBases;
  my $q          = $transversions / $wellCharacterizedBases;
  my $percIns    = ( ( $insertionInits + $insertionExtns ) * 100 ) / $sbjBases;
  my $percDel    = ( ( $deletionInits + $deletionExtns ) * 100 ) / $qryBases;
  my $logOperand = ( ( 1 - ( 2 * $p ) - $q ) * ( 1 - ( 2 * $q ) )**0.5 );

  my $kimura = 100.00;
  if ( $logOperand > 0 ) {
    $kimura = ( abs( ( -0.5 * log( $logOperand ) ) ) * 100 );
  }

  if ( $DEBUG ) {
    print "#raw score = $score\n";
    print "#Kimura: $kimura\n";
    print "#Insertions(bp): "
        . ( $insertionInits + $insertionExtns )
        . " ( $percIns )\n";
    print "#Deletions(bp): "
        . ( $deletionInits + $deletionExtns )
        . " ( $percDel )\n";
    print "#ungapped raw score = $ungappedRawScore\n";
    print "#query bases = $qryBases\n";
    print "#subject bases = $sbjBases\n";
    print "#wellCharacterizedBases = $wellCharacterizedBases\n";
    print "#trans = $transitions transv = $transversions\n";
    print "#CpGSites = $CpGSites\n\n";
  }

  if ( exists $nameValueParams{'complexityAdjust'} ) {
    my $t_factor    = 0;
    my $t_sum       = 0;
    my $t_counts    = 0;
    my $n_letters   = 0;
    my $matFreqsRef = $matrix->getMatrixFreqsRef();
    for ( my $mCol = 0 ; $mCol < length( $matAlphabet ) ; $mCol++ ) {
      if ( defined $matCounts[ $mCol ] && $matCounts[ $mCol ] > 0 ) {
        if ( $$matFreqsRef[ $mCol ] > 0 && log( $$matFreqsRef[ $mCol ] ) != 0 )
        {
          my $count = $matCounts[ $mCol ];

          #print "# count $count\n";
          $t_factor += $count * log( $count );
          $t_sum    += $count * log( $$matFreqsRef[ $mCol ] );
          $t_counts += $count;
          $n_letters++;
        }
      }
    }

    print "#t_factor (unadj) = $t_factor\n" if ( $DEBUG );
    $t_factor -= $t_counts * log( $t_counts );
    $t_sum    -= $t_factor;

    my $matLambda = $matrix->getLambda();

    print "#Matrix lambda = $matLambda\n" if ( $DEBUG );
    print "#t_factor = $t_factor\n"       if ( $DEBUG );
    print "#t_sum = $t_sum\n"             if ( $DEBUG );
    my $adj_score = sprintf( "%0.0d", $score + $t_sum / $matLambda + .999 );

    $adj_score = 0 if ( !( $adj_score =~ /\d+/ ) || $adj_score < 0 );

    $score = $adj_score;
  }

  return (
           $score,           $kimura,                 $CpGSites,
           $percIns,         $percDel,                \@positionScores,
           \@xdropFragments, $wellCharacterizedBases, $transitions,
           $transversions
  );

}

##-------------------------------------------------------------------------##
## Private Methods
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my $csvString = $obj->_toCSVFormat();
##
##   TODO: Document
##
##-------------------------------------------------------------------------##
sub _toCSVFormat {
  my $obj           = shift;
  my $displayParams = shift;

  # AAGAA
  #   |
  # AACAA
  #
  # Encode as:  AAG/CAA
  #
  my $qryPos    = 0;
  my $sbjPos    = 0;
  my $insStart  = -1;
  my $delStart  = -1;
  my $indCount  = 0;
  my $indRec    = "";
  my $outSeq    = "";
  my $totIndels = 0;

  my $qrySeq = $obj->getQueryString();
  my $sbjSeq = $obj->getSubjString();
  while ( $qrySeq ne "" ) {
    $qrySeq =~ s/^(\S)//;
    my $qryChar = $1;
    $sbjSeq =~ s/^(\S)//;
    my $sbjChar = $1;
    if ( $qryChar eq "-" ) {
      if ( $insStart == -1 ) {
        $outSeq .= "+";
        $insStart = 1;
      }
      $outSeq .= $sbjChar;
      $sbjPos++;
    }
    elsif ( $sbjChar eq "-" ) {
      if ( $delStart == -1 ) {
        $outSeq .= "-";
        $delStart = 1;
      }
      $outSeq .= $qryChar;
      $qryPos++;
    }
    else {
      if ( $delStart != -1 ) {
        $outSeq .= "-";
        $delStart = -1;
      }
      elsif ( $insStart != -1 ) {
        $outSeq .= "+";
        $insStart = -1;
      }
      if ( $qryChar eq $sbjChar ) {
        $outSeq .= $qryChar;
      }
      else {
        $outSeq .= $qryChar . "/" . $sbjChar;
      }
      $qryPos++;
      $sbjPos++;
    }
  }

  my $cRec =
        $obj->getScore() . ","
      . $obj->getPctDiverge() . ","
      . $obj->getPctDelete() . ","
      . $obj->getPctInsert() . ","
      . $obj->getQueryName() . ","
      . $obj->getQueryStart() . ","
      . $obj->getQueryEnd() . ","
      . $obj->getQueryRemaining() . ","
      . $obj->getSubjName() . ","
      . $obj->getSubjType() . ","
      . $obj->getSubjStart() . ","
      . $obj->getSubjEnd() . ","
      . $obj->getSubjRemaining() . ",";

  if ( $obj->getOrientation =~ /C|c/ ) {
    $cRec .= "1,";
  }
  else {
    $cRec .= "0,";
  }

  $cRec .= $obj->getOverlap() . "," . $obj->getId() . "," . $outSeq;

  return $cRec;
}

##-------------------------------------------------------------------------##
## Use: my $cmString = $obj->_toCrossMatchFormat( $alignmentMode,
##                                                $displayParams );
##
## Build a Crossmatch-like alignment from the data in this object.
## Some extensions include swapping the alignment direction ( query vs
## subject centric ), highlighting subregions of the alignment ( using
## lower case ), or just displaying the header line and no alignment
## data.
##
##-------------------------------------------------------------------------##
sub _toCrossMatchFormat {
  my $obj           = shift;
  my $alignmentMode = shift;
  my $displayParams = shift;

  $alignmentMode = SearchResult::NoAlign
      if ( !defined $alignmentMode );
  croak $CLASS
      . "::toStringFormatted: Unknown alignment mode "
      . "( $alignmentMode )\n"
      if (    $alignmentMode != SearchResult::NoAlign
           && $alignmentMode != SearchResult::AlignWithQuerySeq
           && $alignmentMode != SearchResult::AlignWithSubjSeq
           && $alignmentMode != SearchResult::RangeHighlightedAlignment );

  my $retStr = "";
  my $sbjDir;
  my $sbjDirStr;
  my $alignIndex = 0;
  my $alignCol;

  #
  # Build annotation line
  #
  my ( $qryName ) = ( $obj->{'qryName'} =~ /(\S+).*/ );
  $retStr .=
        "$obj->{'score'} $obj->{percDiv} $obj->{percDel} "
      . "$obj->{'percIns'} $qryName $obj->{'qryBegin'} "
      . "$obj->{'qryEnd'} ($obj->{'qryLeft'}) ";
  my ( $sbjName ) = ( $obj->{'sbjName'} =~ /(\S+).*/ );
  if ( $obj->{sbjOrient} eq "C" ) {
    $retStr .=
        "C $sbjName ($obj->{sbjLeft}) " . "$obj->{sbjEnd} $obj->{sbjBegin}";
  }
  else {
    $retStr .=
        "$sbjName $obj->{sbjBegin} $obj->{sbjEnd} " . "($obj->{sbjLeft})";
  }
  if ( defined $obj->{id} ) {
    $retStr .= " $obj->{id}";
  }
  if ( defined $obj->{lineageId} ) {
    $retStr .= " $obj->{lineageId}";
  }
  if ( defined $obj->{overlap} ) {
    $retStr .= " $obj->{overlap}";
  }
  $retStr .= "\n";

  #
  # Build alignment data ( if requested )
  #
  if (    $alignmentMode != SearchResult::NoAlign
       && defined $obj->{'querySeq'}
       && $obj->{'querySeq'} ne "" )
  {
    my $qMasked;
    my $sMasked;
    my $insertions = 0;
    my $deletions  = 0;

    my $query   = $obj->{'querySeq'};
    my $subject = $obj->{'subjSeq'};

    #
    # Higlight alignment ( if requested )
    #
    if (
         $alignmentMode == SearchResult::RangeHighlightedAlignment
         && (    defined $displayParams->{'qryRangeList'}
              || defined $displayParams->{'sbjRangeList'} )
        )
    {
      if ( defined $displayParams->{'qryRangeList'} ) {
        my $rl   = $displayParams->{'qryRangeList'};
        my $qPos = $obj->{'qryBegin'};
        if ( length( $query ) != length( $subject ) ) {
          print "Bad:\n$query\n$subject\n\n";
        }
        for ( my $i = 0 ; $i < length( $query ) ; $i++ ) {
          last if ( !@{$rl} );
          if ( $qPos >= $rl->[ 0 ] && $qPos <= $rl->[ 1 ] ) {

            # highlight
            substr( $query,   $i, 1 ) = lc( substr( $query,   $i, 1 ) );
            substr( $subject, $i, 1 ) = lc( substr( $subject, $i, 1 ) )
                if ( $i < length( $subject ) );
          }
          elsif ( $qPos > $rl->[ 1 ] ) {

            # Remove range
            shift @{$rl};
            shift @{$rl};
          }
          next if ( substr( $query, $i, 1 ) eq "-" );
          $qPos++;
        }
      }
      else {
        my $rl   = $displayParams->{'sbjRangeList'};
        my $sPos = $obj->{'sbjBegin'};
        for ( my $i = 0 ; $i <= length( $subject ) ; $i++ ) {
          last if ( !@{$rl} );
          next if ( substr( $subject, $i, 1 ) eq "-" );
          if ( $sPos >= $rl->[ 0 ] && $sPos <= $rl->[ 1 ] ) {

            # highlight ( currently by lower casing the typically UC sequence )
            substr( $query,   $i, 1 ) = lc( substr( $query,   $i, 1 ) );
            substr( $subject, $i, 1 ) = lc( substr( $subject, $i, 1 ) );
          }
          elsif ( $sPos > $rl->[ 1 ] ) {

            # Remove range
            shift @{$rl};
            shift @{$rl};
          }
          $sPos++;
        }
      }
    }

    $retStr .= "\n";

    if (    $obj->{'sbjOrient'} eq "C"
         && $alignmentMode == SearchResult::AlignWithSubjSeq )
    {
      $query = reverse $query;
      $query =~
          tr/ACGTYRMKHBVDacgtyrmkhbvd/TGCARYKMDVBHtgcarykmdvbh/;    # complement
      $subject = reverse $subject;
      $subject =~
          tr/ACGTYRMKHBVDacgtyrmkhbvd/TGCARYKMDVBHtgcarykmdvbh/;    # complement
    }

    my $qStart = $obj->{'qryBegin'};
    if (    $obj->{'sbjOrient'} eq "C"
         && $alignmentMode == SearchResult::AlignWithSubjSeq )
    {
      $qStart = $obj->{'qryEnd'};
    }
    my $qEnd   = 0;
    my $sStart = $obj->{'sbjBegin'};
    if (    $obj->{'sbjOrient'} eq "C"
         && $alignmentMode == SearchResult::AlignWithQuerySeq )
    {
      $sStart = $obj->{'sbjEnd'};
    }
    my $sEnd = 0;
    while ( $query ) {
      $query =~ s/^(.{1,50})//;
      my $qSeq = $1;
      $subject =~ s/^(.{1,50})//;
      my $sSeq = $1;

      $insertions = ( $qSeq =~ tr/-/-/ );
      $deletions  = ( $sSeq =~ tr/-/-/ );

      if ( $sEnd > 0 ) {
        my $qIncr = 0;
        my $sIncr = 0;
        $qIncr = 1 if ( length( $qSeq ) > $insertions );
        $sIncr = 1 if ( length( $sSeq ) > $deletions );
        if ( $obj->{'sbjOrient'} eq "C" ) {
          if ( $alignmentMode == SearchResult::AlignWithSubjSeq ) {
            $qStart = $qEnd - $qIncr;
            $sStart = $sEnd + $sIncr;
          }
          else {
            $qStart = $qEnd + $qIncr;
            $sStart = $sEnd - $sIncr;
          }
        }
        else {
          $qStart = $qEnd + $qIncr;
          $sStart = $sEnd + $sIncr;
        }
      }

      # Indicate orientation in alignment data by placing a "C" in front
      # of sequences which have been reverse complemented
      if (    $obj->{'sbjOrient'} eq "C"
           && $alignmentMode == SearchResult::AlignWithSubjSeq )
      {
        $qEnd = $qStart - length( $qSeq ) + 1 + $insertions;
        $retStr .= "C ";
      }
      else {
        $qEnd = $qStart + length( $qSeq ) - 1 - $insertions;
        $retStr .= "  ";
      }
      $qEnd = $qStart if ( length( $qSeq ) == $insertions );

      # Up to 13 characters are allowed from the QueryName/SubjName
      $retStr .= substr( $obj->{'qryName'}, 0, 13 )
          . " " x (
           13 - (
             length( $obj->{'qryName'} ) < 13 ? length( $obj->{'qryName'} ) : 13
           )
          );

      # Required Intermediate Space
      $retStr .= " ";

      # Up to 10 characters for the position, followed by the sequence and
      # end positions
      $retStr .= " " x ( 10 - length( $qStart ) ) . $qStart . " $qSeq $qEnd\n";

      # The modification codes
      $retStr .= " " x 27;
      for ( my $j = 0 ; $j < length( $qSeq ) ; $j++ ) {
        my $qChar = substr( $qSeq, $j, 1 );
        my $sChar = substr( $sSeq, $j, 1 );
        if ( $qChar eq $sChar ) {
          $retStr .= " ";
        }
        elsif ( $qChar eq "-" || $sChar eq "-" ) {
          $retStr .= "-";
        }
        elsif (    ( $qChar . $sChar ) =~ /CT/i
                || ( $qChar . $sChar ) =~ /TC/i
                || ( $qChar . $sChar ) =~ /AG/i
                || ( $qChar . $sChar ) =~ /GA/i )
        {
          $retStr .= "i";
        }
        elsif (    ( $qChar . $sChar ) =~ /GT/i
                || ( $qChar . $sChar ) =~ /TG/i
                || ( $qChar . $sChar ) =~ /GC/i
                || ( $qChar . $sChar ) =~ /CG/i
                || ( $qChar . $sChar ) =~ /CA/i
                || ( $qChar . $sChar ) =~ /AC/i
                || ( $qChar . $sChar ) =~ /AT/i
                || ( $qChar . $sChar ) =~ /TA/i )
        {
          $retStr .= "v";
        }
        elsif (    ( $qChar =~ /[BDHVRYKMSWNX]/i )
                || ( $sChar =~ /[BDHVRYKMSWNX]/i ) )
        {
          $retStr .= "?";
        }
        else {
          $retStr .= " ";
        }
      }
      $retStr .= "\n";

      # Subject orientation/name
      if (    $obj->{'sbjOrient'} eq "C"
           && $alignmentMode == SearchResult::AlignWithQuerySeq )
      {
        $retStr .= "C ";
        $sEnd = $sStart - length( $sSeq ) + 1 + $deletions;
      }
      else {
        $retStr .= "  ";
        $sEnd = $sStart + length( $sSeq ) - 1 - $deletions;
      }
      $retStr .= substr( $obj->{'sbjName'}, 0, 13 )
          . " " x (
           13 - (
             length( $obj->{'sbjName'} ) < 13 ? length( $obj->{'sbjName'} ) : 13
           )
          );

      # Required Intermediate Space
      $retStr .= " ";

      $retStr .= " " x ( 10 - length( $sStart ) ) . $sStart . " $sSeq $sEnd\n";
      $retStr .= "\n";

    }    # while ( $query )

    if ( defined $obj->getMatrixName() ) {
      $retStr .= "Matrix = " . $obj->getMatrixName() . "\n";
    }
    else {
      $retStr .= "Matrix = Unknown\n";
    }
    $retStr .= "Transitions / transversions = ";
    my ( $mismatches, $transitions, $transversions, $numGaps, $totGapLen ) =
        $obj->_getAlignmentStats();
    if ( defined $transitions ) {
      if ( $transversions > 0 ) {
        $retStr .= sprintf( "%0.2f", ( $transitions / $transversions ) );
      }
      else {
        $retStr .= "0.0";
      }
      $retStr .= " ($transitions / $transversions)\n";
      if ( ( $obj->getQueryEnd() - $obj->getQueryStart() ) > 0 ) {
        $retStr .= "Gap_init rate = "
            . sprintf( "%0.2f",
                       $numGaps /
                           ( $obj->getQueryEnd() - $obj->getQueryStart() ) )
            . " ($numGaps / "
            . ( $obj->getQueryEnd() - $obj->getQueryStart() ) . ")";
      }
      else {
        $retStr .= "Gap_init rate = 0.0 ( $numGaps / 0 )";
      }
      if ( $numGaps > 0 ) {
        $retStr .=
              ", avg. gap size = "
            . sprintf( "%0.2f", ( $totGapLen / $numGaps ) )
            . " ($totGapLen / $numGaps)\n\n";
      }
      else {
        $retStr .= ", avg. gap size = 0.0 (0 / 0)\n\n";
      }

    }
    else {
      $retStr .= "Transitions / transversions = Unknown\n";
      $retStr .= "Gap_init rate = Unknown, avg. gap size = Unknown\n\n";
    }
  }

  return $retStr;
}

##-------------------------------------------------------------------------##
## Use: my $outString = $obj->_toOUTFileFormat();
##
##   TODO: Document
##
##-------------------------------------------------------------------------##
sub _toOUTFileFormat {
  my $obj = shift;

  my $outStr    = "";
  my $orient    = "+";
  my $sbjCoord1 = $obj->{sbjBegin};
  my $sbjCoord2 = $obj->{sbjEnd};
  my $sbjCoord3 = "(" . $obj->{sbjLeft} . ")";
  if ( $obj->{sbjOrient} eq "C" ) {
    $orient    = "C";
    $sbjCoord1 = "(" . $obj->{sbjLeft} . ")";
    $sbjCoord2 = $obj->{sbjEnd};
    $sbjCoord3 = $obj->{sbjBegin};
  }
  $outStr = sprintf(
                     "%6d %4.1f %4.1f %4.1f %-17s %8d %8d %8s %1s %-15s %-15s "
                         . "%7s %7s %7s %-5s %3s %3s\n",
                     $obj->{score},   $obj->{percDiv},
                     $obj->{percDel}, $obj->{percIns},
                     $obj->{qryName}, $obj->{qryBegin},
                     $obj->{qryEnd},  "(" . $obj->{qryLeft} . ")",
                     $orient,         $obj->{sbjName},
                     $obj->{sbjType}, $sbjCoord1,
                     $sbjCoord2,      $sbjCoord3,
                     $obj->{id},      $obj->{lineageId},
                     $obj->{overlap}
  );

  return $outStr;

}

##-------------------------------------------------------------------------##
## Use: my ($mismatches, $transitions, $transversion, $numGaps, $totGapLen) =
##                                                 _getAlignmentStats();
##
##   Given an alignment calculate the number of mismatches, transitions
##   transversions, gap count, average gap size etc.
##
##-------------------------------------------------------------------------##
sub _getAlignmentStats {
  my $obj = shift;

  my $querySeq = $obj->getQueryString();
  my $subjSeq  = $obj->getSubjString();

  return if ( !defined $querySeq || $querySeq eq "" );

  my $transitions   = 0;
  my $transversions = 0;
  my $gaps          = 0;
  my $totGapLen     = 0;
  my $mismatches    = 0;
  my $inGap         = 0;

  while ( $querySeq ne "" ) {
    $querySeq =~ s/^(\S)//;
    my $qryChar = $1;
    $subjSeq =~ s/^(\S)//;
    my $sbjChar = $1;
    if ( $qryChar eq $sbjChar ) {
      $inGap = 0;
      next;
    }
    if ( $qryChar eq "-" || $sbjChar eq "-" ) {
      if ( $inGap == 0 ) {
        $gaps++;
        $inGap = 1;
      }
      $totGapLen++;
    }
    else {
      $inGap = 0;
      my $basePair = uc( $qryChar . $sbjChar );
      if ( $basePair =~ /CT|TC|GA|AG/ ) {
        $transitions++;
      }
      elsif ( $basePair =~ /GT|TG|TA|AT|CA|AC|CG|GC/ ) {
        $transversions++;
      }
      $mismatches++;
    }
  }

  return ( $mismatches, $transitions, $transversions, $gaps, $totGapLen );
}

##-------------------------------------------------------------------------##
## Use: my _ucFirst( $string );
##
##   Uppercases the first character in a string and returns it.
##
##-------------------------------------------------------------------------##
sub _ucFirst {
  my $string = shift;

  if ( defined $string && $string ne "" ) {
    substr( $string, 0, 1 ) = uc( substr( $string, 0, 1 ) );
  }
  return $string;
}

##-------------------------------------------------------------------------##
## Serialization & Debug Routines
##-------------------------------------------------------------------------##

##-------------------------------------------------------------------------##
## Use: my $string = toString([$this]);
##
##      $this         : Normally passed implicitly
##
##  Returns
##
##      Uses the Data::Dumper to create a printable reprentation
##      of a data structure.  In this case the object data itself.
##
##-------------------------------------------------------------------------##
sub toString {
  my $this = shift;
  my $data_dumper = new Data::Dumper( [ $this ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  return $data_dumper->Dump();
}

##-------------------------------------------------------------------------##
## Use: my serializeOUT( $filename );
##
##	  $filename	: A filename to be created
##
##  Returns
##
##	Uses the Data::Dumper module to save out the data
##	structure as a text file.  This text file can be
##	read back into an object of this type.
##
##-------------------------------------------------------------------------##
sub serializeOUT {
  my $this     = shift;
  my $fileName = shift;

  my $data_dumper = new Data::Dumper( [ $this ] );
  $data_dumper->Purity( 1 )->Terse( 1 )->Deepcopy( 1 );
  open OUT, ">$fileName";
  print OUT $data_dumper->Dump();
  close OUT;
}

##-------------------------------------------------------------------------##
## Use: my serializeIN( $filename );
##
##	$filename	: A filename containing a serialized object
##
##  Returns
##
##	Uses the Data::Dumper module to read in data
##	from a serialized PERL object or data structure.
##
##-------------------------------------------------------------------------##
sub serializeIN {
  my $this         = shift;
  my $fileName     = shift;
  my $fileContents = "";
  my $oldSep       = $/;
  undef $/;
  my $in;
  open $in, "$fileName";
  $fileContents = <$in>;
  $/            = $oldSep;
  close $in;
  return eval( $fileContents );
}

1;
