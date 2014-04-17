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
#     Revision 1.21  2013/11/06 19:01:26  rhubley
#     Cleanup before a distribution
#
#
###############################################################################
#
# To Do:
#

=head1 NAME

calcDivergenceFromAlign.pl - Recalculate/Summarize the divergences in an align file.

=head1 SYNOPSIS

  calcDivergenceFromAlign.pl [-version] [-s <summary_file>] [-noCpGMod]
                             [-a <new_align_file>] *.align[.gz]

=head1 DESCRIPTION

  A utility script to calculate a new divergence measure on the
  RM alignment files.  Currently we only calculate the Kimura 2-Parameter
  divergence metric.  

  Treat "CG" dinucleotide sites in the consensus sequence as follows:
  Two transition mutations are counted as a single transition, one
  transition is counted as 1/10 of a standard transition, and 
  transversions are counted normally (as the would outside of a CpG
  site).  This modification to the Kimura 2 parameter model accounts
  for the extremely high rate of mutations in at a CpG locus.  


The options are:

=over 4

=item -version

Displays the version of the program

=item -noCpGMod

Do not modify the transition counts at CpG sites. 

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2013 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use Getopt::Long;
use Data::Dumper;
use FileHandle;

## TODO: Remove this
use lib "/home/rhubley/projects/RepeatMasker";
use FindBin;
use lib $FindBin::RealBin;
use lib "$FindBin::Bin/..";
use SearchResult;
use CrossmatchSearchEngine;

#
# Version
#
#  This is a neat trick.  CVS allows you to tag
#  files in a repository ( i.e. cvs tag "2003/12/03" ).
#  If you check out that release into a new
#  directory with "cvs co -r "2003/12/03" it will
#  place this string into the $Name: open-4-0-5 $ space below
#  automatically.  This will help us keep track
#  of which release we are using.  If we simply
#  check out the code as "cvs co Program" the
#  $Name: open-4-0-5 $ macro will be blank so we should default
#  to what the ID tag for this file contains.
#
my $CVSNameTag = '$Name: open-4-0-5 $';
my $CVSIdTag   =
    '$Id: calcDivergenceFromAlign.pl,v 1.21 2013/11/06 19:01:26 rhubley Exp $';
my $Version = $CVSNameTag;
$Version = $CVSIdTag if ( $Version eq "" );

#
# Magic numbers/constants here
#  ie. my $PI = 3.14159;
#
my $DEBUG = 0;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
                    '-version',    # print out the version and exit
                    '-noCpGMod',
                    '-a=s',
                    '-s=s'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n\n";
  exec "pod2text $0";
  exit;
}

if ( $options{'version'} ) {
  print "$Version\n";
  exit;
}

if ( !$options{'s'} && !$options{'a'} ) {
  print
"\n\nError: One or more of the options '-a' or '-s' must be supplied!\n\n";
  usage();
}

my $alignFile  = $ARGV[ 0 ];
my $maxDiv     = 70;
my $cntAlign   = 0;
my %repeatMuts = ();
my %len        = ();
my %seqSizes   = ();

my $searchResultsFH = new FileHandle;

if ( $alignFile =~ /.+\.gz/ ) {
  open $searchResultsFH, "gunzip -c $alignFile|"
      or die
      "RepeatLandscape: Could not open gunzip for reading $alignFile: $!\n";
}
else {
  open $searchResultsFH, "<$alignFile"
      or die "RepeatLandscape: Could not open $alignFile for reading: $!\n";
}

if ( $options{'s'} ) {
  open SOUT, ">$options{'s'}"
      or die "Error: Could not open $options{'s'} for writing!\n";
  if ( $options{'noCpGMod'} ) {
    print SOUT "Jukes/Cantor and Kimura subsitution levels\n";
    print SOUT "==========================================\n";
  }
  else {
    print SOUT
        "Jukes/Cantor and Kimura subsitution levels adjusted for CpG sites\n";
    print SOUT
        "=================================================================\n";
  }
  print SOUT "File: " . $alignFile . "\n";
}

#
# Process the alignment file
#
if ( $options{'a'} ) {
  open COUT, ">$options{'a'}"
      or die "Could not open $options{'a'} for writing!\n";
  CrossmatchSearchEngine::parseOutput( searchOutput => $searchResultsFH,
                                       callback => \&processAlignmentWithOutput
  );
  close COUT;
}
else {
  CrossmatchSearchEngine::parseOutput( searchOutput => $searchResultsFH,
                                       callback     => \&processAlignment );
}

if ( $options{'s'} ) {
  my $genomeSize = 0;
  foreach my $seq ( keys( %seqSizes ) ) {
    $genomeSize += $seqSizes{$seq};
  }
  print SOUT
"Genome Size = $genomeSize bp ( Calculated from alignment data -- may be underestimate )\n\n";
  print SOUT "Class	Repeat	Jukes%	Kimura%\n";
  print SOUT "-----     ------  -----   -------\n";
  foreach my $class ( sort keys %repeatMuts ) {
    foreach my $id ( sort keys %{ $repeatMuts{$class} } ) {
      my $kimura = 100.0;

      # Determine which method was used to obtain divergence
      if ( $repeatMuts{$class}->{$id}->{'transitions'} ) {
        my $p =
            $repeatMuts{$class}->{$id}->{'transitions'} /
            $repeatMuts{$class}->{$id}->{'length'};
        my $q =
            $repeatMuts{$class}->{$id}->{'transversions'} /
            $repeatMuts{$class}->{$id}->{'length'};
        my $logOperand = ( ( 1 - ( 2 * $p ) - $q ) * ( 1 - ( 2 * $q ) )**0.5 );
        if ( $logOperand > 0 ) {
          $kimura = ( abs( ( -0.5 * log( $logOperand ) ) ) * 100 );
        }
        $kimura = sprintf( "%4.2f", $kimura );
      }
      else {
        $kimura = sprintf( "%4.2f",
                           $repeatMuts{$class}->{$id}->{'sumdiv'} /
                               $repeatMuts{$class}->{$id}->{'count'} );
      }
      $kimura = $maxDiv if ( $kimura > $maxDiv );
      print SOUT "$class\t$id\t--\t$kimura\n";
    }
  }
  print SOUT "\n\n";
  print SOUT "Coverage for each repeat class and divergence (Kimura)\n";
  print SOUT "Div ";
  foreach my $class ( sort keys %repeatMuts ) {
    print SOUT "$class ";
  }
  print SOUT "\n";

  my $j = 0;
  while ( $j <= $maxDiv ) {
    print SOUT "$j ";
    foreach my $class ( sort keys %repeatMuts ) {
      my $label = "$class $j";

      # counting is pretty useless without linked fragments
      #  -- can check for modern *.align files that contain the linkage info
      #$cnt{$cat} = 0 unless $cnt{$cat};
      $len{$label} = 0 unless $len{$label};

      #print "$cnt{$cat} $len{$cat} ";
      print SOUT "$len{$label} ";
    }
    print SOUT "\n";
    ++$j;
  }
  close SOUT;
}

exit;

######################## S U B R O U T I N E S ############################

##-------------------------------------------------------------------------##
## Use: my _privateMethod( $parameter => value );
##
##      $parameter       : A parameter to the method
##
##  Returns
##      Methods with the prefix "_" are conventionally considered
##      private.  This bit-o-documentation is not formatted to
##      print out when perldoc is run on this file.
##
##-------------------------------------------------------------------------##
sub processAlignmentWithOutput {
  my $result = shift;

  return if ( !$result );

  my $subjName = $result->getSubjName();
  my $hitname;
  my $class;
  if ( $subjName =~ /(\S+)\#(\S+)/ ) {
    $hitname = $1;
    $class   = $2;
  }
  else {
    $hitname = $subjName;
    $class   = $result->getSubjType();
  }

  my $seqName = $result->getQueryName();
  my $seqSize = $result->getQueryEnd() + $result->getQueryRemaining();
  $seqSizes{$seqName} = $seqSize;

  # TODO: double check sizes are consistent?

  if ( $class =~ /Simple|Low_complexity/ ) {
    print COUT ""
        . $result->toStringFormatted( SearchResult::AlignWithQuerySeq ) . "\n";
    return;
  }

  print STDERR "." if ( $cntAlign++ % 1000 == 0 );

  my ( $div, $transi, $transv, $wellCharBases, $numCpGs );

  ##
  ## TODO: There are two methods encoded here.  The first
  ##       uses the "Kimura .. = #.#" lines from a modern
  ##       ( >= 4.0.4 ) *.align file and averages them by
  ##       count of annotations.  The second method counts
  ##       all transversions/transitions over all aligned
  ##       bases in a class/id and reports that as the
  ##       divergence.  These methods are not the same.
  ##

  # Obtain divergence from modern *.align files directly
  $div = $result->getPctKimuraDiverge();
  my $alen = $result->getQueryEnd() - $result->getQueryStart() + 1;
  $wellCharBases = $alen - int( $alen * ( $result->getPctInsert() / 100 ) );

  if ( $div eq "" ) {

    # Calculate divergence on the fly
    ( $div, $transi, $transv, $wellCharBases, $numCpGs ) =
        $result->calcKimuraDivergence( divCpGMod => 1 );
    $result->setPctKimuraDiverge( sprintf( "%4.2f", $div ) );
    $repeatMuts{$class}->{$hitname}->{'transitions'}   += $transi;
    $repeatMuts{$class}->{$hitname}->{'transversions'} += $transv;
  }

  $repeatMuts{$class}->{$hitname}->{'sumdiv'} += $div;
  $repeatMuts{$class}->{$hitname}->{'length'} += $wellCharBases;
  $repeatMuts{$class}->{$hitname}->{'count'}++;
  $div = int( $div );
  my $label = "$class $div";
  $len{$label} += $wellCharBases;
  print COUT ""
      . $result->toStringFormatted( SearchResult::AlignWithQuerySeq ) . "\n";
}

##-------------------------------------------------------------------------##
## Use: my _privateMethod( $parameter => value );
##
##      $parameter       : A parameter to the method
##
##  Returns
##      Methods with the prefix "_" are conventionally considered
##      private.  This bit-o-documentation is not formatted to
##      print out when perldoc is run on this file.
##
##-------------------------------------------------------------------------##
sub processAlignment {
  my $result = shift;

  return if ( !$result );

  my $subjName = $result->getSubjName();
  my $hitname;
  my $class;
  if ( $subjName =~ /(\S+)\#(\S+)/ ) {
    $hitname = $1;
    $class   = $2;
  }
  else {
    $hitname = $subjName;
    $class   = $result->getSubjType();
  }

  my $seqName = $result->getQueryName();
  my $seqSize = $result->getQueryEnd() + $result->getQueryRemaining();
  $seqSizes{$seqName} = $seqSize;

  # TODO: double check sizes are consistent?

  return if ( $class =~ /Simple|Low_complexity/ );

  print STDERR "." if ( $cntAlign++ % 1000 == 0 );

  my ( $div, $transi, $transv, $wellCharBases, $numCpGs );

  ##
  ## TODO: There are two methods encoded here.  The first
  ##       uses the "Kimura .. = #.#" lines from a modern
  ##       ( >= 4.0.4 ) *.align file and averages them by
  ##       count of annotations.  The second method counts
  ##       all transversions/transitions over all aligned
  ##       bases in a class/id and reports that as the
  ##       divergence.  These methods are not the same.
  ##

  # Obtain divergence from modern *.align files directly
  $div = $result->getPctKimuraDiverge();
  my $alen = $result->getQueryEnd() - $result->getQueryStart() + 1;
  $wellCharBases = $alen - int( $alen * ( $result->getPctInsert() / 100 ) );

  if ( $div eq "" ) {

    # Calculate divergence on the fly
    ( $div, $transi, $transv, $wellCharBases, $numCpGs ) =
        $result->calcKimuraDivergence( divCpGMod => 1 );
    $result->setPctKimuraDiverge( sprintf( "%4.2f", $div ) );
    $repeatMuts{$class}->{$hitname}->{'transitions'}   += $transi;
    $repeatMuts{$class}->{$hitname}->{'transversions'} += $transv;
  }

  $repeatMuts{$class}->{$hitname}->{'sumdiv'} += $div;
  $repeatMuts{$class}->{$hitname}->{'length'} += $wellCharBases;
  $repeatMuts{$class}->{$hitname}->{'count'}++;
  $div = int( $div );
  my $label = "$class $div";
  $len{$label} += $wellCharBases;
}

1;
