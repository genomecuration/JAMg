#!/usr/bin/env perl

=pod

=head1 NAME

#############################################################
# optimize_augustus.pl based on Mario Stanke's script (23.04.2007) 
#
# Mario Stanke, 23.04.2007
#############################################################

=head1 USAGE 

Mandatory parameters:

 --species                prefix of the species name
 --optimize               genbank file for training with bona fide gene structures

Important optional parameters:

 --metapars               File with the names and their ranges of the meta parameters that are subject to optimization (default: generic_metapars.cfg)
 --cpus                   The number of CPUs to use (default: 1)

Useful optional parameters:

 --rounds                 The number of rounds to run the optimization for (default: 5)
 --onlytrain              an optional genbank file that is used in addition to train.gb but only for etrain not for intermediate evaluation of accuracy. These genes may e.g. be incomplete.
 --kfold                  Make a k-fold cross validation (default: 8)
 --config_path            Specify the config directory d if not set as environment variable
 --aug_exec_dir           Path to augustus and etraining executable. If not specified it must be in $PATH environment variable
 --onlyutr                Use this option, if the exon, intron and intergenic models need not be trained. (default: 0)
 
Other optional parameters:

 --pstep                  For integer and floating parameters start with p tests equidistributed in the allowed range of values (default: 5)
 --min_coding             Minimum coding length
 --UTR                    Turn untranslated region model on for training and prediction
 --notrain                Use this option, if the parameters to optimize do not affect training. The training step (etraining) is omitted completely. (default: 0)

Others (expert):

 --genemodel              Augustus genemodel parameter
 --trans_matrix       Optimize the transition matrix file s. s must be the transition file used. e.g. ../species/nt/generic/generic_trans_shadow_partial.pbl
 --matrix_constraints     A file with try list, normed list and bindings


=cut

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use Time::localtime;
use IO::File;
use threads;
use threads::shared;
use FindBin;
use lib ( $FindBin::RealBin . '/../PerlLib' );
use Thread_helper;

my (
     $species,      $metapars,           $rounds,      $onlytrain,
     $kfold,        $pstep,              $config_path, $cpus,
     $trans_matrix, $matrix_constraints, $utr,         $exec_dir,
     $onlyutr,     $notrain,            $trans_table, $genemodel,
     $min_coding,   $output_directory,   $verbose,     $optimize_gb,
     $onlytrain_gb
);

$rounds      = 5;
$kfold       = 8;
$pstep       = 5;
$cpus        = 1;
$exec_dir    = $ENV{'AUGUSTUS_PATH'} . '/bin' if $ENV{'AUGUSTUS_PATH'};
$config_path = $ENV{'AUGUSTUS_PATH'} . '/config' if $ENV{'AUGUSTUS_PATH'};
$config_path = $ENV{'AUGUSTUS_CONFIG_PATH'} if $ENV{'AUGUSTUS_CONFIG_PATH'};

my ( $common_parameters, $modelrestrict ) = ('','');

&GetOptions(
             'species:s'              => \$species,
             'optimize:s'             => \$optimize_gb,
             'metapars:s'             => \$metapars,
             'rounds:i'               => \$rounds,
             'onlytrain:s'            => \$onlytrain_gb,
             'kfold:i'                => \$kfold,
             'pstep:i'                => \$pstep,
             'AUGUSTUS_CONFIG_PATH:s' => \$config_path,
             'cpus|threads:i'         => \$cpus,
             'trans_matrix:s'     => \$trans_matrix,
             'matrix_constraints:s'   => \$matrix_constraints,
             'UTR'                    => \$utr,
             'aug_exec_dir:s'         => \$exec_dir,
             'onlyutr'               => \$onlyutr,
             'notrain'                => \$notrain,
             'translation_table:i'    => \$trans_table,
             'genemodel:s'            => \$genemodel,
             'min_coding_len:i'       => \$min_coding,
             'output:s'               => \$output_directory,
             'verbose'                => \$verbose
);

#globals
$SIG{INT} = \&got_interrupt_signal;
my %storedsnsp;    # hash with the stored sn and sp array references
my $be_silent =
" --/augustus/verbosity=0 --/ExonModel/verbosity=0 --/IGenicModel/verbosity=0 --/IntronModel/verbosity=0 --/UtrModel/verbosity=0 --/genbank/verbosity=0 ";
&check_options();
my ( $augustus_exec, $etrain_exec ) = &check_program( 'augustus', 'etraining' );



##############################################################
# Read in the meta parameters
##############################################################

my @metastartvalues = ();
my @metaparnames    = ();
my @metaparranges   = ();

###############################################
# open the file with the parameters to optimize
#
# either metapars.cfg or opt_transition_matrix
###############################################
my $metaparsfilename;
my $n = 0;    # number of states
my @trans;    # transition matrix, array of array references

if ( !$trans_matrix ) {    # optimize meta parameters
 if ($metapars) {
  $metaparsfilename = $metapars;
 }
 else {
  $metaparsfilename = $config_path . "species/generic/generic_metapars.cfg";
 }
 open( META, $metaparsfilename ) or die("Could not open $metaparsfilename\n");
 print 
"Reading in the meta parameters used for optimization from $metaparsfilename...\n";
 my $orig_sep = $/;
 $/ = "\n";
 while (<META>) {
  my ( $parname, $range );
  if (/^\s*\#/) {
   next;
  }
  if (/^\s*(\S+)\s+(.*)\s*/) {
   $parname = $1;
   $range   = $2;
   if ( idx( \@metaparnames, $1 ) != -1 ) {
    die("Meta parameter $1 occurs twice in $metaparsfilename.");
   }
   else {
    if ( $range =~ /"([^"]+)"-"([^"]+)"/ ) {
     push @metaparnames, $parname;
     push @metaparranges, [ 'intrange', $1, $2 ];
    }
    elsif ( $range =~ /"([^"]+)"_"([^"]+)"/ ) {
     push @metaparnames, $parname;
     push @metaparranges, [ 'floatrange', $1, $2 ];
    }
    else {
     my @tokens = split /\s+/, $2;
     my @list = ();
     foreach (@tokens) {
      s/^"(.*)"$/$1/;
      push @list, $_;
     }
     push @metaparnames, $parname;
     push @metaparranges, [ 'list', @list ];
    }
   }
  }
 }
 $/ = $orig_sep;
}
else {    # read in transition matrix for optimization
 open( TRANS, $trans_matrix )
   or die("Could not open transition matrix file $trans_matrix");
 print "Reading in the transition matrix...\n";
 my $orig_sep = $/;
 $/ = "\n";
 while (<TRANS>) {
  my ( $from, $to, $prob );
  if (/^\s*\#/) {
   next;
  }
  if ( $n == 0 && /(\d+)/ ) {
   $n = $1;
   print "Transition matrix has dimension ${n}x${n}.\n";
  }
  if (/^\s*(\d+)\s+(\d+)\s*(\S+)/) {
   $from = $1;
   $to   = $2;
   $prob = $3;

   #print "trans[$from][$to]=$prob\n";
   if ( $from < 0 || $from >= $n || $to < 0 || $to >= $n ) {
    print
"State of transition matrix out of bounds ($n) for transition $from->$to:$prob\n";
   }
   if ( $prob < 0 ) {
    print "Error: negative probability in transition $from->$to:$prob\n";
   }
   if ( !defined( $trans[$from] ) ) {
    $trans[$from] = [];
   }
   $trans[$from][$to] = $prob;
  }
 }
 $/ = $orig_sep;
}

#printmetaranges(\@metaparnames, \@metaparranges);

# open species_parameters.cfg
my ( @spcfilelines, @transfilelines );
my $speciesdir           = $config_path . "species/$species/";
my $species_cfg_filename = $speciesdir . $species . "_parameters.cfg";
if ( !$trans_matrix ) {

 # make a copy of the original parameter file
 my $y = 1;
 while ( $y < 20
         && sysopen( ORIG, "$species_cfg_filename.orig$y", O_WRONLY | O_EXCL ) )
 {
  $y++;
 }
 if ( $y < 20 ) {
  close(ORIG);
  system("cp $species_cfg_filename $species_cfg_filename.orig$y");
 }
 else {
  die("Too many $species_cfg_filename.orig copies. Please delete some.");
 }

 if ( -e "$species_cfg_filename" ) {
  open( SPCCFG, "<$species_cfg_filename" )
    or die("Could not open $species_cfg_filename");
 }
 else { die "File $species_cfg_filename does not seem to exist!\n"; }
 print "Reading in the starting meta parameters from $species_cfg_filename...\n";
 my $orig_sep = $/;
 $/            = "\n";
 @spcfilelines = <SPCCFG>;
 close(SPCCFG);
 foreach (@spcfilelines) {
  my ( $parname, $value );
  if ( /^\s*\#.*/ || /^\s*$/ ) {    # skip comment lines
   next;
  }
  if (/^\s*(\S+)\s+(\S*)\s*/) {
   $parname = $1;
   $value   = $2;
   my $index = idx( \@metaparnames, $1 );
   if ( $index != -1 ) {
    $metastartvalues[$index] = $value;
   }
  }
 }
 $/ = $orig_sep;

 for ( my $i = 0 ; $i <= $#metaparnames ; $i++ ) {
  if ( !defined $metastartvalues[$i] ) {
   die(
"No start value for parameter $metaparnames[$i] found in file $species_cfg_filename.\n\
Maybe you misspelled this parameter in $metaparsfilename.\n" );
  }
 }
}
else {

 # make a copy of the original transition matrix file
 my $y = 1;
 while ( $y < 40
         && sysopen( ORIG, "$trans_matrix.orig$y", O_WRONLY | O_EXCL ) )
 {
  $y++;
 }
 if ( $y < 40 ) {
  close(ORIG);
  system("cp $trans_matrix $trans_matrix.orig$y");
 }
 else {
  die("Too many $trans_matrix.orig copies. Please delete some.");
 }

 open( TRANS, $trans_matrix ) or die("Could not open $trans_matrix");
 my $orig_sep = $/;
 $/              = "\n";
 @transfilelines = <TRANS>;
 close(TRANS);
 $/ = $orig_sep;
}

print "Started: ". &mytime."\n";

#######################################################################################
# initialize and first test
#######################################################################################
my @curoptmeta = @metastartvalues;
my ( $a, $b, $finished, @testlist, @testlisttargets, $opttarget, $optvalue );
my @snsp = evalsnsp(@curoptmeta);
my $target = sprintf( "%.4f", gettarget(@snsp) );
$opttarget = $target;
print "starting accuracy: "
  . join( ", ", @snsp )
  . ", starting target: $target\n";
my $found_improvement;
my @bindings;

#######################################################################################
# optimization loop for meta parameters
#######################################################################################

if ( !$trans_matrix ) {
 my (@testmeta);
 for ( my $r = 0 ; $r < $rounds ; $r++ ) {
  $found_improvement = 0;
  for ( my $idx = 0 ; $idx <= $#metaparnames ; $idx++ ) {
   print
"improving parameter $metaparnames[$idx] curently set to $curoptmeta[$idx]\n";
   @testmeta = @curoptmeta;

   # set the initial min and max of the range to test
   if ( $metaparranges[$idx][0] ne 'list' ) {
    $a = $metaparranges[$idx][1];
    $b = $metaparranges[$idx][2];
    print "$a-$b\n";
   }
   $finished = 0;
   while ( !$finished ) {
    $finished = 1;

    # generate a list of values to test
    @testlist = ();
    if ( $metaparranges[$idx][0] eq 'list' ) {
     @testlist = @{ $metaparranges[$idx] };
     shift @testlist;
    }
    elsif ( $metaparranges[$idx][0] eq 'floatrange' ) {
     for ( my $n = 0 ; $n < $pstep ; $n++ ) {
      push @testlist, $a + $n * ( $b - $a ) / ( $pstep - 1 );
     }
    }
    else {    # round the values
     for ( my $n = 0 ; $n < $pstep ; $n++ ) {
      my $tv = int( $a + $n * ( $b - $a ) / ( $pstep - 1 ) );
      if ($tv && $testlist[$#testlist] && $tv ne $testlist[$#testlist] ) {
       push @testlist, $tv;
      }
     }
    }
    @testlisttargets = ();
    print "$metaparnames[$idx]: " . join( "\t", @testlist ) . "\n";
    foreach my $testvalue (@testlist) {
     $testmeta[$idx] = $testvalue;    # set the parameter to the testvalue
     @snsp = evalsnsp(@testmeta);
     
     $target = sprintf( "%.4f", gettarget(@snsp) );
     push @testlisttargets, $target;
     if ( $target > $opttarget ) {    # found improvement
      $optvalue          = $testvalue;
      @curoptmeta        = @testmeta;
      $opttarget         = $target;
      $found_improvement = 1;
      print "found improvement: "
        . join( ", ", @snsp )
        . ", optimal target: $target\n";
      print "changing $metaparnames[$idx] to $optvalue\n";
      printmetavalues( \@metaparnames, \@testmeta );
      savenewpars(@curoptmeta);
      $finished = 0;
     }
    }
    print "values  " . join( "\t", @testlist ) . "\n" if @testlist;
    print "targets " . join( "\t", @testlisttargets ) . "\n" if @testlisttargets;
    if ( $finished == 0 ) {

     # determine whether further improvements are possible at all
     # and compute the new range boundaries
     if ( $metaparranges[$idx][0] eq 'list' ) {
      $finished = 1;
     }
     else {
      my ( $newa, $newb );
      $newa = $optvalue - ( $b - $a ) / ( $pstep - 1 );
      $newb = $optvalue + ( $b - $a ) / ( $pstep - 1 );
      $newa = ( $newa < $a ) ? $a : $newa;
      $newb = ( $newb > $b ) ? $b : $newb;
      $a    = $newa;
      $b    = $newb;
      if ( $metaparranges[$idx][0] eq 'intrange' ) {
       $a = int( $a + 1 );
       $b = int($b);
       if ( $b < $a ) {
        $finished = 1;
       }
      }
     }
    }
   }
  }

#	if (!$found_improvement && $r<$rounds-1) {
#	    print "Could not further improve. Skipping last ". ($rounds-$r-1) ." rounds\n";
#	    last;
#	}
 }
}
else {
#######################################################################################
 # optimization loop for transition probabilities
#######################################################################################
 my ( @trylist, @normedlist );
 if ($matrix_constraints) {
  @trylist    = getStateList("TRY");
  @normedlist = getStateList("NORMED");
 }
 else {
  @trylist    = ( 0 .. $n );    # optimize all transitions
  @normedlist = ( 0 .. $n );    # normalize all transitions
 }
 print "Try list: " .    ( join " ", @trylist ) . "\n";
 print "Normed list: " . ( join " ", @normedlist ) . "\n";
 getBindings();

 print "Optimizing transitions from these states.\n";
 my @curopttrans = @trans;
 save_trans_matrix( \@curopttrans, $trans_matrix . '.curopt' );
 my @testtrans;
 for ( my $r = 0 ; $r < $rounds ; $r++ ) {
  print "Improvement round/cycle " . ( $r + 1 ) . "\n";
  $found_improvement = 0;
  foreach my $idx (@trylist) {
   my $normed = 0;
   if ( grep /^$idx$/, @normedlist ) {
    $normed = 1;
   }
   my $normsum;
   if ($normed) {
    $normsum = 0;
    for ( my $j = 0 ; $j < $n ; $j++ ) {
     $normsum += $curopttrans[$idx][$j];
    }
   }
   my @tolist;
   my @transvec;
   for ( my $j = 0 ; $j < $n ; $j++ ) {
    if ( $curopttrans[$idx][$j] > 0 ) {
     push @tolist,   $j;
     push @transvec, $curopttrans[$idx][$j];
    }
   }

# skip state if it is normed and just one transition out of this state is possible.
   next unless ( !$normed || @tolist > 1 );

   print "Improving transitions out of state $idx\n";
   print "Nonzero transitions from state $idx: "
     . join( " ", grep ( $_ > 0, @transvec ) ) . "\n";

   # make a list with all the varied probability vectors to try
   my @tryvectors = getVariedTransVectors( \@transvec, $normsum, $normed );
   print "Trying "
     . scalar(@tryvectors) . " "
     . ( $normed ? "normed" : "unnormed" )
     . " variations of transition vector.\n";

   # change each transition probability
   # evaluate the accuracy
   foreach my $varieddist (@tryvectors) {
    print "Try varied distribution " . join( " ", @{$varieddist} ) . " ";

    # create varied transition matrix
    copyMatrix( \@testtrans, \@curopttrans, $n );
    for ( my $k = 0 ; $k < @tolist ; $k++ ) {
     $testtrans[$idx][ $tolist[$k] ] = $varieddist->[$k];
    }
    realizeBindings( $idx, \@testtrans, 0 );

    # save it to the file
    save_trans_matrix( \@testtrans, $trans_matrix );

    # start an evaluation run
    @snsp = evalsnsp();
    $target = sprintf( "%.4f", gettarget(@snsp) );
    print "\ttarget=$target\n";
    if ( $target > $opttarget ) {    # found improvement
     $opttarget         = $target;
     $found_improvement = 1;
     print "*** Found improvement: "
       . join( ", ", @snsp )
       . ", optimal target: $target\n";
     print "changing trans. probs out of state $idx from "
       . ( join " ", ( grep ( $_ > 0, @{ $curopttrans[$idx] } ) ) ) . " to "
       . ( join " ", ( grep ( $_ > 0, @{ $testtrans[$idx] } ) ) ) . "\n";
     copyMatrix( \@curopttrans, \@testtrans, $n );
     save_trans_matrix( \@curopttrans, $trans_matrix . ".curopt" );
     $finished = 0;
    }
    else {    # no improvement
              # restore file with currently optimal transition matrix
     save_trans_matrix( \@curopttrans, $trans_matrix );
    }
   }
  }
  if ( !$found_improvement && $r < $rounds - 1 ) {
   print "Could not further improve. Skipping last "
     . ( $rounds - $r - 1 )
     . " rounds\n";
   last;
  }
 }
}

#######################################################################################
# final training
#######################################################################################

if ( !$notrain ) {

 my @todelete = glob("$speciesdir/*tmp*pbl");
 foreach (@todelete) { unlink($_) }
 @todelete = glob("$output_directory/curtrain-*");
 foreach (@todelete) { unlink($_) }
 @todelete = glob("$output_directory/predictions-*");
 foreach (@todelete) { unlink($_) }

 # make the joint training file (train.gb and onlytrain.gb)
 unlink("$output_directory/curtest");
 unlink("$output_directory/curtrain");
 print "Making final training with the optimized parameters.\n";
 if ($onlytrain_gb) {
  system("cp $optimize_gb $output_directory/curtrain");
  system("cat $onlytrain_gb >> $output_directory/curtrain");
  &process_cmd(
"$etrain_exec --species=$species --AUGUSTUS_CONFIG_PATH=$config_path --/genbank/verbosity=0 $output_directory/curtrain $common_parameters $modelrestrict >/dev/null 2>/dev/null"
  );
  unlink("$output_directory/curtrain");
 }
 else {
  &process_cmd(
"$etrain_exec --species=$species --AUGUSTUS_CONFIG_PATH=$config_path --/genbank/verbosity=0 $optimize_gb $common_parameters $modelrestrict  >/dev/null 2>/dev/null"
  );
 }
}

&process_cmd("rm -rf $output_directory");

sub evalsnsp {
################################################
# evalsnsp: determine the values
# base sn, base sp, exon sn, exon sp, gene sn, gene sp, tss medianDiff, tts medianDiff
# sn: sensitivity, sp: specificity
# given a set of metaparameter values
################################################
 my @values = @_;
 my ( $cbsn, $cbsp, $cesn, $cesp, $cgsn, $cgsp, $csmd, $ctmd )
   ;    # accuracy values of current bucket
 my ( $gbsn, $gbsp, $gesn, $gesp, $ggsn, $ggsp, $gsmd, $gtmd )
   ;    # total accuracy values
 $gbsn = $gbsp = $gesn = $gesp = $ggsn = $ggsp = $gsmd = $gtmd = 0;
 my $argument = '';
 if ( !$trans_matrix ) {

  # make the parameters string for the command line
  for ( my $i = 0 ; $i <= $#metaparnames ; $i++ ) {
   $argument = $argument . " --" . $metaparnames[$i] . "=" . $values[$i];
  }

  #print "argument:$argument\n";
  # check if accuracy has already been computed for this parameter combination
  if ( exists( $storedsnsp{$argument} ) ) {
   return @{ $storedsnsp{$argument} };
  }
 }

 # Loop over the buckets and chose bucket k as the one for testing.
 # All other buckets are taken for training if appropriate
 print "bucket ";
 my $thread_helper = new Thread_helper($cpus);
 for ( my $k = 1 ; $k <= $kfold ; $k++ ) {
  my $thread = threads->create( 'start_prediction', $k, $argument );
  $thread_helper->add_thread($thread);
 }
 $thread_helper->wait_for_all_threads_to_complete();
 my @failed_threads = $thread_helper->get_failed_threads();
 if (@failed_threads) {
  die "Error, " . scalar(@failed_threads) . " threads failed.\n";
  exit(1);
 }

 # compute accuracy
 for ( my $k = 1 ; $k <= $kfold ; $k++ ) {
  open( PRED, "<$output_directory/predictions-$k.txt" );
  while (<PRED>) {
   if (/nucleotide level \| +(\S+) \| +(\S+) \|$/) {
    ( $cbsn, $cbsp ) = ( $1, $2 );
   }
   if (/exon level \|.*-- \| +(\S+) \| +(\S+) \|$/) {
    ( $cesn, $cesp ) = ( $1, $2 );
   }
   if (/gene level \|.* \| +(\S+) \| +(\S+) \|$/) {
    ( $cgsn, $cgsp ) = ( $1, $2 );
   }
   if (/TSS \|.* \|.* \|.* \| +(\S+) \|$/) {
    $csmd = $1;
   }
   if (/TTS \|.* \|.* \|.* \| +(\S+) \|$/) {
    $ctmd = $1;
   }
  }
  close(PRED);

  if (    !$cbsn
       || !$cbsp
       || !$cesn
       || !$cesp
       || !$cgsn
       || !$cgsp )
  {
   warn(
"Could not read the accuracy values out of predictions.txt when processing bucket $k."
   );
   $cbsn = int(0);
   $cbsp = int(0);
   $cesn = int(0);
   $cesp = int(0);
   $cgsn = int(0);
   $cgsp = int(0);
   $csmd = int(0);
   $ctmd = int(0);
  }

#print "accuracy on bucket$k: $cbsn, $cbsp, $cesn, $cesp, $cgsn, $cgsp, $csmd, $ctmd\n";
  $gbsn += $cbsn;
  $gbsp += $cbsp;
  $gesn += $cesn;
  $gesp += $cesp ? $cesp : int(0);
  $ggsn += $cgsn ? $cgsn : int(0);
  $ggsp += $cgsp ? $cgsp : int(0);
  $gsmd += $csmd ? $csmd : int(0);
  $gtmd += $ctmd ? $ctmd : int(0);
 }

 #print "\n";
 $gbsn = sprintf( "%.4f", $gbsn / $kfold );
 $gbsp = sprintf( "%.4f", $gbsp / $kfold );
 $gesn = sprintf( "%.4f", $gesn / $kfold );
 $gesp = sprintf( "%.4f", $gesp / $kfold );
 $ggsn = sprintf( "%.4f", $ggsn / $kfold );
 $ggsp = sprintf( "%.4f", $ggsp / $kfold );
 $gsmd = sprintf( "%.2f", $gsmd / $kfold );
 $gtmd = sprintf( "%.2f", $gtmd / $kfold );

 my @returnarray = ( $gbsn, $gbsp, $gesn, $gesp, $ggsn, $ggsp, $gsmd, $gtmd );
 $storedsnsp{$argument} = \@returnarray;
 print "\n";
 return @returnarray;
}

sub start_prediction() {
 my $k        = shift;
 my $argument = shift;
 my $pbloutfiles =
"--/ExonModel/outfile=exon-tmp$k.pbl --/IntronModel/outfile=intron-tmp$k.pbl --/IGenicModel/outfile=igenic-tmp$k.pbl --/UtrModel/outfile=utr-tmp$k.pbl";
 my $pblinfiles =
   $onlyutr
   ? ''
   : "--/ExonModel/infile=exon-tmp$k.pbl --/IntronModel/infile=intron-tmp$k.pbl --/IGenicModel/infile=igenic-tmp$k.pbl --/UtrModel/infile=utr-tmp$k.pbl";

 if ( !$notrain ) {
  &process_cmd(
"$etrain_exec --species=$species --AUGUSTUS_CONFIG_PATH=$config_path $argument $common_parameters $be_silent $modelrestrict $pbloutfiles $output_directory/curtrain-$k  >/dev/null 2>/dev/null"
  );
 }
 else {
  $pblinfiles = '';
  ; # training did not take place, so the $pbloutfiles have not beeen created and cannot be used for prediction
 }

 &process_cmd(
"$augustus_exec --genemodel=complete --species=$species --AUGUSTUS_CONFIG_PATH=$config_path $argument $common_parameters $pblinfiles $output_directory/bucket$k.gb > $output_directory/predictions-$k.txt 2>/dev/null"
 );
 print "$k ";

}

sub gettarget {
######################################################################################
# gettarget: get an optimization target value from
# base sn, base sp, exon sn, exon sp, gene sn, gene sp, tss medianDiff, tts medianDiff
# feel free to change the weights
######################################################################################
 my ( $bsn, $bsp, $esn, $esp, $gsn, $gsp, $smd, $tmd ) = @_;
 return ( 3 * $bsn +
          3 * $bsp +
          4 * $esn +
          4 * $esp +
          2 * $gsn +
          2 * $gsp +
          40 / ( $smd + 40 ) +
          40 / ( $tmd + 40 ) ) / 20;

 #   return (3*$bsn + 9*$bsp + 4*$esn + 12*$esp + 2*$gsn + 6*$gsp)/36;
}

sub savenewpars {
################################################
 # savenewpars: replace parameters in the file $species_cfg_filename
################################################
 open( SPCCFG, ">$species_cfg_filename" )
   or die("Could not open $species_cfg_filename");
 print "Writing new parameters to $species_cfg_filename...\n";
 my $parname;
 foreach my $line (@spcfilelines) {
  if ( $line =~ /^\s*\#.*/ || $line =~ /^\s*$/ ) {
   print SPCCFG $line;    # print unchanged line
   next;
  }

  # format:
  # parname   value   # comment
  # or
  # parname   value
  if ( $line =~ /^(\s*)(\S+)(\s+)(\S+)(.*)$/ ) {
   $parname = $2;
   my $index = idx( \@metaparnames, $parname );
   if ( $index != -1 ) {
    print SPCCFG $1 . $parname . $3 . $_[$index] . $5 . "\n";
   }
   else {
    print SPCCFG $line;
   }
  }
  else {
   print SPCCFG $line;
  }
 }
 close(SPCCFG);
}

sub save_trans_matrix {
################################################
 # save_trans_matrix: replace existing transition
 # probabilities with the given new ones
################################################
 my $newtransref = shift;
 my $filename    = shift;
 open( TRANS, ">$filename" ) or die("Could not open $filename for writing.");
 foreach my $line (@transfilelines) {
  if ( $line =~ /^(\s*)(\d+)(\s+)(\d+)(\s+)(\S+)/ ) {
   my $value = 0;
   $value = $newtransref->[$2][$4] unless ( !defined $newtransref->[$2][$4] );
   print TRANS $1 . $2 . $3 . $4 . $5 . $value . "\n";
  }
  else {
   print TRANS $line;
  }
 }
 close TRANS;
}

sub printmetavalues {
 my ( $names, $values ) = @_;
 for ( my $i = 0 ; $i <= $#$names ; $i++ ) {
  print $names->[$i] . "\t" . $values->[$i] . "\n";
 }
}

sub printmetaranges {
 my ( $names, $ranges ) = @_;
 print "metapar ranges $#$names:\n";
 for ( my $i = 0 ; $i <= $#$names ; $i++ ) {
  print $names->[$i] . "\t" . join( " ", @{ $ranges->[$i] } ) . "\n";
 }
}

#
# idx: find the index to an element in an array of strings
#

sub idx {
 my $arrayref = shift;
 my $element  = shift;

 for ( my $i = 0 ; $i <= $#$arrayref ; $i++ ) {
  if ( $arrayref->[$i] eq $element ) {
   return $i;
  }
 }
 return -1;
}

sub norm {

 #
 # norm a vector to sum up to $normsum
 # if $normed is true
 my ( $vecref, $normsum, $normed ) = @_;
 if ($normed) {
  my $sum = 0;
  foreach my $item ( @{$vecref} ) {
   $sum += $item;
  }
  if ( $sum > 0 && $normsum > 0 ) {
   foreach my $item ( @{$vecref} ) {
    $item *= $normsum / $sum;
   }
  }
 }

 # round to 6 places
 foreach my $item ( @{$vecref} ) {
  $item = sprintf( "%.6f", $item );
  $item =~ s/(\.\d*[1-9])0+$/$1/;    # remove trailing zeros
 }
}

sub getStateList {

 #
 # read try list or normed list from matrix_constraints file
 #
 my $identifyer = shift;             # 'TRY' or 'NORMED'
 if ( !-s $matrix_constraints ) {
  print
"Could not find the file with the transition matrix constraints: $matrix_constraints\n";
  return;
 }
 my @statelist;
 open CONSTR, $matrix_constraints or die;
 my $scanning = 0;
 my $all      = 0;
 while (<CONSTR>) {
  if (/\[$identifyer\]/) {
   $scanning = 1;
  }
  elsif (/\s*\[.*\]/) {
   $scanning = 0;
  }
  if ( $scanning && /^\s*(\d+)\s*/ ) {
   push @statelist, $1;
  }
  if ( $scanning && /^\s*all/ ) {
   $all = 1;
  }
 }
 close CONSTR;
 if ($all) {
  @statelist = ( 0 .. $n );
 }
 return @statelist;
}

sub getBindings {

 #
 # read bindings list from matrix_constraints file
 # into global variables
 #
 #
 if ( !-s $matrix_constraints ) {
  print
"Could not find the file with the transition matrix constraints: $matrix_constraints\n";
  return;
 }
 my @bindingsstr;
 open( BINDINGS, $matrix_constraints ) or die;
 my $scanning = 0;
 while (<BINDINGS>) {
  if (/\[BINDINGS\]/) {
   $scanning = 1;
  }
  elsif (/\s*\[.*\]/) {
   $scanning = 0;
  }
  if ($scanning) {
   $_ =~ s/#.*//;
   if ( /^\s*(\(.*\))/ || /^(MC\S*)/ ) {
    push @bindingsstr, $1;
   }
  }
 }
 close BINDINGS;
 print "bindings: " . ( join "\n", @bindingsstr ) . "\n";
 my $btype;

# @bindings is a list of references of bindings.
# Each binding is a list of
# - the binding equation string
# - a type, (either '+' or '/')
# - leftTrans
# - rightTrans
# - allStates, the list of states the transitions originate in
# leftTrans and rightTrans are the transition lists on the left-hand side and right-hand side of the equation
# They each are lists of transitions.
# lt1 + lt2 + ... + ltA = rt1 + rt2 + ... + rtB
# A and B can be 1.
# example (0,24)+(0,25)=(0,65)+(0,70)
#
 foreach my $bstr (@bindingsstr) {
  if ( $bstr =~ '/' ) {
   $btype = '/';
  }
  elsif ( $bstr =~ /^MC/ ) {
   $btype = 'M';
  }
  else {
   $btype = '+';
  }
  my @eqsides = split /=/, $bstr;
  if ( @eqsides != 2 ) {
   print
"Error: Wrong format in $matrix_constraints. Each binding must be an equation.\n";
   return;
  }
  my @leftTrans  = parseTransList( $eqsides[0] );  # left hand side of equation
  my @rightTrans = parseTransList( $eqsides[1] );  # right hand side of equation
  my @allStates  = ();
  my %seen       = ();
  foreach my $trans ( ( @leftTrans, @rightTrans ) ) {
   push @allStates, $trans->[0]
     unless $seen{ $trans->[0] }++;                # push the originating state
  }
  push @bindings, [ $bstr, $btype, \@leftTrans, \@rightTrans, \@allStates ];
 }
}

sub realizeBindings {

 #
 # realizeBindings
 #
 # returns error message string if there are any
 my $state = shift;    # line that was just changed
 my $trans = shift;    # transition matrix that may need to be adjusted
 my $doNothingButComplain = shift;    # only error msgs
 my $errmsg;
 print "realizeBindings $state\n";
 foreach my $binding (@bindings) {
  my ( $bstr, $btype, $leftTrans, $rightTrans, $allStates ) = @{$binding};
  if ( grep /^$state$/, @{$allStates} ) {
   print "binding applies: $bstr\n";
   if ( $btype eq '+' ) {
    if ( @{$allStates} == 1 ) {

     # binding applies just to $state
     # lt1 + lt2 + ... + ltA   rt1 + rt2 + ... + rtB
     # ---------------------   ---------------------
     #       lhs                      rhs
     my $lhs = 0;
     my $rhs = 0;
     foreach my $t ( @{$leftTrans} ) {
      $lhs += $trans->[$state][ $t->[1] ];
     }
     foreach my $t ( @{$rightTrans} ) {
      $rhs += $trans->[$state][ $t->[1] ];
     }

     #print "actual value lhs=$lhs, rhs=$rhs\n";
     if ( $lhs != $rhs ) {
      if ($doNothingButComplain) {
       $errmsg .= "Binding $bstr not satisfied: $lhs != $rhs\n";
      }
      else {

       # rescale the probabilities
       if ( $lhs > 0 && $rhs > 0 ) {
        foreach my $t ( @{$leftTrans} ) {
         $trans->[$state][ $t->[1] ] *= ( $lhs + $rhs ) / 2 / $lhs;
        }
        foreach my $t ( @{$rightTrans} ) {
         $trans->[$state][ $t->[1] ] *= ( $lhs + $rhs ) / 2 / $rhs;
        }
        roundVector( @{ $trans->[$state] } );
        print "After applying binding: "
          . ( join " ", grep { $_ != 0 } @{ $trans->[$state] } ) . "\n";
       }
      }
     }
    }
    else {

# binding applies to other states as well
# change only the other states and not the transition probabilities from this state
#
     my $restsum
       ; #sum of of all transition probabilities except the ones from $state. Those are fixed.
     foreach my $t ( @{$rightTrans} ) {
      if ( $t->[0] == $state ) {
       $restsum += $trans->[ $t->[0] ][ $t->[1] ];
      }
     }
     foreach my $t ( @{$leftTrans} ) {
      if ( $t->[0] == $state ) {
       $restsum -= $trans->[ $t->[0] ][ $t->[1] ];
      }
     }

     # compute the actual restsum
     my $actualrestsum;
     foreach my $t ( @{$rightTrans} ) {
      if ( $t->[0] != $state ) {
       $actualrestsum -= $trans->[ $t->[0] ][ $t->[1] ];
      }
     }
     foreach my $t ( @{$leftTrans} ) {
      if ( $t->[0] != $state ) {
       $actualrestsum += $trans->[ $t->[0] ][ $t->[1] ];
      }
     }
     if ( $restsum * $actualrestsum > 0 ) {

      # rescale the other states' transition probabilities
      my %added = ();
      foreach my $t ( @{$rightTrans}, @{$leftTrans} ) {
       if ( $t->[0] != $state ) {
        $added{ $t->[0] } +=
          $trans->[ $t->[0] ][ $t->[1] ] * ( $restsum / $actualrestsum - 1 );
        $trans->[ $t->[0] ][ $t->[1] ] *= $restsum / $actualrestsum;
       }
      }

      # now adjust the other states transition vectors
      foreach my $ostate ( @{$allStates} ) {
       if ( $ostate != $state ) {
        print "trans probs of other state $ostate: "
          . ( join " ", grep { $_ != 0 } @{ $trans->[$ostate] } ) . "\n";

        # rescale all other transition probabilities, such
        # that their sum is decreased by $added{$ostate}
        # print "sum must be decreased by $added{$ostate}\n";
        my %otochanged = ();
        foreach my $t ( @{$rightTrans}, @{$leftTrans} ) {
         if ( $t->[0] == $ostate ) {
          $otochanged{ $t->[1] } = 1;
         }
        }
        my $osum;
        for ( my $j = 0 ; $j < $n ; $j++ ) {
         if ( !$otochanged{$j} ) {
          $osum += $trans->[$ostate][$j];
         }
        }

      # now renorm so that the other trans probs sum up to $osum-$added{$ostate}
        if ( $osum > 0 && $osum - $added{$ostate} > 0 ) {
         for ( my $j = 0 ; $j < $n ; $j++ ) {
          if ( !$otochanged{$j} ) {
           if ( $trans->[$ostate][$j] > 0 ) {
            print "($ostate,$j) is changed from $trans->[$ostate][$j] to ";
            $trans->[$ostate][$j] *= ( $osum - $added{$ostate} ) / $osum;
            print "$trans->[$ostate][$j]\n";
           }
          }
         }

        }
        roundVector( @{ $trans->[$ostate] } );
        print "fixed trans probs of other state $ostate: "
          . ( join " ", grep { $_ != 0 } @{ $trans->[$ostate] } ) . "\n";
       }
      }
     }
     else {

      # Don't fix distribution matrix. Should not happen.
     }
    }

   }
   elsif ( $btype eq '/' ) {

    # not yet implemented
   }
   elsif ( $btype eq 'M' ) {

# markov chain binding
# make sure that one 3x3 Markov chain transition matrix is from the time-reversed
# process of another Markov chain. Assume stationarity.

    my $leftm = [
                  [
                    $trans->[ $leftTrans->[0]->[0] ][ $leftTrans->[0]->[1] ],
                    $trans->[ $leftTrans->[1]->[0] ][ $leftTrans->[1]->[1] ],
                    $trans->[ $leftTrans->[2]->[0] ][ $leftTrans->[2]->[1] ]
                  ],
                  [
                    $trans->[ $leftTrans->[3]->[0] ][ $leftTrans->[3]->[1] ],
                    $trans->[ $leftTrans->[4]->[0] ][ $leftTrans->[4]->[1] ],
                    $trans->[ $leftTrans->[5]->[0] ][ $leftTrans->[5]->[1] ]
                  ],
                  [
                    $trans->[ $leftTrans->[6]->[0] ][ $leftTrans->[6]->[1] ],
                    $trans->[ $leftTrans->[7]->[0] ][ $leftTrans->[7]->[1] ],
                    $trans->[ $leftTrans->[8]->[0] ][ $leftTrans->[8]->[1] ]
                  ]
    ];
    my $rightm = [
                   [
                     $trans->[ $rightTrans->[0]->[0] ][ $rightTrans->[0]->[1] ],
                     $trans->[ $rightTrans->[1]->[0] ][ $rightTrans->[1]->[1] ],
                     $trans->[ $rightTrans->[2]->[0] ][ $rightTrans->[2]->[1] ]
                   ],
                   [
                     $trans->[ $rightTrans->[3]->[0] ][ $rightTrans->[3]->[1] ],
                     $trans->[ $rightTrans->[4]->[0] ][ $rightTrans->[4]->[1] ],
                     $trans->[ $rightTrans->[5]->[0] ][ $rightTrans->[5]->[1] ]
                   ],
                   [
                     $trans->[ $rightTrans->[6]->[0] ][ $rightTrans->[6]->[1] ],
                     $trans->[ $rightTrans->[7]->[0] ][ $rightTrans->[7]->[1] ],
                     $trans->[ $rightTrans->[8]->[0] ][ $rightTrans->[8]->[1] ]
                   ]
    ];

    my $revLeft = 0;    # 1 iff the left Markov chain needs to be adjusted
    for ( my $i = 0 ; $i < 9 ; $i++ ) {
     if ( $rightTrans->[$i]->[0] == $state ) {
      $revLeft = 1;
     }
    }
    my ( $m, $rm );
    if ($revLeft) {
     $m  = $rightm;
     $rm = $leftm;
    }
    else {
     $m  = $leftm;
     $rm = $rightm;
    }

    reverseMatrix( $m, $rm );
    if ($revLeft) {
     $trans->[ $leftTrans->[0]->[0] ][ $leftTrans->[0]->[1] ] = $rm->[0][0];
     $trans->[ $leftTrans->[1]->[0] ][ $leftTrans->[1]->[1] ] = $rm->[0][1];
     $trans->[ $leftTrans->[2]->[0] ][ $leftTrans->[2]->[1] ] = $rm->[0][2];
     $trans->[ $leftTrans->[3]->[0] ][ $leftTrans->[3]->[1] ] = $rm->[1][0];
     $trans->[ $leftTrans->[4]->[0] ][ $leftTrans->[4]->[1] ] = $rm->[1][1];
     $trans->[ $leftTrans->[5]->[0] ][ $leftTrans->[5]->[1] ] = $rm->[1][2];
     $trans->[ $leftTrans->[6]->[0] ][ $leftTrans->[6]->[1] ] = $rm->[2][0];
     $trans->[ $leftTrans->[7]->[0] ][ $leftTrans->[7]->[1] ] = $rm->[2][1];
     $trans->[ $leftTrans->[8]->[0] ][ $leftTrans->[8]->[1] ] = $rm->[2][2];
    }
    else {
     $trans->[ $rightTrans->[0]->[0] ][ $rightTrans->[0]->[1] ] = $rm->[0][0];
     $trans->[ $rightTrans->[1]->[0] ][ $rightTrans->[1]->[1] ] = $rm->[0][1];
     $trans->[ $rightTrans->[2]->[0] ][ $rightTrans->[2]->[1] ] = $rm->[0][2];
     $trans->[ $rightTrans->[3]->[0] ][ $rightTrans->[3]->[1] ] = $rm->[1][0];
     $trans->[ $rightTrans->[4]->[0] ][ $rightTrans->[4]->[1] ] = $rm->[1][1];
     $trans->[ $rightTrans->[5]->[0] ][ $rightTrans->[5]->[1] ] = $rm->[1][2];
     $trans->[ $rightTrans->[6]->[0] ][ $rightTrans->[6]->[1] ] = $rm->[2][0];
     $trans->[ $rightTrans->[7]->[0] ][ $rightTrans->[7]->[1] ] = $rm->[2][1];
     $trans->[ $rightTrans->[8]->[0] ][ $rightTrans->[8]->[1] ] = $rm->[2][2];
    }
   }
  }
 }
}

sub parseTransList {

 #
 # parseTransList
 # parses a string like (0,24)+(0,25) or (26,28)/(26,29)
 # into a list of pairs
 my $tstr      = shift;
 my @transList = ();
 foreach my $tok ( split /\)[\(\)\+\*,]*\(/, $tstr ) {
  $tok =~ /(\d+),\s*(\d+)/;
  push @transList, [ $1, $2 ];
 }
 return @transList;
}

sub getVariedTransVectors {

 #
 # getVariedTransVectors
 # parameters: \@transvec, $normsum, $normed
 # returns a list of references to varied transition probability vectors
 #
 my $transvec   = shift;
 my $normsum    = shift;
 my $normed     = shift;
 my @tryvectors = ();      #list of references to transition probability vectors
 my $end;                  # vary components up to this index
 my $maxfactor  = 2;       # factor is between 1 and this
 my $variations = 3;

 if ( !$normed || @{$transvec} > 2 ) {
  $end = @{$transvec};
 }
 else {
  $end = 1;
 }

 # simply enlarge and decrease each single transition for now
 for ( my $k = 0 ; $k < $end ; $k++ ) {
  for ( my $r = 0 ; $r < $variations ; $r++ ) {

   # change element k
   my @vartransvec = @{$transvec};
   my $factor      = 1 + rand( $maxfactor - 1.0 );
   if ( rand(2) < 1 ) {
    $factor = 1 / $factor;
   }
   $vartransvec[$k] = $transvec->[$k] * $factor;
   norm( \@vartransvec, $normsum, $normed );
   push @tryvectors, \@vartransvec;
  }
 }
 return @tryvectors;
}

sub copyMatrix {

 #
 # copyMatrix
 # copy a tranistion matrix
 #
 my $new  = shift;
 my $old  = shift;
 my $size = shift;
 for ( my $i = 0 ; $i < $size ; $i++ ) {
  $new->[$i] = [];
  for ( my $j = 0 ; $j < $size ; $j++ ) {
   push @{ $new->[$i] }, $old->[$i][$j];
  }
 }
}

sub roundVector {

 #
 # roundVector
 #
 my $sum = 0;
 foreach my $v (@_) {
  $sum += $v;
  $v = sprintf( "%.6f", $v );
 }
 my $roundedsum = sprintf( "%.6f", $sum );    # is usually 1
 my $newsum     = 0;
 my $maxel      = 0;
 for ( my $i = 0 ; $i < @_ ; $i++ ) {
  $newsum += $_[$i];
  if ( $_[$i] > $_[$maxel] ) {
   $maxel = $i;
  }
 }

# in case the rounding changed the sum a little, adjust the largest element only
 if ( $newsum != $roundedsum ) {
  $_[$maxel] += $roundedsum - $newsum;
 }
}

sub matrixSquare {

 #
 # matrixSquare
 #
 my $m = shift;    # reference to quadratic matrix

 my $n = @{$m};    # dimension nxn
 my @ret;
 for ( my $i = 0 ; $i < $n ; $i++ ) {
  $ret[$i] = [];
  for ( my $j = 0 ; $j < $n ; $j++ ) {
   for ( my $k = 0 ; $k < $n ; $k++ ) {
    $ret[$i][$j] += $m->[$i][$k] * $m->[$k][$j];
   }
  }
 }
 return \@ret;
}

sub printMatrix {

 #
 # printMatrix
 #
 my $m = shift;
 my $n = @{$m};    # dimension nxn
 for ( my $i = 0 ; $i < $n ; $i++ ) {
  for ( my $j = 0 ; $j < $n ; $j++ ) {
   print $m->[$i][$j] . "\t";
  }
  print "\n";
 }
}

sub reverseMatrix {

 #
 # reverseMatrix
 #
 #
 my $m  = shift;
 my $rm = shift;

 # first scale both matrices so they are real transition matrices
 # and remember the scaling factor
 #
 my ( @scale, @scaler );
 for ( my $i = 0 ; $i < @{$m} ; $i++ ) {
  for ( my $j = 0 ; $j < @{ $m->[$i] } ; $j++ ) {
   $scale[$i] += $m->[$i][$j];
  }
  for ( my $j = 0 ; $j < @{ $m->[$i] } ; $j++ ) {
   $m->[$i][$j] /= $scale[$i];
  }
 }
 for ( my $i = 0 ; $i < @{$rm} ; $i++ ) {
  for ( my $j = 0 ; $j < @{ $rm->[$i] } ; $j++ ) {
   $scaler[$i] += $rm->[$i][$j];
  }
  for ( my $j = 0 ; $j < @{ $rm->[$i] } ; $j++ ) {
   $rm->[$i][$j] /= $scaler[$i];
  }
 }
 my $pm;
 copyMatrix( \@{$pm}, $m, scalar( @{$m} ) );
 for ( my $i = 0 ; $i < 5 ; $i++ ) {
  $pm = matrixSquare($pm);
 }
 my @pi = @{ $pm->[0] };    # stationary distribution
      #print "stationary distribution: " . join(" ", @pi) . "\n";
      # rtime-everse other matrix rm
      # according to the formula
      # Q_ij = pi_j * P_ji / pi_i
      # Q: transition matrix of time-reversed chain
      # P: transition matrix of normal chain
      # pi: stationary distribution

 for ( my $i = 0 ; $i < @{$rm} ; $i++ ) {

  for ( my $j = 0 ; $j < @{ $rm->[$i] } ; $j++ ) {
   $rm->[$i][$j] = $pi[$j] * $m->[$j][$i] / $pi[$i];
  }
  roundVector( @{ $rm->[$i] } );
 }
 print "reversed chain:\n";
 printMatrix($rm);

 # rescale both chains so the sums are like before
 for ( my $i = 0 ; $i < @{$rm} ; $i++ ) {
  for ( my $j = 0 ; $j < @{ $rm->[$i] } ; $j++ ) {
   $rm->[$i][$j] *= $scaler[$i];
  }
  for ( my $j = 0 ; $j < @{ $m->[$i] } ; $j++ ) {
   $m->[$i][$j] *= $scale[$i];
  }
 }
}

sub got_interrupt_signal {

 print STDERR "$0 was interrupted.\n";
 if ($trans_matrix) {
  if ( -s $trans_matrix . ".curopt" ) {
   system("cp $trans_matrix.curopt $trans_matrix");
   print STDERR "I replaced the transition matrix in $trans_matrix "
     . "with the currently optimal matrix.\n";
  }
 }
 else {
  print STDERR
    "\nPlease retrain augustus with the new parameters using etraining.\n";
 }
 exit(1);
}

sub check_options() {
 $optimize_gb = shift if !$optimize_gb;
 pod2usage("Optimization training file missing\n")
   unless ( $optimize_gb && -s $optimize_gb );
 pod2usage("No species specified\n") if ( !$species );
 pod2usage("--kfold must be at least two fold cross validation")
   if ( $kfold < 1 );
 $cpus = 1 if !$cpus || $cpus < 1;
 $ENV{'PATH'} .= ':' . $exec_dir if $exec_dir && -d $exec_dir;
 die(
"Augustus config directory not in environment (\$AUGUSTUS_CONFIG_PATH) and not specified.\n"
 ) unless $config_path && -d $config_path;
 $config_path .= '/' unless $config_path =~ /\/$/;
 $notrain = 1 if ($trans_matrix);

 if ($onlyutr) {
  $modelrestrict =
"--/EHMMTraining/statecount=2 --/EHMMTraining/state00=intronmodel --/EHMMTraining/state01=utrmodel --/IntronModel/outfile=/dev/null";
  $utr = 1;
 }

 $common_parameters .= "--UTR=on"                          if $utr;
 $common_parameters .= " --translation_table=$trans_table" if $trans_table;
 $common_parameters .= " --genemodel=$genemodel"           if $genemodel;
 $common_parameters .= " --/Constant/min_coding_len=" . $min_coding
   if $min_coding;

 $output_directory = $utr ? "tmp_opt_utr_$species" : "tmp_opt_$species";
 &process_cmd("rm -rf $output_directory;mkdir $output_directory");
 die unless -d $output_directory;

 &split_genbank( $optimize_gb, $kfold );

 &create_cross_validation();
}

sub create_cross_validation() {

 # create training sets for cross-validation (parallel version)
 if ( $cpus > 1 ) {
  for ( my $k = 1 ; $k <= $kfold ; $k++ ) {

   # make the temporary training and testing files
   system("rm -f $output_directory/curtrain-$k");
   for ( my $m = 1 ; $m <= $kfold ; $m++ ) {
    if ( $m != $k ) {
     system(
          "cat $output_directory/bucket$m.gb >> $output_directory/curtrain-$k");
    }
   }
   if ($onlytrain_gb) {
    system("cat $onlytrain_gb >> $output_directory/curtrain-$k");
   }
  }
 }
}

sub split_genbank() {
 my $gb_file    = shift;
 my $partitions = shift;

 # Split $gb_file into $partitions equal portions
 print "Splitting training file into $partitions buckets...\n";
 open( TRAINGB, $gb_file ) or die("Could not open $gb_file");

 my @seqlist = ();
 push @seqlist, "12";
 @seqlist = <TRAINGB>;
 my @namelines = grep /^LOCUS   +/, @seqlist;
 my @names = ();

 if ( @namelines < $partitions ) {
  print "Number of training sequences is too small\n";
  exit;
 }

 foreach (@namelines) {
  /LOCUS +([^ ]+) */;
  push @names, $1;
 }

 my $bucket    = 0;
 my %bucketmap = ();
 srand(88);
 while ( $#names >= 0 ) {
  my $rand = rand(@names);
  $bucketmap{ $names[$rand] } = $bucket;
  $bucket++;
  if ( $bucket == $partitions ) {
   $bucket = 0;
  }
  splice @names, $rand, 1;    # delete array element
 }
 my $handle;
 my @fh = ();
 for ( $bucket = 1 ; $bucket <= $partitions ; $bucket++ ) {
  $handle = IO::File->new(">$output_directory/bucket$bucket.gb")
    or die("Could not open bucket$bucket.gb");
  push @fh, $handle;
 }
 my $orig_sep = $/;
 $/ = "\n//\n"
   ; # this causes a huge single chunk when DOS carriage returns are used at line breaks
 seek( TRAINGB, 0, 0 );
 my $nloci = 0;
 while (<TRAINGB>) {
  my $gendaten = $_;
  m/^LOCUS +(\S+) .*/;
  my $genname = $1;

  $bucket = $bucketmap{$genname};
  my $handle = $fh[$bucket];
  print $handle $gendaten;
  $nloci++;
 }
 $/ = $orig_sep;
 foreach my $handle (@fh) {
  close $handle;
 }

 if ( $nloci < @namelines ) {
  die( "Genbank input file appears to have fewer records than expected.\n"
   . "This could be a consequence of using DOS (Windows) carriage return symbols at line breaks."
  );
 }
}

sub process_cmd {
 my ( $cmd, $dir ) = @_;
 print &mytime . "CMD: $cmd\n" if $verbose;
 chdir($dir) if $dir;
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  die "Error, cmd died with ret $ret\n";
 }
 return;
}

sub mytime() {
 my @mabbr =
   qw(January February March April May June July August September October November December);
 my @wabbr = qw(Sunday Monday Tuesday Wednesday Thursday Friday Saturday);
 my $sec   = localtime->sec() < 10 ? '0' . localtime->sec() : localtime->sec();
 my $min   = localtime->min() < 10 ? '0' . localtime->min() : localtime->min();
 my $hour =
   localtime->hour() < 10 ? '0' . localtime->hour() : localtime->hour();
 my $wday = $wabbr[ localtime->wday ];
 my $mday = localtime->mday;
 my $mon  = $mabbr[ localtime->mon ];
 my $year = localtime->year() + 1900;
 return "$wday, $mon $mday, $year: $hour:$min:$sec\t";
}

sub check_program() {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  pod2usage "Error, path to a required program ($prog) cannot be found\n\n"
    unless $path =~ /^\//;
  chomp($path);
  $path = readlink($path) if -l $path;
  push( @paths, $path );
 }
 return @paths;
}
