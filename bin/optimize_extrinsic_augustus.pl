#!/usr/bin/env perl

=pod

=head1 TODO

print out new hint file with only the src that we are testing

=head1 NAME

Optimize extrinsic configuration for Augustus. Based on optimize_augustus by Mario Stanke, 23.04.2007 

=head1 USAGE 

Mandatory parameters:

 --species                prefix of the species name
 --optimize               genbank file for training with bona fide gene structures
 --hints                  hints file to use
 --extrinsic              File with the names and their ranges of the meta parameters that are subject to optimization (default: generic_extrinsic.cfg)

Optional parameters:

 --config_path            Specify the config directory d if not set as environment variable
 --aug_exec_dir           Path to augustus and etraining executable. If not specified it must be in $PATH environment variable
 --cpus                   The number of CPUs to use (default: 1)
 --training_set           an optional genbank file that is used in addition to train.gb but only for etrain not for intermediate evaluation of accuracy. These genes may e.g. be incomplete.
 --utr                    Switch on UTR. Recommended to avoid noncoding exons being wrongly predicted as coding.
 --pstep                  For integer and floating parameters start with p tests equidistributed in the allowed range of values (default: 5)
 --min_coding             Minimum coding length
 --genemodel              Augustus genemodel parameter
 --notrain                The training step (etraining) is omitted. Use this if the HMM is already trained. This program does not optimize the HMM


=cut

use strict;
use warnings;
use Data::Dumper;
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
     $species,     $optimize_gb,      $extrinsic,
     $rounds,      $training_set,     $pstep,
     $config_path, $cpus,             $utr,
     $exec_dir,    $output_directory, $species_config_dir,
     $notrain,     $trans_table,      $genemodel,
     $min_coding,  $hint_file,        $verbose
);
$|           = 1;
$rounds      = 1;
$pstep       = 5;
$cpus        = 1;
$exec_dir    = $ENV{'AUGUSTUS_PATH'} . '/bin' if $ENV{'AUGUSTUS_PATH'};
$config_path = $ENV{'AUGUSTUS_PATH'} . '/config' if $ENV{'AUGUSTUS_PATH'};
$config_path = $ENV{'AUGUSTUS_CONFIG_PATH'} if $ENV{'AUGUSTUS_CONFIG_PATH'};

my ($common_parameters) = '';

&GetOptions(
             'species:s'              => \$species,
             'optimize:s'             => \$optimize_gb,
             'extrinsic:s'            => \$extrinsic,
             'rounds:i'               => \$rounds,
             'training_set:s'         => \$training_set,
             'pstep:i'                => \$pstep,
             'AUGUSTUS_CONFIG_PATH:s' => \$config_path,
             'cpus|threads:i'         => \$cpus,
             'UTR'                    => \$utr,
             'aug_exec_dir:s'         => \$exec_dir,
             'notrain'                => \$notrain,
             'translation_table:i'    => \$trans_table,
             'genemodel:s'            => \$genemodel,
             'min_coding_len:i'       => \$min_coding,
             'hints:s'                => \$hint_file,
             'output:s'               => \$output_directory,
             'verbose'                => \$verbose
);

# globals
&check_options();
my ( $augustus_exec, $etrain_exec ) = &check_program( 'augustus', 'etraining' );
my $cwd = `pwd`;
chomp($cwd);
my $be_silent =
" --/augustus/verbosity=0 --/ExonModel/verbosity=0 --/IGenicModel/verbosity=0 --/IntronModel/verbosity=0 --/UtrModel/verbosity=0 --/genbank/verbosity=0 ";
my @feature_order =
  qw/start stop tss tts ass dss exonpart exon intronpart intron CDSpart CDS UTRpart UTR irpart nonexonpart genicpart/;

print "Started: " . &mytime . "\n";

# we only need to train once
system("cat $optimize_gb $training_set > train.set 2>/dev/null");
&process_cmd(
"$etrain_exec --species=$species --AUGUSTUS_CONFIG_PATH=$config_path $common_parameters $be_silent train.set  >/dev/null 2>/dev/null "
) unless $notrain;
unlink("train.set");

my %sources_that_need_to_be_checked;
my @feature_headers;

my ( $metaextrinsic_hash_ref, $cfg_extra, $basic_cfg ) =
  &parse_extrinsic_meta();

#######################################################################################
# initialize and first test
#######################################################################################
my $extrinsic_iteration_counter : shared;
$extrinsic_iteration_counter = int(0);

my $thread_helper = new Thread_helper($cpus);
my $thread        = threads->create('run_evaluation');
$thread_helper->add_thread($thread);
my $todo = int(1);

foreach my $src ( keys %sources_that_need_to_be_checked ) {
 next if $src eq 'M';
 foreach my $feat (@feature_headers) {
  next if $feat eq 'feature';
  my $bonus_range = $metaextrinsic_hash_ref->{'bonus'}->{$feat}->{'value'};
  my $malus_range = $metaextrinsic_hash_ref->{'malus'}->{$feat}->{'value'};
  foreach my $bonus ( @{$bonus_range} ) {
   foreach my $malus ( @{$malus_range} ) {
    my $parameter_range = $metaextrinsic_hash_ref->{$src}->{$feat}->{'value'};
    if ( scalar(@$parameter_range) > 1 ) {
     foreach my $value (@$parameter_range) {
      $todo++;
     }
    }
   }
  }
 }
}

print "Will search for $todo parameters with $cpus CPUs\n";
sleep(1);

foreach my $src ( keys %sources_that_need_to_be_checked ) {
 next if $src eq 'M';

 # we really expect one line of evidence to be checked at a time...!

 foreach my $feat (@feature_headers) {
  next if $feat eq 'feature';
  my $bonus_range = $metaextrinsic_hash_ref->{'bonus'}->{$feat}->{'value'};
  my $malus_range = $metaextrinsic_hash_ref->{'malus'}->{$feat}->{'value'};

  # for every bonus range
  foreach my $bonus ( @{$bonus_range} ) {

   # for every malus range
   foreach my $malus ( @{$malus_range} ) {
    my $parameter_range = $metaextrinsic_hash_ref->{$src}->{$feat}->{'value'};
    if ( scalar(@$parameter_range) > 1 ) {
     foreach my $value (@$parameter_range) {

      # checks every value against the basic_cfg
      my $cfg = $basic_cfg;
      $cfg->{$src}->{$feat}->{'value'}    = $value;
      $cfg->{'bonus'}->{$feat}->{'value'} = $bonus;
      $cfg->{'malus'}->{$feat}->{'value'} = $malus;
      sleep(3);
      my $thread = threads->create( 'write_extrinsic_cfg', $cfg, $cfg_extra );
      $thread_helper->add_thread($thread);
     }
    }
   }
  }
 }
}
$thread_helper->wait_for_all_threads_to_complete();
my @failed_threads = $thread_helper->get_failed_threads();
if (@failed_threads) {
 die "Error, " . scalar(@failed_threads) . " threads failed.\n";
 exit(1);
}

### FINISHED

my @results = `grep Accuracy $output_directory/*log | sort -rnk 11`;
if ( !$results[0] ) {
 die
"No output files produced. This is weird, maybe you requested it to stop or something else happened? Report please!\n";
}
if ( scalar(@results) == $todo ) {
 print "Done and all good!\n";
}
else {
 print
"Done but something happened, not all results may have been completed ($todo != "
   . scalar(@results) . ")!\n";
}

my $basic_file = "$output_directory/no_hints.prediction.log";
if ( -s $basic_file ) {
 my $first_accuracy = `grep Accuracy $basic_file`;
 if ($first_accuracy) {
  $first_accuracy =~ /Target:\s([\.\d]+)/;
  my $first_target = $1;
  print
"Without any extrinsic evidence the target was $first_target\n$first_accuracy";
 }
}
if ( $results[0] =~ /^([^:]+)/ ) {
 my $file          = $1;
 my $best_accuracy = `grep Accuracy $file`;
 $best_accuracy =~ /Target:\s([\.\d]+)/;
 my $best_target = $1;
 print "Best config file is found in $file with target $best_target\n$best_accuracy\n";
 print "You may however want to run this command and see if there is a config file that performs better in any specific statistic\n";
 print "\tgrep Acc $output_directory/*log | sort -nk11\n";
}

################################################################################################################
################################################################################################################
sub parse_evaluation() {
 my $pred_out = shift;
 my ( $cbsn, $cbsp, $cesn, $cesp, $cgsn, $cgsp, $csmd, $ctmd ) =
   ( int(0), int(0), int(0), int(0), int(0), int(0), int(0), int(0) );

 open( PRED, $pred_out );
 while ( my $ln = <PRED> ) {
  if ( $ln =~ /nucleotide level \| +(\S+) \| +(\S+) \|$/ ) {
   ( $cbsn, $cbsp ) = ( $1, $2 );
  }
  elsif ( $ln =~ /exon level \|.*-- \| +(\S+) \| +(\S+) \|$/ ) {
   ( $cesn, $cesp ) = ( $1, $2 );
  }
  elsif ( $ln =~ /gene level \|.* \| +(\S+) \| +(\S+) \|$/ ) {
   ( $cgsn, $cgsp ) = ( $1, $2 );
  }
  elsif ( $ln =~ /TSS \|.* \|.* \|.* \| +(\S+) \|$/ ) {
   $csmd = $1;
  }
  elsif ( $ln =~ /TTS \|.* \|.* \|.* \| +(\S+) \|$/ ) {
   $ctmd = $1;
  }
 }
 close(PRED);
 $cbsn = $cbsn ? $cbsn : int(0);
 $cbsp = $cbsp ? $cbsp : int(0);
 $cesn = $cesn ? $cesn : int(0);
 $cesp = $cesp ? $cesp : int(0);
 $cgsn = $cgsn ? $cgsn : int(0);
 $cgsp = $cgsp ? $cgsp : int(0);
 $csmd = $csmd ? $csmd : int(0);
 $ctmd = $ctmd ? $ctmd : int(0);

 return ( $cbsn, $cbsp, $cesn, $cesp, $cgsn, $cgsp, $csmd, $ctmd );
}

sub check_prediction_finished(){
	my $file = shift;
	my $todelete = shift;
	my $check = `tail -n 3 $file|head -n 1`;
	if ($check !~/^# total time/){
		unlink($file);
		unlink($todelete);
	}
}

sub run_evaluation {
 my $extrinsic_file = shift;

 my ( $pred_out, $augustus_cmd);
 if ( $extrinsic_file && -s $extrinsic_file ) {
  $pred_out = $extrinsic_file;
  $pred_out =~ s/.cfg/.prediction/;
  &check_prediction_finished($pred_out,$pred_out.'.log');
  $augustus_cmd =
"$augustus_exec --genemodel=complete --species=$species --alternatives-from-evidence=true --AUGUSTUS_CONFIG_PATH=$config_path --extrinsicCfgFile=$extrinsic_file --hintsfile=$hint_file $common_parameters $optimize_gb > $pred_out 2>/dev/null";
 }
 else {
  $pred_out = $output_directory . '/no_hints.prediction';
  &check_prediction_finished($pred_out,$pred_out.'.log');
  $augustus_cmd =
"$augustus_exec --genemodel=complete --species=$species --AUGUSTUS_CONFIG_PATH=$config_path $common_parameters $optimize_gb > $pred_out 2>/dev/null";
 }

 &process_cmd($augustus_cmd) unless -s $pred_out;

 unless ( -s $pred_out.'.log' ) {
  if ( $pred_out && -s $pred_out ) {
   my @eval_results = &parse_evaluation($pred_out);
   my ($target,$results_ref) = &estimate_accuracy(@eval_results) ;
   $target = sprintf( "%.4f", $target );
   open( TLOG, '>' . $pred_out . '.log' );
   print TLOG "#Accuracy: "
     . join( ", ", @$results_ref )
     . "; Target: $target\n";
   if ($extrinsic_file) {
    open( IN, $extrinsic_file );
    while ( my $ln = <IN> ) { print TLOG $ln; }
    close IN;
   }
   close TLOG;
   return $target;
  }
  else {
   warn "Augustus failed to produce $pred_out\n.";
  }
 }
}

sub estimate_accuracy {
 my $switch; # trialling for higher specificity
 my ( $bsn, $bsp, $esn, $esp, $gsn, $gsp, $smd, $tmd ) = @_;
 my $target =int(0);
 if ( !$switch ) {
  $target = (
   3 * $bsn +
     3 * $bsp +
     4 * $esn +
     4 * $esp +
     2 * $gsn +
     2 * $gsp +
     1 * (40 / ( $smd + 40 )) +
     1 * (40 / ( $tmd + 40 ))
  ) / 20;
 }
 else {
  $target = ( 
   3 * $bsn + 
   9 * $bsp + 
   4 * $esn + 
   12 * $esp + 
   2 * $gsn + 
   6 * $gsp 
  ) / 36;
 }
 my @results = ( "Base_SeNsitivity:".$bsn, "Base_SPecificity:".$bsp, "Exon_SeNsitivity:".$esn, "Exon_SPecificity;".$esp, "Gene_SeNsitivity:".$gsn, "Gene_SPecificity:".$gsp, "UTR_5':".$smd, "UTR_3':".$tmd );

 return ($target,\@results);

}

sub printmetavalues {
 my ( $names, $values ) = @_;
 for ( my $i = 0 ; $i <= $#$names ; $i++ ) {
  print $names->[$i] . "\t" . $values->[$i] . "\n";
 }
}

sub check_options() {
 die(
"Augustus config directory not in environment (\$AUGUSTUS_CONFIG_PATH) and not specified.\n"
 ) unless $config_path && -d $config_path;
 $config_path .= '/' unless $config_path =~ /\/$/;

 $optimize_gb = shift if !$optimize_gb;
 pod2usage("Optimizing file missing\n") if ( !$optimize_gb );

 pod2usage "No --extrinsic file with metaparameters provided\n"
   unless $extrinsic && -s $extrinsic;
 pod2usage("Hint file is required!\n") if ( !$hint_file || !-s $hint_file );
 pod2usage("no species specified\n") if ( !$species );
 $common_parameters .= " --species=$species ";
 $ENV{'PATH'} .= ':' . $exec_dir if $exec_dir && -d $exec_dir;
 $common_parameters .= " --UTR=on"                         if ($utr);
 $common_parameters .= " --translation_table=$trans_table" if ($trans_table);
 $common_parameters .= " --genemodel=$genemodel"           if ($genemodel);
 $common_parameters .= " --/Constant/min_coding_len=" . $min_coding
   if ($min_coding);
 $cpus = 1 if !$cpus || $cpus < 1;
 $output_directory = "tmp_opt_extrinsic_$species";
 mkdir($output_directory) unless -d $output_directory;
 die unless -d $output_directory;

}

sub write_extrinsic_cfg() {
 my $extrinsic_ref          = shift;
 my $cfg_extra              = shift;
 my $species_extrinsic_file = $output_directory . "/extrinsic.cfg";
 if ($extrinsic_ref) {
  lock($extrinsic_iteration_counter);
  $species_extrinsic_file .= ".$extrinsic_iteration_counter";
  $extrinsic_iteration_counter++;
 }

 unless ( -s $species_extrinsic_file ) {

  open( OUT, ">$species_extrinsic_file" ) || die;
  print OUT "[SOURCES]\n";
  my $prn;

  foreach my $src ( keys %{$extrinsic_ref} ) {
   $prn .= $src . ' ' unless $src eq 'bonus' || $src eq 'malus';
  }
  chop($prn);
  $prn .= "\n";
  $prn .= $cfg_extra . "\n\n" if $cfg_extra;
  $prn .= "[GENERAL]\n";

  foreach my $feature (@feature_order) {
   my $bonus =
       $extrinsic_ref->{'bonus'}->{$feature}->{'value'}
     ? $extrinsic_ref->{'bonus'}->{$feature}->{'value'}
     : 1;
   my $malus =
       $extrinsic_ref->{'malus'}->{$feature}->{'value'}
     ? $extrinsic_ref->{'malus'}->{$feature}->{'value'}
     : 1;
   $prn .= $feature . "\t$bonus $malus ";
   $prn .= "\t" if length($feature) <= 7;
   $prn .= "\tM 1 " . $extrinsic_ref->{'M'}->{$feature}->{'value'};

   foreach my $source ( keys %{$extrinsic_ref} ) {
    next if $source eq 'bonus' || $source eq 'malus' || $source eq 'M';
    my $multiplier = $extrinsic_ref->{$source}->{$feature}->{'value'}
      || die "No bonus found for $source $feature\n";
    $prn .= "\t$source 1 $multiplier";
   }
   $prn .= "\n";
  }
  print OUT $prn . "\n";
  close OUT;
 }
 my $target = &run_evaluation($species_extrinsic_file);
 print "\r$extrinsic_iteration_counter/$todo    ";
}

sub parse_extrinsic_meta() {
 my ( %extrinsic_hash, %basic_cfg, $cfg_extra );
 my $species_meta_cfg_filename =
     $extrinsic
   ? $extrinsic
   : $species_config_dir . "/" . $species . "_metaextrinsic.cfg";
 open( IN, $species_meta_cfg_filename )
   || die("Can't open $species_meta_cfg_filename\n");
 @feature_headers = split( "\t", <IN> );

 die
   "CFG meta file $species_meta_cfg_filename has less than 18 header columns ("
   . scalar(@feature_headers) . ")\n"
   unless scalar(@feature_headers) == 18;
 chomp( $feature_headers[-1] );
 for ( my $f = 0 ; $f < (@feature_headers) ; $f++ ) {
  $feature_headers[$f] = 'exonpart'    if $feature_headers[$f] eq 'ep';
  $feature_headers[$f] = 'nonexonpart' if $feature_headers[$f] eq 'nep';
  $feature_headers[$f] = 'intronpart'  if $feature_headers[$f] eq 'ip';
 }

 while ( my $ln = <IN> ) {
  next if $ln =~ /^\s*$/ || $ln =~ /^#/;

  #add the extra data, these appear at the end of the file.
  if ( $ln =~ /^EXTRA/ ) {
   while ( my $ln2 = <IN> ) {
    next if $ln =~ /^\s*$/;
    $cfg_extra .= $ln2;
   }
   last;
  }
  chomp($ln);
  my @data = split( "\t", $ln );
  die
"CFG meta file $species_meta_cfg_filename has less than 18 columns for this line:\n$ln\n"
    unless @data == 18;
  my $source = $data[0];
  for ( my $i = 1 ; $i < 18 ; $i++ ) {
   my $score = $data[$i];
   my @parameters;
   my $basic = 1;
   if ( $score =~ /,/ ) {
    @parameters = split( ',', $score );
    $sources_that_need_to_be_checked{$source} = 1
      unless $source eq 'bonus' || $source eq 'malus';
   }
   elsif ( $score =~ /-/ ) {
    my ( $min, $max ) = split( '-', $score );
    if ( $max < $min ) {
     my $t = $max;
     $max = $min;
     $min = $t;
    }
    my $step = ( $max - $min ) / $pstep;
    for ( my $next_level = $min ; $next_level <= $max ; $next_level += $step ) {
     push( @parameters, $next_level );
    }
    $parameters[-1] = $max;
    $sources_that_need_to_be_checked{$source} = 1
      unless $source eq 'bonus' || $source eq 'malus';
   }
   else {
    $basic      = $score;
    @parameters = ($score);
   }

   $basic_cfg{$source}->{ $feature_headers[$i] }->{'value'}      = $basic;
   $extrinsic_hash{$source}->{ $feature_headers[$i] }->{'value'} = \@parameters;

   #$basic_cfg{$source}->{ $header[$i] }      = \@scores;
   #$extrinsic_hash{$source}->{ $header[$i] } = \@parameters;

  }

 }
 close IN;

 foreach my $head (@feature_headers) {
  $basic_cfg{'M'}->{$head}->{'value'} = 1e+100;
 }
 return ( \%extrinsic_hash, $cfg_extra, \%basic_cfg );
}

sub process_cmd {
 my ( $cmd, $dir ) = @_;
 print &mytime . "CMD: $cmd\n" if $verbose;
 chdir($dir) if $dir;
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  chdir($cwd) if $dir;
  die "Error, cmd died with ret $ret\n";
 }
 chdir($cwd) if $dir;
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

###
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
