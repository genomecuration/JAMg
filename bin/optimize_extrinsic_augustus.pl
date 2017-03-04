#!/usr/bin/env perl

=pod

=head1 NAME

Optimize extrinsic configuration for Augustus. Based on optimize_augustus by Mario Stanke, 23.04.2007 
Modified by Jeremy Zucker <djinnome@gmail.com> to include local malus. 
=head1 USAGE 

Mandatory parameters:

 --species                  prefix of the species name
 --optimize                 genbank file for training with bona fide gene structures
 --hints                    hints file to use
 --extrinsic or --metapars  File with the names and their ranges of the meta parameters that are subject to optimization (default: configs/extrinsic_meta.cfg)

Optional parameters:

 --config_path              Specify the config directory d if not set as environment variable
 --aug_exec_dir             Path to augustus and etraining executable. If not specified it must be in $PATH environment variable
 --cpus                     The number of CPUs to use (default: 1)
 --onlytrain                an optional genbank file that is used in addition to train.gb but only for etrain not for intermediate evaluation of accuracy. These genes may e.g. be incomplete.
 --utr                      Switch on UTR. Recommended to avoid noncoding exons being wrongly predicted as coding.
 --pstep                    For integer and floating parameters start with p tests equidistributed in the allowed range of values (default: 5)
 --min_coding               Minimum coding length
 --genemodel                Augustus genemodel parameter
 --notrain                  The training step (etraining) is omitted. Use this if the HMM is already trained. This program does not optimize the HMM
 --specificity              Weight Accuracy statistics in favour of specificity over sensitivity (only affects 'Accuracy' stat in .log files)

=head1 AUTHORS

 Alexie Papanicolaou

        CSIRO Ecosystem Sciences
        alexie@butterflybase.org

 Jeremy Zucker

       Orion Genomics
       djinnome@gmail.com

Based on work by Mario Stanke, 23.04.2007

=head1 DISCLAIMER & LICENSE

Copyright 2012-2014 the Commonwealth Scientific and Industrial Research Organization. 
See LICENSE file for license info
It is provided "as is" without warranty of any kind.

=cut

use strict;
use warnings;
use Data::Dumper;
use File::Basename;
use File::Copy;
use Pod::Usage;
use Getopt::Long;
use Time::localtime;
use IO::File;
use threads;
use threads::shared;
use FindBin qw($RealBin);
use Digest::MD5 qw(md5_hex);
use Storable 'dclone';
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";
use Thread_helper;

my (
     $species,     $optimize_gb,      $extrinsic,
     $rounds,      $training_set,     $pstep,
     $config_path, $cpus,             $utr,
     $exec_dir,    $output_directory, $species_config_dir,
     $notrain,     $trans_table,      $genemodel,
     $min_coding,  $hint_file,        $verbose, $prefer_specificity,
     $debug
);
$|           = 1;
$rounds      = 1;
$pstep       = 5;
$cpus        = 1;
$exec_dir    = $ENV{'AUGUSTUS_PATH'} . '/bin' if $ENV{'AUGUSTUS_PATH'};
$config_path = $ENV{'AUGUSTUS_PATH'} . '/config' if $ENV{'AUGUSTUS_PATH'};
$config_path = $ENV{'AUGUSTUS_CONFIG_PATH'} if $ENV{'AUGUSTUS_CONFIG_PATH'};
my $extrinsic_iteration_counter : shared;
my %cfg_track : shared; # hack to prevent same config from being rerun

my ($common_parameters) = '';

pod2usage $! unless &GetOptions(
             'species:s'              => \$species,
             'optimize:s'             => \$optimize_gb,
             'metapars|extrinsic:s'   => \$extrinsic,
             'rounds:i'               => \$rounds,
             'onlytrain|training_set:s'         => \$training_set,
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
             'verbose'                => \$verbose,
             'specificity'            => \$prefer_specificity,
             'debug'                  => \$debug,
);

# globals
&check_options();
my ( $augustus_exec, $etrain_exec ) = &check_program( 'augustus', 'etraining' );
my $cwd = `pwd`;
chomp($cwd);
my $be_silent =
" --/augustus/verbosity=0 --/ExonModel/verbosity=0 --/IGenicModel/verbosity=0 --/IntronModel/verbosity=0 --/UtrModel/verbosity=0 --/genbank/verbosity=0 "
;
my @feature_order =
  qw/start stop tss tts ass dss exonpart exon intronpart intron CDSpart CDS UTRpart UTR irpart nonexonpart genicpart/;

print "Started: " . &mytime . "\n";

# we only need to train once
system("cat $optimize_gb $training_set > train.set 2>/dev/null") if $training_set ;
&process_cmd("$etrain_exec --species=$species --AUGUSTUS_CONFIG_PATH=$config_path $common_parameters $be_silent train.set  >/dev/null 2>/dev/null ") unless $notrain;
unlink("train.set");

my ($sources_that_need_to_be_checked_ref,$feature_headers_ref, $metaextrinsic_hash_ref, $cfg_extra, $basic_cfg_ref ) =  &parse_extrinsic_meta();
my %sources_that_need_to_be_checked = %{$sources_that_need_to_be_checked_ref};
my @feature_headers = @$feature_headers_ref;

#######################################################################################
# initialize and first test
#######################################################################################

$extrinsic_iteration_counter = int(0);

my $thread_helper = new Thread_helper($cpus);
my $empty_eval_thread        = threads->create('run_evaluation');
$thread_helper->add_thread($empty_eval_thread);
my $todo = int(0);

foreach my $src ( keys %sources_that_need_to_be_checked ) {
 foreach my $feat (@feature_headers) {
  next if $feat eq 'feature';
  my $bonus_range = $metaextrinsic_hash_ref->{'bonus'}->{$feat}->{'value'};
  my $malus_range = $metaextrinsic_hash_ref->{'malus'}->{$feat}->{'value'};
  my $local_malus_range = $metaextrinsic_hash_ref->{'local malus'}->{$feat}->{'value'};
  foreach my $bonus ( @{$bonus_range} ) {
   foreach my $malus ( @{$malus_range} ) {
     foreach my $local_malus (@{$local_malus_range} ) {
       my $parameter_range = $metaextrinsic_hash_ref->{$src}->{$feat}->{'value'};
       if ( scalar(@$parameter_range) > 1 ) {
	 foreach my $value (@$parameter_range) {
	  $todo++;
  	  $sources_that_need_to_be_checked{$src}++;
	 }
       }
     }
   }
  }
 }
}

print "Will search for $todo parameters with $cpus CPUs\n";

my %file_handles;
unlink("$output_directory/track_evidence_lines.txt") if -s "$output_directory/track_evidence_lines.txt";

warn Dumper \%sources_that_need_to_be_checked;

#NB we really expect one line of evidence to be checked at a time...!
foreach my $src ( keys %sources_that_need_to_be_checked ) {
 print "\nProcessing $src\n";
 open (my $track_fh,">$output_directory/track_evidence_lines.$src.txt");
 $file_handles{$src} = $track_fh;
 print $track_fh "$src\t";
 foreach my $feat (@feature_headers) {
  next if $feat eq 'feature';
  my $bonus_range = $metaextrinsic_hash_ref->{'bonus'}->{$feat}->{'value'};
  my $malus_range = $metaextrinsic_hash_ref->{'malus'}->{$feat}->{'value'};
  my $local_malus_range = $metaextrinsic_hash_ref->{'local malus'}->{$feat}->{'value'};

  # for every bonus range
  foreach my $bonus ( @{$bonus_range} ) {

   # for every malus range
   foreach my $malus ( @{$malus_range} ) {

     # for every local malus range
     foreach my $local_malus ( @{$local_malus_range} ) {
       my $parameter_range = $metaextrinsic_hash_ref->{$src}->{$feat}->{'value'};
       if ( scalar(@$parameter_range) > 1 ) {
	 foreach my $value (@$parameter_range) {
	   # checks every value against the basic_cfg
	   my $cfg_hash = dclone $basic_cfg_ref;
	   $cfg_hash->{$src}->{$feat}->{'value'}    = $value;
	   $cfg_hash->{'bonus'}->{$feat}->{'value'} = $bonus;
	   $cfg_hash->{'malus'}->{$feat}->{'value'} = $malus;
	   $cfg_hash->{'local malus'}->{$feat}->{'value'} = $local_malus;
	   my $cfg_serial = Dumper $cfg_hash;
	   if (!$cfg_track{$cfg_serial}){
	     my $thread = threads->create( 'write_extrinsic_cfg', $cfg_hash, $cfg_extra,$track_fh,$src );
	     if ($debug){
		$thread_helper->add_thread($thread,1);
	     }else{
		$thread_helper->add_thread($thread,30);
	    }
	   }
	 }
       }
     }
   }
  }
 }
}


my $no_hint_accuracy = &get_accuracy_without_xtn($empty_eval_thread);
$thread_helper->wait_for_all_threads_to_complete();
sleep(10) unless $debug;
&shutdown_fh();

my @failed_threads = $thread_helper->get_failed_threads();
if (@failed_threads) {
  warn "Error, " . scalar(@failed_threads) . " threads failed.\nThis is probably if Augustus didn't like some combinations of search parameters but performed ok in the majority of searches. In that case that's ok to proceed or you can opt to retry (I won't delete existing). Any other failure is indicative of something else wrong.\n\n";
}

# add here a routine to check between columns
#die Dumper $no_hint_accuracy if $debug;

# this is post-processing
### FINISHED
my ($result_hash,@ordered_results) = &post_process;
#die Dumper $result_hash;



################################################################################################################
################################################################################################################
sub parse_best_predictions(){
 my $tracking = shift;
 my $return_output = "Evidence type\ttop accuracy\ttop file\n";
 open (IN,$tracking);
 while (my $ln=<IN>){
 	chomp($ln);
	next if $ln=~/^\s*$/;
	my @log_files = split("\t",$ln);
	my $source_evidence = shift(@log_files);
	my $top_accuracy = int(0);
 	my $top_accuracy_file;
	foreach my $log (@log_files){
		next unless -s $log;
		open (LOG,$log);
		while (my $ln2 = <LOG>){
			if ($ln2 =~/Accuracy:\s+([\d\.]+);/){
				if ($top_accuracy < $1){
					$top_accuracy = $1;
					$top_accuracy_file = $log;
				}
				last;
			}
		}
		close LOG;
	}
	$return_output .= "$source_evidence\t$top_accuracy\t$top_accuracy_file\n" if $top_accuracy > 0;
 }
 close IN;
 open (OUT,">$tracking.best");
 print OUT $return_output;
 close OUT;
 return $return_output;
}

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
 $cbsn = $cbsn ? sprintf( "%.3f",$cbsn) : int(0);
 $cbsp = $cbsp ? sprintf( "%.3f",$cbsp) : int(0);
 $cesn = $cesn ? sprintf( "%.3f",$cesn) : int(0);
 $cesp = $cesp ? sprintf( "%.3f",$cesp) : int(0);
 $cgsn = $cgsn ? sprintf( "%.3f",$cgsn) : int(0);
 $cgsp = $cgsp ? sprintf( "%.3f",$cgsp) : int(0);
 $csmd = $csmd ? $csmd : int(0);
 $ctmd = $ctmd ? $ctmd : int(0);

 return ( $cbsn, $cbsp, $cesn, $cesp, $cgsn, $cgsp, $csmd, $ctmd );
}

sub check_prediction_finished(){
	my $file = shift;
	my @todeletes = @_;
	my $check = `tail -n 3 $file|head -n 1 `;
	if ($check !~/^# total time/){
		unlink($file);
		foreach my $todelete (@todeletes){
			warn "File $todelete had not finished. Deleting\n" if -s $todelete; 
			unlink($todelete);
		}
	}
}

sub run_evaluation {
 my $extrinsic_file = shift;
 my $track_being_optimized = shift;
 my ( $pred_out, $augustus_cmd);
 if ( $extrinsic_file && -s $extrinsic_file && $track_being_optimized) {
  $pred_out = $extrinsic_file;
  $pred_out =~ s/.cfg/.prediction/;
  my $tested_hint_file = $pred_out.'.hints';
  &check_prediction_finished($pred_out,$pred_out.'.log') if -s $pred_out;
  system("grep 'src=$track_being_optimized;' $hint_file > $tested_hint_file") unless -s $tested_hint_file;
  if (-s $tested_hint_file){
	  $augustus_cmd ="$augustus_exec --genemodel=complete --species=$species --alternatives-from-evidence=true --AUGUSTUS_CONFIG_PATH=$config_path --extrinsicCfgFile=$extrinsic_file --hintsfile=$tested_hint_file $common_parameters $optimize_gb >> $pred_out 2>/dev/null";
  }else{
	unlink($tested_hint_file);
	return;
  }
 }
 else {
  $pred_out = $output_directory . '/no_hints.prediction';
  &check_prediction_finished($pred_out,$pred_out.'.log') if -s $pred_out;
  $augustus_cmd ="$augustus_exec --genemodel=complete --species=$species --AUGUSTUS_CONFIG_PATH=$config_path $common_parameters $optimize_gb >> $pred_out 2>/dev/null";
 }
 unless (-s $pred_out){
	open (OUT,">$pred_out");
	print OUT "# $augustus_cmd \n";
	close OUT;
	&process_cmd($augustus_cmd);
 }
 unless ( -s $pred_out.'.log' ) {
  if ( $pred_out && -s $pred_out ) {
   my @eval_results = &parse_evaluation($pred_out);
   my ($accuracy,$results_ref) = &estimate_accuracy(@eval_results) ;
   open( TLOG, '>' . $pred_out . '.log' );
   print TLOG "#Accuracy: $accuracy; "
     . join( ", ", @$results_ref )
     . "\n";
   if ($extrinsic_file && -s $extrinsic_file) {
    open( IN, $extrinsic_file );
    while ( my $ln = <IN> ) { print TLOG $ln; }
    close IN;
   }
   close TLOG;
   return $accuracy;
  }
  else {
   warn "Augustus failed to produce $pred_out\n";
  }
 }
}

sub estimate_accuracy {
 my ( $bsn, $bsp, $esn, $esp, $gsn, $gsp, $smd, $tmd ) = @_;
 my $accuracy =int(0);
 if ( !$prefer_specificity ) {
  if ($utr){
    $accuracy = (
     3 * $bsn +
     3 * $bsp +
     4 * $esn +
     4 * $esp +
     2 * $gsn +
     2 * $gsp +
     1 * (40 / ( $smd + 40 )) +
     1 * (40 / ( $tmd + 40 ))
    ) / 20;
  }else{
    $accuracy = (
     3 * $bsn +
     3 * $bsp +
     4 * $esn +
     4 * $esp +
     2 * $gsn +
     2 * $gsp
    ) / 18;
  }
 }
 else {
  if ($utr){
    $accuracy = ( 
     3 * $bsn + 
     9 * $bsp + 
     4 * $esn + 
     12 * $esp + 
     2 * $gsn + 
     6 * $gsp +
     1 * (40 / ( $smd + 40 )) +
     1 * (40 / ( $tmd + 40 ))
    ) / 38;
  }else{
    $accuracy = ( 
     3 * $bsn + 
     9 * $bsp + 
     4 * $esn + 
     12 * $esp + 
     2 * $gsn + 
     6 * $gsp 
    ) / 36;
  }
 }

 my @results = ( "Base_SeNsitivity: $bsn", "Base_SPecificity: $bsp", "Exon_SeNsitivity: $esn", "Exon_SPecificity: $esp", "Gene_SeNsitivity: $gsn", "Gene_SPecificity: $gsp", "UTR_5_median_diff: $smd", "UTR_3_median_diff: $tmd" );
 $accuracy = sprintf( "%.4f", $accuracy );
 return ($accuracy,\@results);
}

sub printmetavalues {
 my ( $names, $values ) = @_;
 for ( my $i = 0 ; $i <= $#$names ; $i++ ) {
  print $names->[$i] . "\t" . $values->[$i] . "\n";
 }
}

sub check_options() {
 die "-onlytrain file not found\n" if $training_set && !-s $training_set;
 die "-opt file not found\n" if $optimize_gb && !-s $optimize_gb; 

 $exec_dir = "$RealBin/../3rd_party/augustus/bin/" if !$exec_dir;
 $exec_dir.='/' if $exec_dir!~/\/$/;
 $ENV{'PATH'} .= ':' . $exec_dir if $exec_dir && -d $exec_dir;
 die("Augustus exec directory not in environment (\$AUGUSTUS_PATH) and not specified.\n") unless $exec_dir && -d $exec_dir;
 $config_path = "$RealBin/../3rd_party/augustus/config/" if !$config_path;
 $config_path .= '/' unless $config_path =~ /\/$/;
 die("Augustus config directory not in environment (\$AUGUSTUS_CONFIG_PATH) and not specified.\n") unless $config_path && -d $config_path;

 $optimize_gb = shift if !$optimize_gb;
 pod2usage("Optimisation file missing\n") if ( !$optimize_gb );
 $extrinsic = "$RealBin/../configs/extrinsic_meta.cfg" unless $extrinsic;
 pod2usage "No --extrinsic file with metaparameters provided\n"
   unless $extrinsic && -s $extrinsic;
 pod2usage("Hint file is required!\n") if ( !$hint_file || !-s $hint_file );
 pod2usage("no species specified\n") if ( !$species );
 $species=~s/\s+/_/g;
 unless (-d $config_path.'species/'.$species){
        mkdir ($config_path.'species/'.$species);
        die "Cannot create $config_path/species/$species\n" unless -d $config_path.'species/'.$species;
        &copy_search_replace_file($config_path.'species/generic/generic_parameters.cfg',$config_path.'species/'.$species.'/'.$species.'_parameters.cfg','generic_',$species.'_');
        copy($config_path.'species/generic/generic_weightmatrix.txt',$config_path.'species/'.$species.'/'.$species.'_weightmatrix.txt');
        copy($config_path.'species/generic/generic_exon_probs.pbl',$config_path.'species/'.$species.'/'.$species.'_exon_probs.pbl');
        copy($config_path.'species/generic/generic_igenic_probs.pbl',$config_path.'species/'.$species.'/'.$species.'_igenic_probs.pbl');
        copy($config_path.'species/generic/generic_intron_probs.pbl',$config_path.'species/'.$species.'/'.$species.'_intron_probs.pbl');
	copy($extrinsic,$config_path.'species/'.$species.'/'.$species.'_extrinsic_metapars.cfg');
 }
 $common_parameters .= " --species=$species ";
 $common_parameters .= " --UTR=on"                         if ($utr);
 $common_parameters .= " --translation_table=$trans_table" if ($trans_table);
 $common_parameters .= " --genemodel=$genemodel"           if ($genemodel);
 $common_parameters .= " --/Constant/min_coding_len=" . $min_coding
   if ($min_coding);
 $cpus = 1 if !$cpus || $cpus < 1;
 
 $output_directory = fileparse($extrinsic, (".cfg")); # "tmp_opt_extrinsic_$species";
 mkdir($output_directory) unless -d $output_directory;
 die "Cannot create $output_directory" unless -d $output_directory;

}

sub write_extrinsic_cfg() {
 my ($extrinsic_ref, $cfg_extra, $track_fh, $track_being_tested) = @_;

 my $cfg_serial = Dumper $extrinsic_ref;
 $cfg_track{$cfg_serial}=-1;

 my $main_species_extrinsic_file = $output_directory . "/extrinsic.cfg";
 my $local_species_extrinsic_file = $main_species_extrinsic_file;
  my $prn = "[SOURCES]\nM";
  foreach my $src ( keys %{$extrinsic_ref} ) {
   $prn .= " $src" if $src eq $track_being_tested;
  }
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
   my $local_malus = 
     $extrinsic_ref->{'local malus'}->{$feature}->{'value'}
       ? $extrinsic_ref->{'local malus'}->{$feature}->{'value'}
	 :1;

   $prn .= $feature . "\t$bonus $malus $local_malus";
   $prn .= "\t" if length($feature) <= 7;
   $prn .= "\tM 1 " . $extrinsic_ref->{'M'}->{$feature}->{'value'};

   foreach my $source ( keys %{$extrinsic_ref} ) {
      next unless $source eq $track_being_tested;
      my $multiplier = $extrinsic_ref->{$source}->{$feature}->{'value'} || die "No bonus found for $source $feature\n";
      $prn .= "\t$source 1 $multiplier";
   }
   $prn .= "\n";
  }
  $prn .= "\n";
  
  my $check_md5 =  md5_hex($prn);
  if ($check_md5){
	#lock needs a scope
	  lock($extrinsic_iteration_counter);
	  $extrinsic_iteration_counter++;
	  $local_species_extrinsic_file .= ".$check_md5";
  }
  unless (-s $local_species_extrinsic_file){
    open( OUT, ">$local_species_extrinsic_file" ) || die;
    print OUT $prn . "\n";
    close OUT;
  }
  my $pred_out = $local_species_extrinsic_file;
  $pred_out =~ s/.cfg/.prediction/;
  &check_prediction_finished($pred_out,$pred_out.'.log',$local_species_extrinsic_file) if -s $pred_out;
  print $track_fh " $pred_out.log";
  my $accuracy = &run_evaluation($local_species_extrinsic_file,$track_being_tested);
  if ($accuracy){
	  lock(%cfg_track);
	  $cfg_track{$cfg_serial}=$accuracy;
  }

  print STDERR "\r$extrinsic_iteration_counter/$todo    ";
}

sub parse_extrinsic_meta() {
 my ( %extrinsic_hash, %basic_cfg, $cfg_extra );
 my (%sources_that_need_to_be_checked,@feature_headers);
 my $species_meta_cfg_filename =
     $extrinsic
   ? $extrinsic
   : $species_config_dir . "/" . $species . "_metaextrinsic.cfg";
 open( IN, $species_meta_cfg_filename ) || die("Can't open $species_meta_cfg_filename\n");
 @feature_headers = split( "\t", <IN> );

 die "CFG meta file $species_meta_cfg_filename has less than 18 header columns ("
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
  die "CFG meta file $species_meta_cfg_filename has less than 18 columns for this line:\n$ln\n"
    unless @data == 18;
  my $source = $data[0];
  for ( my $i = 1 ; $i < 18 ; $i++ ) {
   my $score = $data[$i];
   my @parameters;
   my $basic = 1;
   if ( $score =~ /,/ ) {
    @parameters = split( ',', $score );
    $sources_that_need_to_be_checked{$source} = 1
      unless $source eq 'bonus' || $source eq 'malus' || $source eq 'local malus' || $source eq 'M';
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
      unless $source eq 'bonus' || $source eq 'malus' || $source eq 'local malus' || $source eq 'M';
   }
   else {
    $basic      = $score;
    @parameters = ($score);
   }

   $basic_cfg{$source}->{ $feature_headers[$i] }->{'value'}      = $basic;
   $extrinsic_hash{$source}->{ $feature_headers[$i] }->{'value'} = \@parameters;

  }

 }
 close IN;
 foreach my $head (@feature_headers) {
  $basic_cfg{'M'}->{$head}->{'value'} = 1e+100;
 }
 if (!exists $basic_cfg{'local malus'}) {
   my $basic = 1;
   my @parameters = ($basic);
   foreach my $head (@feature_headers) {
     
        $basic_cfg{'local malus'}->{$head}->{'value'} = $basic;
	$extrinsic_hash{'local malus'}->{ $head }->{'value'} = \@parameters;
      }
 }
 return (\%sources_that_need_to_be_checked,\@feature_headers, \%extrinsic_hash, $cfg_extra, \%basic_cfg );
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

sub copy_search_replace_file(){
	my ($filein,$fileout,$search,$replace,$noglobal) = @_;
	return unless $filein && -s $filein && $fileout;
	open (IN,$filein);
	open (OUT,">$fileout");
	while (my $ln=<IN>){
		if ($search){
			$replace = '' if !$replace; #delete
			$ln=~s/$search/$replace/ if $noglobal;
			$ln=~s/$search/$replace/g if !$noglobal;
			print OUT $ln;
		}
	}
	close IN;
	close OUT;
	
}


sub get_accuracy_without_xtn(){
  my $thread = shift;
  $thread->join() if $thread->is_running();
  if (my $error = $thread->error()) {
     die "ERROR, Failed to finish evaluation without hints.\n";
  }
  my $pred_out = $output_directory . '/no_hints.prediction';
  &check_prediction_finished($pred_out,$pred_out.'.log');
  die "Failed to finish evaluation without hints.\n" unless -s $pred_out;
  my @eval_results = &parse_evaluation($pred_out);
  my ($accuracy,$results_ref) = &estimate_accuracy(@eval_results) ;
  return $accuracy;
}


sub shutdown_fh(){
	foreach my $src (sort keys %sources_that_need_to_be_checked ) {
	 my $track_fh = $file_handles{$src};
	 print $track_fh "\n";
	 close $track_fh;
	 system("cat $output_directory/track_evidence_lines.$src.txt >> $output_directory/track_evidence_lines.txt");
	}
}

sub post_process(){
	my (%result_hash,%better_that_nohints);
	my $no_hint_file = "$output_directory/no_hints.prediction.log";
	my $no_hint_log = `grep Accuracy $no_hint_file`;
	print "Without any extrinsic evidence the accuracy was $no_hint_accuracy\n$no_hint_log";
	my @results = `grep -H Accuracy $output_directory/*log | sort -rnk 2`; # sorting based on accuracy so no glob
	if ( !$results[0] ) { die "No output files produced. This is weird, maybe you requested it to stop or something else happened? Report please!\n";}
	else {print "Done and all good!\n";}
	for (my $i=0;$i<scalar(@results);$i++){
		if ($results[$i] =~ /^([^:]+)/ ) {
			my $f = $1;
			if ($f && -s $f){
				$results[$i] = $f;
				$result_hash{$f}{'acc_line'} = `head -n 1 $f`;
				$result_hash{$f}{'src_line'} = `head -n 3 $f|tail -n 1`;
				if ($result_hash{$f}{'acc_line'} =~ /Accuracy:\s([\.\d]+)/){
					$result_hash{$f}{'value'} = $1;
					chomp($result_hash{$f}{'acc_line'});
					chomp($result_hash{$f}{'src_line'});
					$result_hash{$f}{'src_line'} =~s/^M\s+//;
					if ($result_hash{$f}{'value'} > $no_hint_accuracy){
						#already sorted
						push(@{$better_that_nohints{$result_hash{$f}}{'all'}{'src_line'}}, $result_hash{$f}{'value'});
						$better_that_nohints{$result_hash{$f}}{'best'}{'src_line'}  = $result_hash{$f}{'value'} if !$better_that_nohints{$result_hash{$f}}{'best'}{'src_line'}  || $better_that_nohints{$result_hash{$f}}{'best'}{'src_line'} < $result_hash{$f}{'value'};
						#print "File $f tested for ".$result_hash{$f}{'src_line'}." and has a better than base ($no_hint_accuracy) accuracy of ".$result_hash{$f}{'value'}."\n";
					}
				}else{
					warn "File $f failed to produce any accuracy statistics\n";
				}
			}
		}
	}
	my $best_log = $results[0];
	my $best_accuracy = $result_hash{$best_log}{'value'};
	my $best_line = $result_hash{$best_log}{'acc_line'};
	print "Best config file (across all evidence tracks) is found in $best_log with accuracy $best_accuracy\n$best_line\n";
	print "You could run this command and see if there is a config file that performs better in any specific statistic:\n";
	print "\$  grep Acc $output_directory/*log | sort -nk 2\n\n";
	my $best_output = &parse_best_predictions("$output_directory/track_evidence_lines.txt");
	if (-s "$output_directory/track_evidence_lines.txt"){
	  print "Remember that this is for a single line of evidence at a time. This is the best file for each evidence you have provided (cf. $output_directory/track_evidence_lines.*):\n$best_output\n";
	}

	return (\%result_hash,@results,\%better_that_nohints);
}

