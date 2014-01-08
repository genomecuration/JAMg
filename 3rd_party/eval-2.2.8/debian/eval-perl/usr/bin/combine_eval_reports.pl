#!/usr/bin/perl

use strict;
use Eval;
use Getopt::Std;

use vars qw($opt_e $opt_h $opt_s);
getopts('ehs');

my $MIN_DISP_LEN = 25;
my $MIN_VAL_LEN = 10;
my @GENERAL_SKIP = qw(Correct
		      Specificity
		      Sensitivity
		      Matched);

my %GENERAL_SKIP = (stat => {Correct => 1,
			     Specificity => 1,
			     Sensitivity => 1,
			     Ann_Matched => 1,
			     Correct_Nucleotides => 1,
			     Nucleotide_Specificity => 1,
			     Nucleotide_Sensitivity => 1,
			     Ann_Matched_Nucleotides => 1},
		    type => {Exact => 1,
			     Partial => 1,
			     Overlap => 1,
			     Nuc_Overlap => 1,
			     All_Introns => 1,
			     Start_Stop => 1,
			     Splice_Acceptor => 1,
			     Splice_Donor => 1,
			     Overlap_80p => 1});
my $Precision = "%.2f";
my %General;
my %Display;

my $EVAL = 0;
my $STATS = 1;
my $Mode = $EVAL;

my @reports;
my @report_names;
my @ann_reports;
my @ann_report_names;

my $usage = "$0 [-hes] <report 1> <report 2> [report 3] ...
This script combines several Eval reports or general statistics reports into one. 
  Options:
    -e: Eval report mode [default]
    -s: General statistics report mode.  Cannot be used with -e.
    -h: Display this help message and exit.
";

if($opt_h){
    print $usage;
    exit(0);
}
if($opt_s){
    if($opt_e){
	die "Cannot give both -e and -s flags.\n";
    }
    $Mode = $STATS;
}
elsif($opt_e){
    $Mode = $EVAL;
}

die $usage unless @ARGV >= 2;

init_hashes();

foreach my $file (@ARGV){
    my $fh;
    open($fh,$file) or die "Couldn't open $file.";
    while(my $line = <$fh>){
	if($line =~ /\*\*General Stats\*\*/){
	    <$fh>;
	    my $line = <$fh>;
	    chomp $line;
	    $line =~ /^\s+(.+)$/;
	    my $pred_text = $1;
	    my @preds = split /\t/, $pred_text;

	    my @data = readdata($fh);
	    if($Mode == $STATS){ 
		push @reports, @data;
		push @report_names, @preds;
	    }
	    elsif($Mode == $EVAL){
		push @ann_reports, $data[0];
		push @ann_report_names, $preds[0];		
				
		my $line = <$fh>;
		chomp $line;
		$line =~ /Annotation:\s+(\S+)/;
		my $ann = $1;
		my $line = <$fh>;
		chomp $line;
		$line =~ /Predictions:\s+(.+)$/;
		my $pred_text = $1;
		my @preds = split /\t/, $pred_text;
		
		my @data = readdata($fh);

		push @reports, @data;
		push @report_names, @preds;
	    }
	}
    }
    close($fh);
}

if($Mode == $STATS){
    print get_general_stats_text(\@reports,\@report_names);    
}
elsif($Mode == $EVAL){
    my $first_name = $ann_report_names[0];
    foreach my $name (@ann_report_names){
	unless($name eq $first_name){
	    print STDERR 
		"WARNING: combinding reports with different annotation set names.\n";
	    last;
	}
    }
    my @eval_reports = ($ann_reports[1]);
    push @eval_reports, @reports;
    my @eval_report_names = ($ann_report_names[1]);
    push @eval_report_names, @report_names;
    print_eval_output(\@reports,\@report_names);
    #print_eval_output(\@eval_reports,\@eval_report_names);
}

exit(0);



sub readdata{
    my ($fh) = @_;
    my @data;
    my $level = "none";
    my $type = "none";
    my $stat = "none";
    while(my $line = <$fh>){
	if($line =~ /\*\*Detailed Stats\*\*/){
	    last;
	}
	#print $line;
	chomp $line;
	if($line =~ /^(\S.*\S)\s*$/){
	    $level = $1;
	    $level =~ s/ /_/;
	}
	elsif($line =~ /^\t(\S.*\S)\s*$/){
	    $type = $1;
	    $type =~ s/ /_/;
	}
	elsif($line =~ /^\t\t(\S.*)$/){
	    $line = $1;
	    my @stats;
	    while($line =~ /^(\S+)\s(.+)$/){
		push @stats, $1;
		$line = $2;
	    }
	    $stat = $stats[0];
	    for(my $i = 1;$i <= $#stats;$i++){
		$stat .= "_".$stats[$i];
	    }
	    $line =~ /^\s*(.+)$/;
	    my $val_text = $1;
	    my @vals = split /\s+/, $val_text;
	    for(my $i = 0;$i <= $#vals;$i++){
		unless(defined($data[$i])){
		    $data[$i] = {};
		}
		unless(defined($data[$i]{$level})){
		    $data[$i]{$level} = {}
		}
		unless(defined($data[$i]{$level}{$type})){
		    $data[$i]{$level}{$type} = {};
		}
		#print "load data[$i]{$level}{$type}{$stat} = $vals[$i];\n";
		$data[$i]{$level}{$type}{$stat} = $vals[$i];
	    }
	}
    }
    
    return @data;
}



sub print_eval_output{
    my ($data,$names) = @_;
    my $report = "";
    # add summary stats
    $report .= "\n**Summary Stats**\n";    
    $report .= "Annotation:\t$$names[0]\n";
    $report .= pad_string("Predictions:",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$names[$i]);
    }
    $report .= "\n";
    $report .= "\n".pad_string("Gene Sensitivity",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$data[$i]{Gene}{All}{Consistent_Sensitivity});
    }
    $report .= "\n".pad_string("Gene Specificity",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$data[$i]{Gene}{All}{Consistent_Specificity});
    }
    $report .= "\n".pad_string("Transcript Sensitivity",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$data[$i]{Transcript}{All}{Consistent_Sensitivity});
    }
    $report .= "\n".pad_string("Transcript Specificity",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$data[$i]{Transcript}{All}{Consistent_Specificity});
    }
    $report .= "\n".pad_string("Exon Sensitivity",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$data[$i]{Exon}{All}{Correct_Sensitivity});
    }
    $report .= "\n".pad_string("Exon Specificity",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$data[$i]{Exon}{All}{Correct_Specificity});
    }
    $report .= "\n".pad_string("Nucleotide Sensitivity",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$data[$i]{Nuc}{All}{Correct_Sensitivity});
    }
    $report .= "\n".pad_string("Nucleotide Specificity",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= "\t".pad_string($$data[$i]{Nuc}{All}{Correct_Specificity});
    }
    # add general stats
    $report .= get_general_stats_text($data,$names);
    # add detailed stats
    $report .= "\n\n**Detailed Stats**\n";
    $report .= "Annotation:\t$$names[0]\n";
    $report .= "Predictions:\t".pad_string("",$MIN_DISP_LEN);
    for(my $i = 1; $i <= $#$data;$i++){
	$report .= pad_string($$names[$i]);
	if($i == $#$data){
	    $report .= "\n";
	}
	else{
	    $report .= "\t";
	}
    }
    if($#$data >= 0){
	my %order = Eval::get_list_struct();	
	foreach my $level (@{$order{Levels}}){
	    my $level_disp = $level;
	    $level_disp =~ s/_/ /g;
	    $report .= "$level_disp\n";
	    foreach my $type (@{$order{$level}{Type}}){
		if($Display{$level}{Type}{$type} == 1){
		    my $type_disp = $type;
		    $type_disp =~ s/_/ /g;
		    $report .= "\t$type_disp\n";
		    foreach my $stat (@{$order{$level}{Stat}}){
			if($Display{$level}{Stat}{$stat} == 1){			
			    my $stat_disp = $stat;
			    $stat_disp =~ s/_/ /g;
			    $report .= "\t\t".pad_string($stat_disp,$MIN_DISP_LEN);
			    for(my $i = 1; $i <= $#$data;$i++){
				$report .= "\t".
				    pad_string($$data[$i]{$level}{$type}{$stat});
			    }
			    $report .= "\n";
			}
		    }
		}
	    }
	}
    }
    print $report;
}

sub get_general_stats_text{
    my ($data,$names) = @_;
    my $report = "";
    $report .= "\n\n**General Stats**\n";
    $report .= "Predictions:\n\t\t";
    for(my $i = 0;$i < $MIN_DISP_LEN;$i++){
	$report .= " ";
    }
    for(my $i = 0; $i <= $#$data;$i++){
	my $val = "$$names[$i]";
	my $len = length($val);
	for(my $i = $len;$i < $MIN_VAL_LEN;$i++){
	    $val .= " ";
	}
	$report .= "\t$val";
    }
    $report .= "\n";
    my %order = Eval::get_general_list_struct();
    if($#$data >= 0){
	foreach my $level (@{$order{Levels}}){
	    my $level_disp = $level;
	    $level_disp =~ s/_/ /g;
	    $report .= "$level_disp\n";
	    foreach my $type (@{$order{$level}{Type}}){		
		if(($General{$level}{Type}{$type} == 1) &&
		   ($Display{$level}{Type}{$type} == 1)){
		    my $type_disp = $type;
		    $type_disp =~ s/_/ /g;
		    my $len = length($type_disp);
		    for(my $i = $len;$i < $MIN_DISP_LEN;$i++){
			$type_disp .= " ";
		    }
		    $report .= "\t$type_disp\n";
		    foreach my $stat (@{$order{$level}{Stat}}){
			if(($General{$level}{Stat}{$stat} == 1) &&
			   ($Display{$level}{Stat}{$stat} == 1) &&
			   ($stat !~ /Ann_/)){
			    my $stat_disp = $stat;
			    $stat_disp =~ s/_/ /g;
			    $len = length($stat_disp);
			    for(my $i = $len;$i < $MIN_DISP_LEN;$i++){
				$stat_disp .= " ";
			    }
			    $report .= "\t\t$stat_disp";
			    for(my $i = 0; $i <= $#$data;$i++){
				my $val = $$data[$i]{$level}{$type}{$stat};
				$len = length($val);
				for(my $i = $len;$i < $MIN_VAL_LEN;$i++){
				    $val .= " ";
				}
				$report .= "\t$val";
			    }
			    $report .= "\n";
			}
		    }
		}
	    }
	}
    }
	return $report;
}

sub pad_string{
    my ($string,$min_len) = @_;
    unless(defined($min_len)){
	$min_len = $MIN_VAL_LEN;
    }
    my $len = length($string);
    for(my $i = $len;$i < $min_len;$i++){
	$string .= " ";
    }
    return $string;
}

# init the General and Display hashes
sub init_hashes{
    my %order = Eval::get_list_struct();
    foreach my $level (@{$order{Levels}}){
	foreach my $type (@{$order{$level}{Type}}){
	    $Display{$level}{Type}{$type} = 1;
	    $General{$level}{Type}{$type} = 1;
	}
	foreach my $stat (@{$order{$level}{Stat}}){
	    $Display{$level}{Stat}{$stat} = 1;
	    $General{$level}{Stat}{$stat} = 1;
	    foreach my $skip (@GENERAL_SKIP){
		if($stat =~ /$skip/){
		    $General{$level}{Stat}{$stat} = 0;
		}
	    }
	}
    }
}

sub adjust_data_precision{
    my ($data) = @_;
    my %order = Eval::get_list_struct();
    for(my $i = 0; $i <= $#$data;$i++){
	foreach my $level (@{$order{Levels}}){
	    foreach my $type (@{$order{$level}{Type}}){
		foreach my $stat (@{$order{$level}{Stat}}){
		    $$data[$i]{$level}{$type}{$stat} = 
			sprintf $Precision, $$data[$i]{$level}{$type}{$stat};
		    if(($stat =~ /Sensitivity/) ||
		       ($stat =~ /Specificity/)){
			$$data[$i]{$level}{$type}{$stat} .= "%";
		    }
		    else{
			$$data[$i]{$level}{$type}{$stat} .= " ";
		    }
		}
	    }
	}
    }    
}
