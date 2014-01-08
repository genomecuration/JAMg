#!/usr/bin/perl
use strict;
use lib '.';
use Eval;
use GTF;
use Getopt::Std;
use vars qw($opt_g $opt_h $opt_v $opt_q $opt_A);
getopts('ghvqA');
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
my $usage = "$0 <annotation list> <prediction list 1> [prediction list 2] ...
Run the evaluation code in text mode.
Options:
  -g: Input files are gtf not lists
  -q: Quick load the gtf file.  Do not check them for errors.
  -A: Do not evaluate for alternative splicing events. (Faster)
  -v: Verbose mode
  -h: Display this help message and exit
";
if($opt_h){
    print $usage;
    exit(0);
}
my $Use_ase = 1;
if($opt_A){
    $Use_ase = 0;
}
die $usage unless (@ARGV >= 2);
print STDERR "@ARGV\n";
system("date");
my $verbose = $opt_v;
my $Quick_Load = $opt_q;
init_hashes();
my @gtfs;
my @names;
if($verbose){
    print STDERR "loading...";
}
my $gtf_mode = 1;
foreach my $arg (@ARGV){
    unless($arg =~ /\.g[tf]f$/){
	$gtf_mode = 0;
    }
}
if($gtf_mode){
    $opt_g = 1;
}
foreach my $arg (@ARGV){
    my @this_list;
    if($opt_g){
	push @this_list, GTF::new({gtf_filename => $arg,
				   no_check => $Quick_Load,
                   mark_ase => $Use_ase});
    }
    else{
	open(IN,$arg) or die "Couldn't open $arg";
	while(my $in = <IN>){
	    chomp $in;
	    my @data = split /\t/, $in;
	    $in = $data[0];
	    if($in =~ /\S/){
		push @this_list, GTF::new({gtf_filename => $in,
					   no_check => $Quick_Load,
                       mark_ase => $Use_ase});
	    }
	}
    }
    push @gtfs, \@this_list;
    push @names, $arg;
}
my $start_time = time;
if($verbose){
    print STDERR "done\n";
    print STDERR "running...";
}
my $ann = $gtfs[0];
my @preds;
for(my $i = 1;$i <= $#gtfs;$i++){
    push @preds, $gtfs[$i];
}
my @data = Eval::evaluate($ann,\@preds);
if($verbose){
    my $end_time = time;
    my $total_time = $end_time - $start_time + 1;
    Eval::print_time($total_time);
}
adjust_data_precision(\@data);
print_eval_output(\@data,\@names);
system("date");
exit;

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

__END__
