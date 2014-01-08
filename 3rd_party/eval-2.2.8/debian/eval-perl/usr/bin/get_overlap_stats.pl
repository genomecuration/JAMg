#!/usr/bin/perl
use strict;
use lib '.';
use Eval;
use GTF;
use Getopt::Std;
use vars qw($opt_m $opt_g $opt_h $opt_v $opt_q);
my $Precision = "%.2f";
getopts('m:ghvq');
my $usage = "$0 [-ghv] [m mode] <list 1> <list 2> ...
Computes overlap statistics using the Eval package.  Inputs are gtf list files.
Options:
  -m <mode>: Specify overlap mode.  Must be a number selected from the list below.
      Default is mode 1.
  -g: Input files are in GTF format.
  -q: Quick load the gtf file.  Do not check them for errors.
  -v: Verbose mode.
  -h: Display this help message and exit.
Overlap Modes:
";
my @modes = Eval::get_overlap_mode_list();
for(my $i = 0;$i <= $#modes;$i++){
    $usage .= "  ".($i+1).") $modes[$i]\n";
}
if($opt_h){
    print $usage;
    exit(0);
}
die $usage unless (@ARGV >= 1);
my $mode = 1;
if($opt_m){
    $mode = $opt_m - 1;    
    unless(($mode >= 0) && ($mode <= $#modes)){
	die "Bad mode value, $opt_m.  Give -h option for help.\n";
    }
}
my $verbose = $opt_v;
my $Quick_Load = $opt_q;
my @gtfs;
my @names;
system("date");
if($verbose){
    print STDERR "@ARGV\n";
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
    if($opt_g){
	push @gtfs, [GTF::new({gtf_filename => $arg,
			       no_check => $Quick_Load})];
	push @names, $arg;
    }
    else{
	open(IN,$arg) or die "Could not open $arg for read.\n";
	my @this_list;
	while(my $in = <IN>){
	    chomp $in;
	    if($in =~ /\S/){
		push @this_list, GTF::new({gtf_filename => $in,
					   no_check => $Quick_Load});
	    }
	}
	close(IN);
	push @gtfs, \@this_list;
	push @names, $arg;
    }
}
if($verbose){
    print STDERR "done\n";
    print STDERR "running...";
}
my $start_time = time;
my %data = Eval::get_overlap_statistics(\@gtfs,$modes[$mode],$verbose);
print get_overlap_stats_text(\%data,\@names,$modes[$mode]);
system("date");
if($verbose){
    my $end_time = time;
    my $total_time = $end_time - $start_time + 1;
    Eval::print_time($total_time);
}
exit;

sub get_overlap_stats_text{
    my ($data,$pred_names,$overlap_type) = @_;
    $overlap_type =~ s/_/ /g;
    my $report = "Overlap Type: $overlap_type\n";
    my %total = (all => 0);
    my @labels = Eval::get_overlap_labels($#$pred_names);
    foreach my $label (@labels){
	$total{$label} = 0;
    }    
    foreach my $group (keys %$data){
	$total{all} += $$data{$group}{total};
	foreach my $label (@labels){
	    $total{$label} += $$data{$group}{$label};
	}
    }
    for(my $i = 0;$i <= $#labels;$i++){
	$report .= "$labels[$i] -> $$pred_names[$i] ($total{$labels[$i]})\n";
    }
    $report .= "\nname\tcount";
    for(my $i = 0;$i <= $#labels;$i++){
	$report .= "\t$labels[$i]%";
    }
    $report .= "\ttotal\n";
    foreach my $group (sort {length($a) <=> length($b) ||
				 $a cmp $b} (keys %$data)){
	$report .= $group ."\t".$$data{$group}{total};
	for(my $i = 0;$i <= $#labels;$i++){
	    $report .= sprintf "\t".$Precision, (100*$$data{$group}{$labels[$i]}/
						 $total{$labels[$i]});
	    $report .= "%";
	}
	$report .= sprintf "\t".$Precision, (100*$$data{$group}{total}/$total{all});
	$report .= "%\n";
    }
    return $report;
}
