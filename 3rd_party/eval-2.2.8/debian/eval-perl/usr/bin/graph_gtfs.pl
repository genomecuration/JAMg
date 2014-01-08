#!/usr/bin/perl
use Eval;
use GTF;
use strict;
use Getopt::Std;
use vars qw($opt_G $opt_g $opt_h $opt_r $opt_q);
getopts('Gghr:q');
my $usage = "$0 [-gGh] [-r <file>] <graph file> <ann> <pred 1> [pred 2] ...
Takes a graph file (see below), and annotation and one or more predictions and 
creates each graph specified by the graph file for each pred.
  Options:
    -G: Display list of possible x and y values for graphs
    -g: Load GTFs instead of lists of GTFs
    -q: Quick load the gtf file.  Do not check them for errors.
    -r <resolution file>:  Load resolution from this file
        instead of users .eval.rc or default 
    -h: Display this help message
  Graph file format:
    Each line should be of the format:
    \"y_level::y_type::y_stat vs x_type::x_level\"
    where options for y_level,y_type_,y_stat, and x_type can be found by 
    giving the -G option
  Resolution file format:
    Each line should be in one of the following formats (all fields are  
    separated by tabs):
    1)\"x_type User # # #\"
      where values for \'x_type\' can be found by using the -G option and bins of 
      values for \'x_type\' are made from each \'#\' to the next \'#\'
    2)\"x_type Uniform min size count\"
      where \'min\' is the minimum value of any bin, bins are of size \'size\',
      and there are a total of \'count\' bins, and \'x_type\'is as above
";
if($opt_h){
    print $usage;
    exit(0);
}
if($opt_G){
    print_graph_types();
    exit(0);
}
my $Quick_Load = $opt_q;
my $gtf_mode = 1;
for(my $i = 1;$i <= $#ARGV;$i++){
    my $arg = $ARGV[$i];
    unless($arg =~ /\.g[tf]f$/){
	$gtf_mode = 0;
    }
}
if($gtf_mode){
    $opt_g = 1;
}
my $List_Mode = 1;
if($opt_g){
    $List_Mode = 0;
}
die $usage unless @ARGV >= 3;
my($graph_file,$ann_file,@pred_files) = @ARGV;
# get graph resolution
my %resolution = get_default_resolution();
my $rcfile = "~`whoami`/.evalrc";
if($opt_r){
    open(RES,$opt_r) or die "Couldn't open resolution file, $opt_r.\n";
    my @text = <RES>;
    close(RES);
    parse_resolution_text(\%resolution,@text);
}
elsif(-e $rcfile){
    my @text;
    open(RES,$rcfile) or die "Couldn't open rc file, $rcfile.\n";
    my $in_graph = 0;
    while(my $in = <RES>){
	if($in =~ /^\@Graph/){
	    $in_graph = 1;
	}
	elsif($in =~ /^\@Done/){
	    $in_graph = 0;
	}
	elsif($in_graph){
	    push @text,$in;
	}
    }
    close(RES);
    parse_resolution_text(\%resolution,@text);
}

#get graphs
my @graphs;
my %x_splits;
my %data_types = Eval::get_stats_struct();
foreach my $split (Eval::get_graph_x_types()){
    $x_splits{$split} = {};
    foreach my $level (Eval::get_graph_x_levels()){
	$x_splits{$split}{$level} = 0;
    }
}
open(IN,$graph_file) or die "Couldn't open graph file, $graph_file.\n";
while(my $in = <IN>){
    chomp $in;
    if(($in eq '') ||
       ($in =~ /^\s*$/)){
	next;
    }
    if($in =~ /^\s*(\S+)\s+vs\s+(\S+)\s*$/){
	my $y_val = $1;
	my $x_val = $2;
	my @x_vals = split /::/, $x_val;
	if(defined($x_splits{$x_vals[1]})){
	    if(defined($x_splits{$x_vals[1]}{$x_vals[0]})){
		my @y_vals = split /::/, $y_val;
		if(defined($data_types{$y_vals[0]})){
		    if(defined($data_types{$y_vals[0]}{$y_vals[1]})){
			if(defined($data_types{$y_vals[0]}{$y_vals[1]}{$y_vals[2]})){
			    push @graphs, {y => {level => $y_vals[0],
						 type => $y_vals[1],
						 stat => $y_vals[2]},
					   x => {level => $x_vals[0],
						 split => $x_vals[1]}};
			    $x_splits{$x_vals[1]}{$x_vals[0]} = 1;
			}
			else{
			    die "Bad y stat, $y_vals[2], in graph file.\n";
			}
		    }
		    else{
			die "Bad y type, $y_vals[1], in graph file.\n";
		    }
		}
		else{
		    die "Bad y level, $y_vals[0], in graph file.\n";
		}
	    }
	    else{
		die "Bad x level, $x_vals[1], in graph file.\n";
	    }
	}
	else{
	    die "Bad x axis, $x_vals[0], in graph file.\n";
	}
    }
    else{
	die "Bad graph file format:\n\t\"$in\"\n";
    }
}
close(IN);
my @graph_types;
foreach my $split (Eval::get_graph_x_types()){
    foreach my $level (Eval::get_graph_x_levels()){
	if($x_splits{$split}{$level}){
	    push @graph_types, {split => $split,
			       level => $level};
	}
    }
}
#load the gtf sets
my $ann_gtf = load_func($ann_file);
my @pred_gtfs;
foreach my $pred_file (@pred_files){
    push @pred_gtfs, load_func($pred_file);
}
#make the graphs
my @data = Eval::make_graphs($ann_gtf,\@pred_gtfs,\@graph_types,\%resolution,1);
#report the results
foreach my $graph (@graphs){
    my $x_split = $$graph{x}{split};
    my $x_level = $$graph{x}{level};
    my $level = $$graph{y}{level};
    my $type = $$graph{y}{type};
    my $stat = $$graph{y}{stat};
    print "$level $type $stat vs $x_level $x_split\n";

    for(my $bin = 0;$bin <= $#{$data[0]{$x_split}{$x_level}};$bin++){
	print "".$data[0]{$x_split}{$x_level}[$bin]{min}." - ".
	    $data[0]{$x_split}{$x_level}[$bin]{max};
	for(my $i = 0;$i <= $#pred_gtfs;$i++){
	    print "\t".$data[$i]{$x_split}{$x_level}[$bin]{data}{$level}{$type}{$stat};
	}
	print "\n";
    }
}
exit(0);

sub parse_resolution_text{
    my ($res,@text) = @_;
    foreach my $line (@text){
	chomp $line;
	if($line =~ /^\s+(.+)$/){
	    $line = $1;
	}
	my ($graph_type,$split_type,@info) = split /\t/,$line;
	my $ok = 0;
	foreach my $type (Eval::get_graph_x_types()){
	    if($type eq $graph_type){
		$ok = 1;
		last;
	    }
	}
	unless($ok){
	    die "Bad graph type, $graph_type.  Give -G flag to see valid options.\n";
	}
	if($split_type eq "User"){
	    @info = sort {$a <=> $b} @info;
	    my @bins;
	    for(my $i = 0;$i < $#info;$i++){
		$$res{$graph_type}{user}[$i] = {min => $info[$i],
						max => $info[$i+1]};
	    }
	}
	elsif($split_type eq "Uniform"){
	    for(my $i = 0;$i <= 2;$i++){
		unless(defined($info[$i])){
		    $info[$i] = "";
		}
	    }
	    $$res{$graph_type}{uniform} = {};
	    if($info[0] ne ''){
		$$res{$graph_type}{uniform}{min} = $info[0];
	    }
	    if($info[1] ne ''){
		$$res{$graph_type}{uniform}{size} = $info[1];
	    }
	    if($info[2] ne ''){
		$$res{$graph_type}{uniform}{count} = $info[2];
	    }
	}
	else{
	    die "Bad graph resolution type, $split_type.  Should be \"User\"".
		"or \"Uniform\".\n";
	}
    }
}

sub get_default_resolution{
    my %res;
    foreach my $x_split (Eval::get_graph_x_types()){
	$res{$x_split}{uniform} = {min => 0,
				   max => 1,
				   size => .1,
				   count => 10};
    }
    return %res;
}

sub print_graph_types{
    print "X-axis:\n";
    print "Splits:\n";
    foreach my $x_type (Eval::get_graph_x_types()){
	print "\t$x_type\n";
    }
    print "Levels:\n";
    foreach my $x_level (Eval::get_graph_x_levels()){
	print "\t$x_level\n";
    }
    print "\n\n\n";
    my %y_types = Eval::get_graph_y_types();
    print "Y-axis:\n";
    foreach my $level (@{$y_types{Levels}}){
	print "\t$level level\n";
	print "\t\tTypes:\n";
	foreach my $type (@{$y_types{$level}{Type}}){
	    print "\t\t\t$type\n";
	}
	print "\t\tStats:\n";
	foreach my $stat (@{$y_types{$level}{Stat}}){
	    print "\t\t\t$stat\n";
	}
	print "\n";
    }
}

sub load_func{
    my ($filename) = @_;
    print STDERR "Loading $filename\n";
    unless(-e $filename){
	error_func("File $filename does not exist\n");
    }
    if(($filename =~ /^.*\.list/) || 
       ($List_Mode && ($filename !~/^.*\.g[tf]f$/))){
	open(LIST,$filename) 
	    or error_func("Unable to open $filename for input.",1);
	my $gtf = [];
	while(my $file = <LIST>){
	    if($file =~ /^\s+(\S+)\s+$/){
		$file = $1;
	    }
	    elsif($file =~ /^\s+(\S+)$/){
		$file = $1;
	    }
	    elsif($file =~ /^(\S+)\s+$/){
		$file = $1;
	    }
	    if(($file =~ /\S/) &&
	       ($file !~ /^\#/)){
		chomp $file;
		my @info = split /\t/, $file;	
		my %hash = (gtf_filename => $info[0],
			    no_check => $Quick_Load);
		if((defined($info[1])) && 
		   ($info[1] ne "")){
		    $hash{seq_filename} = $info[1];
		}
		if((defined($info[2])) && 
		   ($info[2] ne "")){
		    $hash{conseq_filename} = $info[2];
		}
		push @$gtf, GTF::new(\%hash);
	    }
	}
	close(LIST);
	return $gtf;
    }
    else{
	return [GTF::new({gtf_filename => $filename})];
    }	
}

sub error_func{
    my ($message) = 0;
    die $message."\n";
}
