#!/usr/bin/perl
use Eval;
use GTF;
use strict;
use Getopt::Std;
use vars qw($opt_f $opt_g $opt_h $opt_q $opt_A);
getopts('fghqA');
my $usage = "$0 [fg] <filter file> <ann gtf> <pred gtf 1> [pred gtf 2] ...
Takes a filter file (see below) a annotation gtf and one or more 
prediction gtfs and filters them according to the filter file.
Options: 
  -f: List filter types
  -g: Inputs are gtf files instead of list files
  -A: Do not check for alternative splices. (Faster)
  -q: Quick load the gtf file.  Do not check them for errors.  
  -h: Display this help message
Filter File Format:
  A list of filter types with a single character label for each:
    A - Gene Correct
    B - Transcript All_Introns
    C - Exon Correct
  This list is followed by one or more empty lines then the filter string:
    (A&&B)||!C
";
if($opt_h){
    print $usage;
    exit(0);
}
my $Use_ase = 1;
if($opt_A){
    $Use_ase = 0;
}
if($opt_f){
    print_filter_types();
    exit(0);
}
die $usage unless (@ARGV >= 3);
my $Quick_Load = $opt_q;
my ($filter_file,$ann_file,@pred_files) = @ARGV;
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
open(FILT,$filter_file) or die "Couldn't open filter file $filter_file.\n";
my %filter_keys;
my $filter_string;
while(my $in = <FILT>){
    chomp $in;
    if($in eq ""){
	while(($in eq "") && ($in = <FILT>)){
	    chomp $in;
	}
	$filter_string = $in;
	last;
    }
    else{
	my @line = split /\s+/, $in;
	unless($line[1] eq '-'){
	    die "Bad filter file format\n";
	}
	check_filter($line[2],$line[3]);
	$filter_keys{$line[0]} = [$line[2],$line[3]];
    }
}
close(FILT);
my $filter = [];
parse_filter_string(\%filter_keys,$filter_string,$filter);

my $ann_gtf = load_func($ann_file);
my @pred_gtfs;
foreach my $pred_file (@pred_files){
    push @pred_gtfs, load_func($pred_file);
}

my @new_gtfs = Eval::filter_predictions($ann_gtf,\@pred_gtfs,$filter,1);

for(my $i = 0;$i <= $#new_gtfs;$i++){
    my $pred_file = $pred_files[$i];
    if($List_Mode){
	my $list_file = $pred_file;
	if($list_file =~ /^(.+)\.list$/){
	    $list_file = "$1.filtered.list";
	}
	else{
	    $list_file = "$list_file.filtered.list";
	}
	open(NEWLIST,">$list_file") or errzor_func("Couldn't open file $list_file.");
	for(my $j = 0;$j <= $#{$pred_gtfs[$i]};$j++){
	    my $new_file = $pred_gtfs[$i][$j]->filename;
	    if($new_file =~ /^(.+)\.gtf$/){
		$new_file = "$1.filtered.gtf";
	    }
	    else{
		$new_file = "$new_file.filtered.gtf";
	    }
	    open(NEW,">$new_file") or error_func("Couldn't open file $new_file.");
	    my $genes = $new_gtfs[$i][0]->genes;
	    $new_gtfs[$i][$j]->output_gtf_file(\*NEW);
	    close(NEW);
	    print NEWLIST "$new_file\n";
	}
	close(NEWLIST);
    }
    else{
	my $new_file = $pred_file;
	if($new_file =~ /^(.+)\.gtf$/){
	    $new_file = "$1.filtered.gtf";
	}
	else{
	    $new_file = "$new_file.filtered.gtf";
	}
	open(NEW,">$new_file") or error_func("Couldn't open file $new_file.");
	$new_gtfs[$i][0]->output_gtf_file(\*NEW);
	close(NEW);
    }
}
exit(0);


sub print_filter{
    my ($filter) = @_;
    if($$filter[0] eq 'Not'){
	print "!(";
	print_filter($$filter[1]);
	print ")";
    }
    elsif($$filter[0] eq 'Check'){
	print "($$filter[1] $$filter[2])";
    }
    else{
	print "(";
	print_filter($$filter[1]);
	print " $$filter[0] ";
	print_filter($$filter[2]);
	print ")";
    }
}

sub print_filter_types{
    my %filter_types = Eval::get_filter_types();
    foreach my $level (@{$filter_types{Levels}}){
	print "$level Filters:\n";
	foreach my $type (@{$filter_types{$level}}){
	    print "$level $type\n";
	}
	print "\n";
    }
}

sub check_filter{
    my ($level,$type) = @_;
    my %filter_types = Eval::get_filter_types();
    unless(defined($filter_types{$level})){
	error_func("Bad filter level $level\n");
    }
    my $ok = 0;
    foreach my $f_type (@{$filter_types{$level}}){
	if($type eq $f_type){
	    $ok = 1;
	    last;
	}
    }
    unless($ok){
	error_func("Bad filter type $type\n");
    }
}

sub parse_filter_string{
    my ($filter_keys,$filter_text,$filter) = @_;
    $filter_text =~ s/ //g;
    $filter_text =~ s/AND/&/gi;
    $filter_text =~ s/OR/|/gi;
    while($filter_text =~ /&&/){
	$filter_text =~ s/&&/&/g;
    }
    while($filter_text =~ /\|\|/){
	$filter_text =~ s/\|\|/\|/g;
    }
    return parse_filter_helper($filter_keys,$filter_text,$filter);
}

sub parse_filter_helper{
    my ($keys,$text,$filter) = @_;
    my @chars = split //, $text;
    my $trim = 0;
    if($#chars == -1){
	error_func("Bad filter string foramt: Unmatched connector.");
	return 0;
    }
    elsif($chars[0] eq '('){
	my $in = 1;
	my $index = 1;
	while(($in > 0) && ($index <= $#chars)){
	    if($chars[$index] eq '('){
		$in++;
	    }
	    elsif($chars[$index] eq ')'){
		$in--;
	    }
	    $index++;
	}
	$index--;
	if($in == 0){
	    if($index == $#chars){
		return(parse_filter_helper($keys,join("",@chars[1..($index-1)]),
					   $filter));
	    }
	    else{
		my $start = $index;
		if($chars[$index+1] eq '&'){
		    $$filter[0] = 'And';
		    $start+=2;		
		}
		elsif($chars[$index+1] eq '|'){
		    $$filter[0] = 'Or';
		    $start+=2;		
		}
		else{
		    $$filter[0] = 'And';
		    $start = $index+1;
		}
		$$filter[1] = [];
		$$filter[2] = [];
		return(parse_filter_helper($keys,join("",@chars[0..$index]),
					   $$filter[1]) &&
		       parse_filter_helper($keys,join("",@chars[$start..$#chars]),
					   $$filter[2]));
	    }
	}
	else{
	    error_func("Bad filter string foramt: mismatched parentheses.");
	    return 0;
	}
    }
    else{
	my $f1_stop = -1;
	my $f2_start = -1;
	my $open = 0;
	my $close = 0;
	for(my $index = 1;$index < $#chars;$index++){
	    if($chars[$index] eq '('){
		$open++;
	    }
	    elsif($chars[$index] eq ')'){
		$close++;
	    }
	    if(($open - $close) == 0){
		if($chars[$index] eq '&'){
		    $$filter[0] = 'And';
		    $f1_stop = $index-1;
		    $f2_start = $index+1;
		    last;
		}
		elsif($chars[$index] eq '|'){
		    $$filter[0] = 'Or';
		    $f1_stop = $index-1;
		    $f2_start = $index+1;
		    last;
		}
	    }
	}
	# this string has a connector (ie AB&&CD, !(A&B)||C)
	if($f1_stop != -1){
	    $$filter[1] = [];
	    $$filter[2] = [];
	    return(parse_filter_helper($keys,join("",@chars[0..$f1_stop]),$$filter[1]) &&
		   parse_filter_helper($keys,join("",@chars[$f2_start..$#chars]),
				       $$filter[2]));
	}
	# this is a series of and connected filters (ie ABC, !ABC, !(ABC), !(AB)C)
	if($chars[0] eq '!'){
	    if($#chars == 0){
		error_func("Bad filter string foramt: unmatched negation.");
		return 0;
	    }
	    # this is a grouped negation (ie !(AB&C), !(AB)C)
	    if($chars[1] eq '('){
		my $in = 1;
		my $index = 2;
		while(($in > 0) && ($index <= $#chars)){
		    if($chars[$index] eq '('){
			$in++;
		    }
		    elsif($chars[$index] eq ')'){
			$in--;
		    }
		    $index++;
		}
		$index--;
		if($in == 0){
		    # this is a whole group (ie !(ABC))
		    if($index == $#chars){
			$$filter[0] = 'Not';
			$$filter[1] = [];
			return(parse_filter_helper($keys,join("",@chars[1..$#chars]),
						   $$filter[1]));
		    }
		    # this is a group negation and something more (ie !(AB)C, !(AB)!C)
		    else{
			$$filter[0] = 'And';
			$$filter[1] = [];
			$$filter[2] = [];
			return(parse_filter_helper($keys,join("",@chars[0..$index]),
						   $$filter[1]) &&
			       parse_filter_helper($keys,
						   join("",@chars[$index+1..$#chars]),
						   $$filter[2]));
		    }
		}
	    }
	    # this is a single filter negation
	    else{
		# this is asingle negation and more (ie !AB !A!B)
		if($#chars > 1){
		    $$filter[0] = 'And';
		    $$filter[1] = [];
		    $$filter[2] = [];
		    return(parse_filter_helper($keys,$chars[0..1],$$filter[1]) &&
			   parse_filter_helper($keys,join("",@chars[2..$#chars]),
					       $$filter[2]));
		}
		# this is a single negation (ie !A)
		else{
		    $$filter[0] = 'Not';
		    $$filter[1] = [];
		    return(parse_filter_helper($keys,join("",$chars[1]),
					       $$filter[1]));
		}
	    }
	}
	else{
	    # This is a string of letters (ie ABC)
	    if($#chars >= 1){
		$$filter[0] = 'And';
		$$filter[1] = [];
		$$filter[2] = [];
		return(parse_filter_helper($keys,$chars[0],$$filter[1]) &&
		       parse_filter_helper($keys,join("",@chars[1..$#chars]),
					   $$filter[2]));
	    }
	    # this is a single filter (ie A)
	    elsif($#chars == 0){
		unless(defined($$keys{$chars[0]})){
		    if(($chars[0] eq '&') ||
		       ($chars[0] eq '|') ||
		       ($chars[0] eq '!')){
			error_func("Bad filter string format: misplaced character, ".
				   "$chars[0].");
			return 0;
		    }		
		    else{
			error_func("Bad filter string format: bad character, ".
				   "$chars[0].");
		    }
		    return 0;
		}
		$$filter[0] = "Check";
		$$filter[1] = $$keys{$chars[0]}[0];
		$$filter[2] = $$keys{$chars[0]}[1];
		return 1;
	    }
	}
    }
    error_func("Bad filter string format.");
    return 0;
}

sub error_func{
    my ($message) = @_;
    die "$message\n";
}

sub load_func{
    my ($filename) = @_;
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
			    no_check => $Quick_Load,
                mark_ase => $Use_ase);
		if((defined($info[1])) && 
		   ($info[1] ne "")){
		    #$hash{seq_filename} = $info[1];
		}
		if((defined($info[2])) && 
		   ($info[2] ne "")){
		    #$hash{conseq_filename} = $info[2];
		}
		push @$gtf, GTF::new(\%hash);
	    }
	}
	close(LIST);
	return $gtf;
    }
    else{
	return [GTF::new({gtf_filename => $filename,
                      no_check => $Quick_Load, 
                      mark_ase => $Use_ase})];
    }	
}

