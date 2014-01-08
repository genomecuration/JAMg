#!/usr/bin/perl
use strict;
use lib '.';
use Eval;
use GTF;
use Getopt::Std;
use vars qw($opt_b $opt_B $opt_g $opt_h $opt_m $opt_v $opt_q $opt_x);
getopts('b:B:ghm:vqx:');
my $MIN_DISP_LEN = 25;
my $MIN_VAL_LEN = 10;

my $Precision = "%.3f";

my $usage = "$0 <annotation list> <prediction list 1> [prediction list 2] ...
Create a graph of intron performance vs intron length 
Options:
  -m <min_bin_start>: Sets the minimum bin start [default: min intron length];
  -x <max_bin_stop>: Sets the maximum bin end [default: max intron length];
  -b <bin_size>: Sets the bin size [default: 1/10 length range] 
                 Cannot be used with -B
  -B <bin_count>: Sets the number of bins [default: 10] 
                  Cannot be used with -b
  -g: Input files are gtf not lists
  -q: Quick load the gtf file.  Do not check them for errors.
  -v: Verbose mode
  -h: Display this help message and exit
";
if($opt_h){
    print $usage;
    exit(0);
}
die $usage unless (@ARGV >= 2);

my $verbose = $opt_v;
my $Quick_Load = $opt_q;

my @gtfs;
my @names;
if($verbose){
    print STDERR "Loading...";
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
my $Min_Bin_Start = -1;
my $Max_Bin_Stop = -1;
my $Bin_Size = -1;
my $Bin_Count = 10;
if($opt_m){
    $Min_Bin_Start = $opt_m;
    if($Min_Bin_Start < 0){
	die "Min bin start (-m) must be greater than 0\n";
    } 
}
if($opt_x){
    $Max_Bin_Stop = $opt_x;
    if($Max_Bin_Stop <= $Min_Bin_Start){
	die "Max bin stop (-x) must be greater than min bin start (-m)\n";
    } 
}
if($opt_b){
    if($opt_B){
	die "Cannot give -b and -B options\n";
    }
    $Bin_Size = $opt_b;
    $Bin_Count = -1;
    if($Bin_Size <= 0){
	die "Bin size (-b) must be greater than 0\n";
    }
}
if($opt_B){
    $Bin_Count = $opt_B;
    $Bin_Size = -1;
    if($Bin_Count <= 0){
	die "Bin count (-B) must be greater than 0\n";
    }
}

my $start_time = time;
my @names;
foreach my $arg (@ARGV){
    my $name = $arg;
    while($name =~ /\/(\S+)$/){
	$name = $1;
    }
    my @this_list;
    if($opt_g){
	push @this_list, GTF::new({gtf_filename => $arg,
				   no_check => $Quick_Load});
    }
    else{
	open(IN,$arg) or die "Couldn't open $arg";
	while(my $in = <IN>){
	    chomp $in;
	    my @data = split /\t/, $in;
	    $in = $data[0];
	    if($in =~ /\S/){
		push @this_list, GTF::new({gtf_filename => $in,
					   no_check => $Quick_Load});
	    }
	}
    }
    push @gtfs, \@this_list;
    push @names, $name;
}
my $end_time = time;
my $total_time = $end_time - $start_time + 1;
$start_time = time;
if($verbose){
    print STDERR "done\n    ";
    Eval::print_time($total_time);
    print STDERR "Running...";
}
my $ann = $gtfs[0];
my $ann_name = $names[0];
my @preds;
my @pred_names;
for(my $i = 1;$i <= $#gtfs;$i++){
    push @preds, $gtfs[$i];
    push @pred_names, $names[$i];
}
my @data = make_graphs($ann,\@preds);
if($verbose){
    $end_time = time;
    $total_time = $end_time - $start_time + 1;
    print STDERR "done\n    ";
    Eval::print_time($total_time);
}
print_graph(\@data,$ann_name,\@pred_names);
exit;

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

sub make_graphs{
    my ($ann,$preds) = @_;
    my @a_introns;
    my @p_introns;
    #get all ann introns
    for(my $f = 0; $f <= $#$ann; $f++){
	$a_introns[$f] = [];
	foreach my $a_tx (@{$$ann[$f]->transcripts}){
	    push @{$a_introns[$f]}, @{$a_tx->introns};
	}
	$a_introns[$f] = [sort {$a->start <=> $b->start || $a->stop <=> $b->stop} 
			  @{$a_introns[$f]}];
    }
    my $max_length;
    my $min_length;
    my $f = 0;
    until(defined($max_length)){
	foreach my $a_tx (@{$$ann[$f]->transcripts}){
	    my $ais = $a_tx->introns;
	    if($#$ais > 0){
		$max_length = $$ais[0]->length;
		$min_length = $$ais[0]->length;
		last;
	    }
	}
    }
    #get all pred introns
    for(my $f = 0; $f <= $#$ann; $f++){
	@p_introns[$f] = [];
	for(my $i = 0; $i <= $#preds;$i++){
	    foreach my $p_tx (@{$$preds[$i][$f]->transcripts}){
		#print "add ".$p_tx->id."\n";
		foreach my $int (@{$p_tx->introns}){
		    $int->set_tag(get_tag());
		    push @{$p_introns[$f]}, {intron => $int,
					     set => $i};
		    if($int->length < $min_length){
			$min_length = $int->length;
		    }
		    if($int->length > $max_length){
			$max_length = $int->length;
		    }
		}
	    }
	}
	$p_introns[$f] = [sort {$$a{intron}->start <=> $$b{intron}->start || 
				    $$a{intron}->stop <=> $$b{intron}->stop} 
			  @{$p_introns[$f]}];
    }
    #do the comparisons
    for(my $f = 0; $f <= $#$ann; $f++){
	my $first_p = 0;
	for(my $a = 0; $a <= $#{$a_introns[$f]};$a++){
	    #load this ann introns
	    my $ai = $a_introns[$f][$a];
	    if($ai->length < $min_length){
		$min_length = $ai->length;
	    }
	    if($ai->length > $max_length){
		$max_length = $ai->length;
	    }
	    $ai->set_tag(get_ann_tag($#preds));
	    my $a_start = $ai->start;
	    my $a_stop = $ai->stop;
	    my $a_strand = $ai->strand;
	    if($first_p > $#{$p_introns[$f]}){
		next;
	    }
	    #check min/max ann length
	    my $a_length = $ai->length;
	    if($a_length > $max_length){
		$max_length = $a_length;
	    }
	    if($a_length < $min_length){
		$min_length = $a_length;
	    }
	    #move to next overlapping pred
	    while(($first_p <= $#{$p_introns[$f]}) && 
		  ($a_start > $p_introns[$f][$first_p]{intron}->stop)){
		$first_p++;
	    }
	    #compare overlapping preds
	    my $p = $first_p;
	    while(($p <= $#{$p_introns[$f]}) && 
		  ($a_stop >= $p_introns[$f][$p]{intron}->start)){
		if($a_strand eq $p_introns[$f][$p]{intron}->strand){
		    compare_overlapping_introns($ai,$p_introns[$f][$p]{intron},
						$p_introns[$f][$p]{set});
		}
		$p++;
	    }
	}
    }
    #resort by length
    for(my $f = 0; $f <= $#$ann; $f++){
	$a_introns[$f] = [sort {$a->length <=> $b->length} @{$a_introns[$f]}];
	$p_introns[$f] = [sort {$$a{intron}->length <=> $$b{intron}->length} 
			  @{$p_introns[$f]}];
    }
    #get bin info
    my $graph_min = $Min_Bin_Start;
    if($graph_min < 0){
	$graph_min = $min_length;
    }
    my $graph_max = $Max_Bin_Stop;
    if($graph_max < $graph_min){
	$graph_max = $max_length;
    }
    my $bin_size = $Bin_Size;
    if($bin_size < 0){
	$bin_size = ($graph_max-$graph_min)/$Bin_Count;
    }
    #prepare graph data structure
    my @data;
    my $bin_start = $graph_min; 
    my $bin_num = 1;
    my $bin_stop = int($graph_min + $bin_num*$bin_size);
    while($bin_stop <= $graph_max){
	my $info = {start => $bin_start,
		    stop => $bin_stop,
		    data => {ann => [],
			     pred => []}};
	for(my $i = 0; $i <= $#preds;$i++){
	    $$info{data}{ann}[$i] = get_tag();
	    $$info{data}{ann}[$i]{Count} = 0;
	    $$info{data}{pred}[$i] = get_tag();
	    $$info{data}{pred}[$i]{Count} = 0;
	}
	push @data,$info;
	$bin_start = $bin_stop+1;
	$bin_num++;
	$bin_stop = int($graph_min + $bin_num*$bin_size);
    }
    $data[$#data]{stop} = $graph_max;
    my @stats = get_ordered_stats();
    for(my $f = 0; $f <= $#$ann; $f++){
	#collect ann info 
	my $a = 0;
	my $d = 0;
	while(($a <= $#{$a_introns[$f]}) &&
	      ($a_introns[$f][$a]->length < $graph_min)){
	    $a++;
	}
	while(($a <= $#{$a_introns[$f]}) &&
	      ($a_introns[$f][$a]->length <= $graph_max)){
	    while($a_introns[$f][$a]->length > $data[$d]{stop}){
		$d++;
	    }
	    my $at = $a_introns[$f][$a]->tag;
	    for(my $i = 0; $i <= $#preds;$i++){
		$data[$d]{data}{ann}[$i]{Count}++;
		foreach my $s (@stats){
		    $data[$d]{data}{ann}[$i]{$s} += $$at[$i]{$s};
		}
	    }
	    $a++;
	}    
	#colect pred data
	my $p = 0;
	$d = 0;
	while(($p <= $#{$p_introns[$f]}) &&
	      ($p_introns[$f][$p]{intron}->length < $graph_min)){
	    $p++;
	}
	while(($p <= $#{$p_introns[$f]}) &&
	      ($p_introns[$f][$p]{intron}->length <= $graph_max)){
	    while($p_introns[$f][$p]{intron}->length > $data[$d]{stop}){
		$d++;
	    }
	    $data[$d]{data}{pred}[$p_introns[$f][$p]{set}]{Count}++;
	    my $pt = $p_introns[$f][$p]{intron}->tag;
	    foreach my $s (@stats){
		$data[$d]{data}{pred}[$p_introns[$f][$p]{set}]{$s} += $$pt{$s};
	    }
	    $p++;
	}
    }
    return @data;
}

sub print_graph{
    my ($graph,$ann_name,$pred_names) = @_;
    my @stats = get_ordered_stats();
    print "Annotation:\t$ann_name\n";
    print "Prediction:";
    foreach my $name (@$pred_names){
	print "\t$name";
    }
    print "\n";
    print "Graph from $$graph[0]{start} to $$graph[$#$graph]{stop}\n\n";
    print "Max Length";
    foreach my $bin (@$graph){
	print "\t$$bin{stop}";
    }
    print "\n\n";
    print "$ann_name Count";
    foreach my $bin (@$graph){
	my $total = $$bin{data}{ann}[0]{Count};
	print "\t$total";
    }
    print "\n";
    for(my $i = 0; $i <= $#$pred_names;$i++){
	print "$$pred_names[$i] Count";
	foreach my $bin (@$graph){
	    my $total = $$bin{data}{pred}[$i]{Count};
	    print "\t$total";
	}
	print "\n";
    }
    print "\n";
    foreach my $stat (@stats){
	for(my $i = 0; $i <= $#$pred_names;$i++){
	    print "$$pred_names[$i] $stat";
	    foreach my $bin (@$graph){
		my $count = $$bin{data}{pred}[$i]{$stat};
		print "\t$count";
	    }
	    print "\n";
	    print "$$pred_names[$i] Ann $stat";
	    foreach my $bin (@$graph){
		my $count = $$bin{data}{ann}[$i]{$stat};
		print "\t$count";
	    }
	    print "\n";
	    print "$$pred_names[$i] $stat Sp";
	    foreach my $bin (@$graph){
		my $count = $$bin{data}{pred}[$i]{$stat};
		my $total = $$bin{data}{pred}[$i]{Count};
		if($total == 0){
		    printf("\t$Precision",0);
		}
		else{
		    printf("\t$Precision", $count/$total);
		}
	    }
	    print "\n";
	    print "$$pred_names[$i] $stat Sn";
	    foreach my $bin (@$graph){
		my $count = $$bin{data}{ann}[$i]{$stat};
		my $total = $$bin{data}{ann}[$i]{Count};
		if($total == 0){
		    printf("\t$Precision",0);
		}
		else{
		    printf("\t$Precision", $count/$total);
		}
	    }
	    print "\n";
	}
	print "\n";
    }
}

sub compare_overlapping_introns{
    my ($a,$p,$set) = @_;
    my $at = $a->tag;
    my $pt = $p->tag;
    $$at[$set]{Overlap} = 1;
    $$pt{Overlap} = 1;
    my $lm = 0;
    my $hm = 0;
    if($a->start == $p->start){
	$lm = 1;
    }
    if($a->stop == $p->stop){
	$hm = 1;
    }
    if($lm == 1 && $hm == 1){
	$$at[$set]{Correct} = 1;
	$$at[$set]{Overlap80p} = 1;
	$$at[$set]{Match5} = 1;
	$$at[$set]{Match3} = 1;
	$$pt{Correct} = 1;
	$$pt{Overlap80p} = 1;
	$$pt{Match5} = 1;
	$$pt{Match3} = 1;
    }
    else{
	if($a->strand eq '+'){
	    if($lm){
		$$at[$set]{Match5} = 1;
		$$pt{Match5} = 1;
	    }
	    if($hm){
		$$at[$set]{Match3} = 1;
		$$pt{Match3} = 1;
	    }
	}
	else{
	    if($hm){
		$$at[$set]{Match5} = 1;
		$$pt{Match5} = 1;
	    }
	    if($lm){
		$$at[$set]{Match3} = 1;
		$$pt{Match3} = 1;
	    }
	}
	my $overlap = min($a->stop,$p->stop)-max($a->start,$p->start)+1;
	if(($overlap/$a->length >= .8) && ($overlap/$p->length > .8)){
	    $$at[$set]{Overlap80p} = 1;
	    $$pt{Overlap80p} = 1;
	}
    }
}

sub max{
    my ($a,$b) = @_;
    if($a > $b){
	return $a;
    }
    else{
	return $b;
    }
}

sub min{
    my ($a,$b) = @_;
    if($a < $b){
	return $a;
    }
    else{
	return $b;
    }
}


sub get_ordered_stats{
    return qw(Correct Overlap Overlap80p Match5 Match3);
}

sub get_tag{
    my %tag;
    foreach my $stat (get_ordered_stats()){
	$tag{$stat} = 0;
    }
    return \%tag;
}

sub get_ann_tag{
    my ($count) = @_;
    my @tags;
    for(my $i = 0; $i < $count; $i++){
	push @tags, get_tag();
    }
    return \@tags;
}

__END__
