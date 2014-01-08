#!/usr/bin/perl
use Eval;
use GTF;
use strict;
use Getopt::Std;
use vars qw($opt_g $opt_h $opt_m $opt_q);
getopts('ghm:q');
my $usage = "$0 [-g] [-m mode] <max val> <bin size> <pred gtf 1> [pred gtf 2] ...
Takes the maximum value to report in the distribution, the size of bins to 
report data in, and one of more gtf sets and creates outputs the distribution 
to standard out.
Options: 
  -m <mode>: Specify distribution mode.  Must be a number selected from the 
      list below.  Default is mode 1.
  -g: Inputs are gtf files instead of list files
  -q: Quick load the gtf file.  Do not check them for errors.
  -h: Display this help message
Distribution Modes:
";
my @modes = Eval::get_distribution_type_list();
for(my $i = 0;$i <= $#modes;$i++){
    $usage .= "  ".($i+1).") $modes[$i]\n";
}
if($opt_h){
    print $usage;
    exit(0);
}
die $usage unless (@ARGV >= 3);
my $Quick_Load = $opt_q;
my ($max,$size,@pred_files) = @ARGV;
my $gtf_mode = 1;
foreach my $arg (@pred_files){
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
my $mode = 1;
if($opt_m){
    $mode = $opt_m - 1;    
    unless(($mode >= 0) && ($mode <= $#modes)){
	die "Bad mode value, $opt_m.  Give -h option for help.\n";
    }
}
my $dist_type = $modes[$mode];
my %dist = Eval::get_distribution_type_hash();
unless(defined($dist{$dist_type})){
    error_func("Bad distribution type.");
}
$dist{$dist_type} = 1;
my @pred_gtfs;
foreach my $pred_file (@pred_files){
    push @pred_gtfs, load_func($pred_file);
}
my @data = Eval::get_distribution(\@pred_gtfs,\%dist,1);
my @bins;
print "$dist_type Distribution\n";
for(my $i = 0;$i <= $#pred_gtfs;$i++){
    print "\t$pred_files[$i]";
}
print "\n";
for(my $i = $size;$i < ($max+$size);$i+=$size){
    if($i > $max){
	$i = $max;
    }
    print "".($i-$size)."-$i";
    $bins[$i] = [];
    for(my $pred = 0;$pred <= $#pred_gtfs;$pred++){
	my $count = 0;
	for(my $j = ($i - $size);$j <= $i;$j++){
	    if(defined($data[$pred]{$dist_type}{$j})){
		$count += $data[$pred]{$dist_type}{$j};
		delete($data[$pred]{$dist_type}{$j});
	    }
	}
	print "\t$count";
    }
    print "\n";
}
print ">$max";
for(my $pred = 0;$pred <= $#pred_gtfs;$pred++){
    my $count = 0;
    foreach my $val (keys %{$data[$pred]{$dist_type}}){
	$count += $data[$pred]{$dist_type}{$val};
    }
    print "\t$count";
}
print "\n";
exit(0);
				  
sub error_func{
    my ($message) = @_;
    die "$message\n";
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
	return [GTF::new({gtf_filename => $filename})];
    }	
}

