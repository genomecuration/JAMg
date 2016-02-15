#!/usr/bin/env perl

use strict;
use warnings;
our $DEBUG = 0;

my $usage = "usage: $0 extFile\n\n";

my $extFile = $ARGV[0] or die $usage;
die unless $extFile && -s $extFile;
my $outfile = $extFile.'Col';

my $cpus = 4;
my $sort_buffer = '5G';  # will run up to two sorts in parallel
my $tmpdir = $ENV{'TMP'};
$tmpdir = $ENV{'TMPDIR'} if !$tmpdir;
$tmpdir = '/tmp' if !$tmpdir;

my $sort_exec = &check_sort_version;
# sort based on hit ID and strand and start
warn "Sorting...\n" if $DEBUG;
open (TMP,"$extFile");
my $header = <TMP>;
close TMP;

system("$sort_exec -nk1 $extFile | $sort_exec -s -nk6,6 | $sort_exec -s -k9,9 > $extFile.sorted") unless -s "$extFile.sorted";
open (EXT, "$extFile.sorted") or die "Cannot open $extFile.sorted\n";
open (OUT1,">$extFile.collapsed1") || die $!;
my $discard = <EXT>;
my $previous_ln;

warn "Processing sorted fiel...\n" if $DEBUG;
while (my $ln = <EXT>) {
    next if $ln=~/^\s*$/;
    chomp($ln);
    $ln =~ s/^\s+//; #rm leading whitespace
    ## using var names as in ext.c
    my ($next_dstart, $next_dend, $next_score, $next_astart, $next_aend, $next_orient, $next_zero1, $next_zero2, $next_acc) = split (/\s+/,$ln);
    next if (!$next_dstart || $next_dstart!~/^\d+$/ || $next_dstart < 1 || !$next_dend || $next_dend!~/^\d+$/ || $next_dend < 1 || !$next_score || $next_score!~/^\d+$/ || $next_score < 1 || !$next_acc);
    
    if (!$previous_ln){
	$previous_ln = $ln if $ln=~/\d+/; 
	next;
    }

    my ($prev_dstart, $prev_dend, $prev_score, $prev_astart, $prev_aend, $prev_orient, $prev_zero1, $prev_zero2, $prev_acc) = split (/\s+/,$previous_ln);
	
    if ($next_acc && $prev_acc && $next_acc eq $prev_acc 
	&& $next_orient == $prev_orient && $next_dstart <= $prev_dend) {
	
	## merge overlapping entry:
	my @dcoords = sort {$a<=>$b} ($prev_dstart, $prev_dend, $next_dstart, $next_dend);
	my $dstart = shift @dcoords;
	my $dend = pop @dcoords;
	my @acoords = sort {$a<=>$b} ($prev_astart, $prev_aend, $next_astart, $next_aend);
	my $astart = shift @acoords;
	my $aend = pop @acoords;
	my @scores = sort {$a<=>$b} ($prev_score, $next_score);
	my $score = pop @scores;
	
	$previous_ln = "$dstart $dend $score $astart $aend $prev_orient $prev_zero1 $prev_zero2 $prev_acc";	
	warn "expanding current chain.\n" if $DEBUG;
    } 
    else {
	$previous_ln = $ln;
	next if ($prev_dstart > 9999999999 || $prev_dend > 9999999999);
	printf OUT1 ("%10d %10d %6d %7d %5d %1d %5d %5d %s\n",
          $prev_dstart, $prev_dend, $prev_score, $prev_astart, $prev_aend, $prev_orient, $prev_zero1, $prev_zero2, $prev_acc) if ($prev_acc && $prev_score && $prev_score >0);
    }
}
#last line
my ($prev_dstart, $prev_dend, $prev_score, $prev_astart, $prev_aend, $prev_orient, $prev_zero1, $prev_zero2, $prev_acc) = split (/\s+/,$previous_ln) if $previous_ln;
next if ($prev_dstart > 9999999999 || $prev_dend > 9999999999);
printf OUT1 ("%10d %10d %6d %7d %5d %1d %5d %5d %s\n",$prev_dstart, $prev_dend, $prev_score, $prev_astart, $prev_aend, $prev_orient, $prev_zero1, $prev_zero2, $prev_acc) if ($prev_acc && $prev_score && $prev_score >0);

close EXT;
close OUT1;

warn "Resorting...\n" if $DEBUG;
open (OUT,">$outfile");
print OUT $header;
close OUT;
# i added uniq here because there is a bug somewhere i cant find
system("$sort_exec -nk1,1 -nk2,2 $extFile.collapsed1|uniq >> $outfile");
system("sed -i '\$d' $outfile");
unlink("$extFile.collapsed1");
unlink("$extFile.sorted");


#########################

sub check_sort_version(){
	my ($sort_exec) = &check_program('sort');
	my @v=`$sort_exec --version`;

	if ($v[0] && $v[0]=~/(\d+)\.(\d+)\s*$/){
		my $major = $1;
		my $minor = $2;
		if ($major >= 8 && $minor >= 6){
			return "$sort_exec -T $tmpdir --parallel $cpus -S $sort_buffer";
		}else{
			return "$sort_exec -S $sort_buffer -T $tmpdir";
		}
	}else{
		die "Sort of coreutils not found!";
	}
}

sub check_program() {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  die "Error, path to required $prog cannot be found\n"
    unless $path =~ /^\//;
  chomp($path);
  $path = readlink($path) if -l $path;
  push( @paths, $path );
 }
 return @paths;
}
