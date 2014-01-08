#!/usr/bin/perl -w

use strict;
use Getopt::Std;
use GTF;

use vars qw($opt_f $opt_l);
getopts('fl');
my $usage = "usage: $0 [-fl] <gtf file> <gene list>

Options:
     -f  Flag to fix the file.
     -l  Output those genes that are in the list.

\n";
die $usage unless (@ARGV == 2);
my ($filename,$genelist) = @ARGV;
my $fix_file = $opt_f;
my $fix_filename;
my $gtf = GTF::new({gtf_filename => $filename,
		    warning_fh   => \*STDERR});
my $genes = $gtf->genes;
my $last = 1;
open(LIST, "<$genelist") 
    or die "Could not open file: $genelist.\n";
my %list;
while(my $line = <LIST>){
    chomp $line;
    $list{$line} =1;
    $list{"$line.a"} =1;
    if($line =~ /^(.+)\.a$/){
	$list{$1} = 1;
    }
}
foreach my $gene (@$genes){
    if ($opt_l) {
	if (defined($list{$gene->gene_id})) {
	    $gene->output_gtf;
	}
    }
    else {
	unless(defined($list{$gene->gene_id})){
	    $gene->output_gtf;
	}
    }
}
exit;

__END__
