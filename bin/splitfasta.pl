#!/usr/bin/perl -w


=head1 NAME

=head1 USAGE

	*  'i|fa|fasta=s'	=> \$file2split,
  	   'depth=i'		=> \$depth,
	   'suffix=s'		=> \$suffix,
	   'dir=s'		=> \$splitfiledir

=cut

$|=1;		# autoflush
use strict;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";

my ($file2split,$splitfiledir);
my $depth=1;
my $suffix="";

GetOptions(
	'i|fa|fasta=s'  => \$file2split,
	'depth=i'	=> \$depth,
	'suffix=s'	=> \$suffix,
	'dir=s'		=> \$splitfiledir
);

if (!$file2split || !-s $file2split){pod2usage;}
if (!$splitfiledir){$splitfiledir="$file2split"."_dir".$depth;}
if ($suffix && $suffix!~/^\./){$suffix=".".$suffix;}
my $label=$file2split;
my ($flag);
my $filecount=0;
my $seqcount=0;
mkdir $splitfiledir;			 # to allow multiple runs of this program in wrapper
open (FILE, "$file2split");
print "Processing $file2split\n";
while (my $line=<FILE>){
	if ($line=~/^\s*$/){next;}	#empty line
	elsif ($line=~/^>(\S+)/)
	 {
		if ($depth==1){
			$label=$1;$label=~s/\|/_/g;
			$label=~s/_+$//;
		}
		$seqcount++;
		if (!$flag){
			$filecount++;
			my $outfile=$label;
			if ($depth>1){$outfile.="_".$filecount;}
			$outfile.=$suffix;
			open (OUT, ">$splitfiledir\/$outfile");
			$flag=1;
		}
		elsif ($seqcount>=$depth){
			$seqcount=0;
			$filecount++;
			my $outfile=$label;
			if ($depth>1){$outfile.="_".$filecount;}
			$outfile.=$suffix;
			close (OUT);
			open (OUT, ">$splitfiledir\/$outfile");
			if ($depth>1){
				print ".";
			}
		 }
		 print OUT $line;
	 }
	else { print OUT $line; }
}
close (FILE);
close (OUT);


print "\nProcessed $filecount files.\n";
