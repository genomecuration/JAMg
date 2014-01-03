#!/usr/bin/env perl

=pod

=head1 USAGE

	-out   :s => Outfile
	-ext      => extension suffix of files in input directory (defaults to seq219)
	-indir :s => Input directory with *.-ext files (defaults to current dir)
	-append   => Append, don't overwrite outfile

=cut


use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;


my ($csfile,$csdir,$a_if_append);
my $csext = 'seq219';

GetOptions(
	'out:s' => \$csfile,
	'indir:s' => \$csdir,
	'append' => \$a_if_append,
	'ext:s'  => \$csext,
);
pod2usage "Give output csfile\n" unless $csfile;
$csdir = `pwd`  if !$csdir;
chomp($csdir);
my @files = glob($csdir."/*$csext");
print "Checking $csdir with ".scalar(@files)." files\n";

my ($numcsfiles,$num_chars);

    if ($a_if_append) {
	open (IN, "<$csfile.sizes") || die("Error: can't open $csfile.sizes: $!");
	my $line = <IN>;
	close IN;
	$line =~ /(\S*)\s+(\S*)/;
	$numcsfiles = $1;
	$num_chars = $2;
        open (OUT, ">>$csfile");
    } else {
        open (OUT, ">$csfile");
    }
    foreach my $seq219file (@files) { 
	next unless -s $seq219file && !-d $seq219file;
        open (IN, "$seq219file");
        my @lines = <IN>;
        close(IN);
        $seq219file =~ s/.*?([^\/]*)\.$csext\s*/$1/ or die ("Error: $seq219file does not have the extension $csext!?\n");
        foreach my $line (@lines) {
            if ($line =~ /^>/) {
                $line = ">".$seq219file."\n";
                $numcsfiles++;
            } else {
                $num_chars += length($line);
            }
            printf(OUT "%s",$line);     
        }
    } 
    close(OUT); 
    
    open (OUT, ">$csfile.sizes");
    print OUT "$numcsfiles $num_chars\n";
    close OUT;

print "Done\n";
