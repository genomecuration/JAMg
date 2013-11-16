#!/usr/bin/env perl

use strict;
use warnings;
use List::Util 'shuffle';

die "Need TRINITY_HOME exported as an environmental variable\n"
  unless $ENV{'TRINITY_HOME'};

print
"Searching for .reads files within a depth of 6 subdirectories and producing new *_trinity_GG commands...\n";
my @files;
@files = `find . -maxdepth 6 -name "*.reads"`;
chomp(@files);
die "No read files found.\n" unless @files && scalar(@files) > 1;

open( OUT,   ">stilltodo.list" );
open( SMALL, ">stilltodo.list.k" );
open( MED,   ">stilltodo.list.m" );
open( LARGE, ">stilltodo.list.g" );
my $minimum_reads = 50;
my $small_cut     = 1024 * 1024 * 1024;
my $medium_cut    = $small_cut * 10;
my $counter;

foreach my $file (@files) {
 my $base = $file;
 $base =~ s/.trinity.reads$//;
 my $sam = $base . '.sam';
 next unless -s $sam;
 unless ( -s $file . ".out.Trinity.fasta" ) {
  my $size = -s $file;
  next unless $size;
  if ( $size < $small_cut ) {
   my $reads = `wc -l < $file`;
   chomp($reads);
   $reads /= 2;
   if ( $reads < $minimum_reads ) {
    warn
"$file: Fewer than $minimum_reads reads. Skipping small readset ($reads reads; $size bytes)\n";
    next;
   }
   else {
    print SMALL $file . "\n";
   }
  }
  elsif ( $size >= $small_cut && $size < $medium_cut ) {
   print MED $file . "\n";
  }
  elsif ( $size >= $medium_cut ) {
   print LARGE $file . "\n";
  }
  print OUT $file . "\n";
  $counter++;
 }
}
close OUT;
close SMALL;
close MED;
close LARGE;

system(
 "rm -f small_trinity_GG.cmds* medium_trinity_GG.cmds* large_trinity_GG.cmds*");

system(
"\$TRINITY_HOME/util/GG_write_trinity_cmds.pl --reads_list_file stilltodo.list.k --paired  > small_trinity_GG.cmds"
) if -s "stilltodo.list.k";
system(
"\$TRINITY_HOME/util/GG_write_trinity_cmds.pl --jaccard_clip --reads_list_file stilltodo.list.m --paired  > medium_trinity_GG.cmds"
) if -s "stilltodo.list.m";
system(
"\$TRINITY_HOME/util/GG_write_trinity_cmds.pl --jaccard_clip --reads_list_file stilltodo.list.g --paired  > large_trinity_GG.cmds"
) if -s "stilltodo.list.g";

unlink("stilltodo.list");
unlink("stilltodo.list.k");
unlink("stilltodo.list.m");
unlink("stilltodo.list.g");

if ( -s "small_trinity_GG.cmds" ) {
 open( IN, "small_trinity_GG.cmds" );
 my @array;
 while ( my $ln = <IN> ) {
  $ln =~ s/JM 2G --CPU 4/JM 1G --CPU 1/;
  chomp($ln);
  $ln .= " >/dev/null\n";
  push( @array, $ln );
 }
 close IN;
 open( OUT, ">small_trinity_GG.cmds." );
 print OUT shuffle(@array);
 close OUT;
 rename( "small_trinity_GG.cmds.", "small_trinity_GG.cmds" );
 system("split -d -a 3 -l 1000 small_trinity_GG.cmds small_trinity_GG.cmds.");
}

if ( -s "medium_trinity_GG.cmds" ) {
 open( IN, "medium_trinity_GG.cmds" );
 my @array;
 while ( my $ln = <IN> ) {
  $ln =~ s/JM 2G --CPU 4/JM 3G --CPU 2/;
  chomp($ln);
  $ln .= " >/dev/null\n";
  push( @array, $ln );
 }
 close IN;
 open( OUT, ">medium_trinity_GG.cmds." );
 print OUT shuffle(@array);
 close OUT;
 rename( "medium_trinity_GG.cmds.", "medium_trinity_GG.cmds" );
 system("split -d -a 3 -l 100 medium_trinity_GG.cmds medium_trinity_GG.cmds.");
}

if ( -s "large_trinity_GG.cmds" ) {
 open( IN,  "large_trinity_GG.cmds" );
 open( OUT, ">large_trinity_GG.cmds." );
 while ( my $ln = <IN> ) {
  $ln =~ s/JM 2G --CPU 4/JM 3G --CPU 6/;
  chomp($ln);
  $ln .= " >/dev/null\n";
  print OUT $ln;
 }
 close IN;
 close OUT;
 rename( "large_trinity_GG.cmds.", "large_trinity_GG.cmds" );
 system("split -d -a 3 -l 1 large_trinity_GG.cmds large_trinity_GG.cmds.");
}
