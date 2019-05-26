#!/usr/bin/perl -w

my @files = @ARGV;

foreach my $file (@files){
	next unless $file && -s $file;
	print "Checking file $file\n";
	&check_file($file);
}

sub check_file(){
  my $file=shift||die;

  open (IN,$file);
  while (my $id=<IN>){
	my $seq=<IN>;
	my $qid=<IN>;
	my $qua=<IN>;
	if ((!$seq || !$qid || !$qua) || $id=~/^\@/){
		unless (length($seq)==length($qua)){
			unlink($file);
			warn  "qual is not same as seq for ID $id\n";
		}
	}
  }
  close (IN);
}
