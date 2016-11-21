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
open (OUT,">$file.out");
my $counter;
my $errors=int(0);
while (my $id=<IN>){
	my $seq=<IN>;
	my $qid=<IN>;
	my $qua=<IN>;
	$counter++;
	if ($id=~/^\@/ && length($id) < 60){
		unless (length($seq)==length($qua)){
			warn  "qual is not same as seq for ID $id\n";
			next;
		}
		print OUT $id.$seq.$qid.$qua;
	}else{
		next;
		$errors++;
		warn "Error at line " . ($counter*4) . "\n".$id.$seq.$qid.$qua."\n";
		die "Too many errors\n" if $errors>100;
		#one of the other ones has it
		if ($seq=~/^\@MP/){
			print OUT $seq.$qid.$qua.<IN>;
		}elsif ($qid=~/^\@MP/){
			print OUT $qid.$qua.<IN>.<IN>;
                }elsif ($qua=~/^\@MP/){
			print OUT $qua.<IN>.<IN>.<IN>;
		}else{
			warn "Unrecoverable error\n";
			my $skip=<IN>;next;
		}
	}
}

close (IN);
close OUT;
if (-s $file == -s $file.'.out'){
	print "All good!\n";
	unlink($file.'.out');
}else{
	print "Errors have been found and they may have been corrected. See $file.out\n";
}
}
