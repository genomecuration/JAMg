#!/usr/bin/env perl

use strict;
use warnings;

my @files = <tmp.iworm.fa.pid_*.thread_*.fa.pslx>;

my %pid_to_files;
foreach my $file (@files) {
    $file =~ /pid_(\d+)/ or die "error, cannnot get pid from filename $file";
    my $pid = $1;
    push (@{$pid_to_files{$pid}}, $file);
}

foreach my $pid (keys %pid_to_files) {
    my @files = @{$pid_to_files{$pid}};

    print "$pid\t" . scalar(@files) . "\n";
    my $num_threads = scalar(@files);
    
    my %FL_acc_to_file;
    
    foreach my $file (@files) {

        my $FL_file_counter = 0;
                
        my $fl_file = "$file.FL_selected";
        open (my $fh, $fl_file) or die $!;
        while (<$fh>) {
            my @x = split(/\t/);
            my $acc = $x[13];
            
            $FL_acc_to_file{$acc}->{$fl_file} = 1;
            $FL_file_counter++;
        }
        close $fh;

        print "Num threads: $num_threads\t$fl_file\t$FL_file_counter\n";
    }



    my $number_FL = scalar(keys %FL_acc_to_file);
    print "Num threads: $num_threads\tTotal:\n$number_FL\n";

    ## generate distribution of FL counts.
    my %count_counter;
    foreach my $acc (keys %FL_acc_to_file) {
        my @files = keys %{$FL_acc_to_file{$acc}};
        my $num_files = scalar @files;
        $count_counter{$num_files}++;
    }
    print "Histogram of reconstructins across threads.\n";
    my @counts = sort {$a<=>$b} keys %count_counter;
    foreach my $count (@counts) {
        print "$count\t$count_counter{$count}\n";
    }
    print "\n"; # spacer
        
    
}


exit(0);

