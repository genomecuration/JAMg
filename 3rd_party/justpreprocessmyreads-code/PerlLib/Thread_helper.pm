package Thread_helper;

use strict;
use warnings;
use Carp;
use threads;

=synopsis


    ## here's how you might use it:

    use threads;
    use Thread_helper;

    my $num_simultaneous_threads = 10;
    my $thread_helper = new Thread_helper($num_simultaneous_threads);

    for (1.1000) {

        $thread_helper->wait_for_open_thread();

        my $thread = threads->create('sleep', int(rand(10)));
        $thread_helper->add_thread($thread);
    }
    $thread_helper->wait_for_all_threads_to_complete();

    my @failures = $thread_helper->get_failed_threads();
    if (@failures) {
        ## examine them...   these are the same threads created above, use use the threads api to access info about them
        ## such as error messages
    }
    else {
        ## all good!
    }


=cut

our $SLEEPTIME = 30;

our $THREAD_MONITORING = 0; # set to 1 to watch thread management


sub new {
    my ($packagename) = shift;
    my ($num_threads,$SLEEPTIME) = @_;
    $SLEEPTIME = 30 if !$SLEEPTIME;
    unless ($num_threads && $num_threads =~ /^\d+$/) {
        confess "Error, need number of threads as constructor param";
    }

    my $self = { 
        max_num_threads => $num_threads,
        current_threads => [],
        
        error_threads => [],
    };

    bless ($self, $packagename);

    return($self);
}


sub add_thread {
    my $self = shift;
    my ($thread,$sleep) = @_;

    my $num_threads = $self->get_num_threads();
    
    if ($num_threads >= $self->{max_num_threads}) {
        
        print STDERR "- Thread_helper: have $num_threads threads running, waiting $sleep sec for a thread to finish...." if $THREAD_MONITORING;
        $self->wait_for_open_thread($sleep);
        print STDERR " done waiting.\n" if $THREAD_MONITORING;
    }
    else {
        print STDERR "- Thread_helper: only $num_threads threads running. Adding another now.\n" if $THREAD_MONITORING;
    }
    
    push (@{$self->{current_threads}}, $thread);

    return;
}


####
sub get_num_threads {
    my $self = shift;
    
    return(scalar @{$self->{current_threads}});
}


####
sub wait_for_open_thread {
    my $self = shift;
    my $sleep = @_;
    $sleep = $SLEEPTIME if !$sleep;

    if ($self->get_num_threads() >= $self->{max_num_threads}) {
        
        my $waiting_for_thread_to_complete = 1;
        
        my @active_threads;
        
        while ($waiting_for_thread_to_complete) {
            
            @active_threads = ();
            
            my @current_threads = @{$self->{current_threads}};
            foreach my $thread (@current_threads) {
                if ($thread->is_running()) {
                    push (@active_threads, $thread);
                }
                else {
                    $waiting_for_thread_to_complete = 0;
                    $thread->join() if $thread->is_joinable;
                    if (my $error = $thread->error()) {
                        my $thread_id = $thread->tid;
                        print STDERR "ERROR, thread $thread_id exited with error $error\n";
                        $self->_add_error_thread($thread);
                    }
                }
            }
            if ($waiting_for_thread_to_complete) {
                sleep($sleep) if $sleep && $sleep >0; 
            }
        }
        
        @{$self->{current_threads}} = @active_threads;
        
        
    }

    return;


}



####
sub wait_for_all_threads_to_complete {
    my $self = shift;
    my @current_threads = @{$self->{current_threads}};
    foreach my $thread (@current_threads) {
        while ($thread->is_running){sleep($SLEEPTIME);}
        $thread->join() if $thread->is_joinable;
        if (my $error = $thread->error()) {
            my $thread_id = $thread->tid;
            print STDERR "ERROR, thread $thread_id exited with error $error\n";
            $self->_add_error_thread($thread);
        }
    }
    
    @{$self->{current_threads}} = (); # clear them out.
    
    return;

}

####
sub get_failed_threads {
    my $self = shift;

    my @failed_threads = @{$self->{error_threads}};
    return(@failed_threads);
}


############################
## PRIVATE METHODS #########
############################


####
sub _add_error_thread {
    my $self = shift;
    my ($thread) = @_;

    push (@{$self->{error_threads}}, $thread);

    return;
}


1; #EOM
