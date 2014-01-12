#!/usr/local/bin/perl

use FindBin;
use lib ($FindBin::Bin);
use Pasa_init;
use lib ("../PerlLib");
use strict;
use CDNA::Splice_graph_assembler;
use CDNA::CDNA_alignment;
use Getopt::Std;

use vars qw ($opt_h $opt_D $opt_p $opt_d $DEBUG $opt_S $opt_M $opt_s);

&getopts ('hD:dp:S:M:s');


$|=1;
our $SEE = $opt_s;

open (STDERR, "&>STDOUT");

my $usage =  <<_EOH_;

usage: $0 < alignment.textfile

The alignment.textfile should have the following format:

cdna_acc,transcribed_orient,segment_coords,...

The transcribed_orient should be '+|-|?'  
Use '?' in cases where you have single exon alignments of ambiguous transcribed orientation.

ie. 
// cluster: 3111
gi|86081823|gb|DR377582.1|DR377582,-,50113-50337,50419-50631,50883-50924
gi|86081670|gb|DR377427.1|DR377427,-,50109-50337,50419-50631,50883-50963,51126-51187
gi|86081223|gb|DR376980.1|DR376980,-,50109-50337,50419-50631,50883-50926
gi|86078045|gb|DR373802.1|DR373802,-,50113-50337,50419-50631,50883-50963,51064-51108
gi|86077457|gb|DR373214.1|DR373214,-,50113-50337,50419-50631,50883-50963,51126-51155
gi|86076745|gb|DR372502.1|DR372502,-,50177-50337,50419-50631,50883-50963,51126-51194
gi|86074966|gb|DR370723.1|DR370723,-,50444-50631,50883-50963,51126-51183
gi|86049609|gb|DR345364.1|DR345364,-,50236-50337,50419-50631,50883-50963,51126-51191
gi|86049607|gb|DR345362.1|DR345362,-,50920-50963,51126-51187
gi|86049606|gb|DR345361.1|DR345361,-,50920-50963,51126-51187    

// cluster ... etc, etc...

Namely,
accession,orientation,coordinates

_EOH_

    ;

if ($opt_h) { die $usage;}


$/ = "\n//";

while (my $input = <STDIN>) {
    print "######################################################\n";
    print $input . "\n";

    
    my @datalines = split (/\n/, $input);
    my $header = "";
    if ($datalines[0] =~ /\/\//) {
        $header = shift @datalines; #lose the // line
    }
    chomp $header;
    if ($datalines[$#datalines] =~ /\/\//) { pop @datalines;} #rid the last //-containing line.
    my @alignments;

    my %seen;

    foreach my $dataline (@datalines) {
        unless ($dataline =~ /,/) { next;}
        $dataline =~ s/\s+$//; #trim terminal whitespace
        my ($acc, $strand, @coordsets) = split (/,/, $dataline);
        
        if ($seen{$acc}) {
            die "Error, cannot report alignments for acc($acc) multiple times in a single input to pasa.\n";
        }
        $seen{$acc} = 1;
        
        my @segments;
        my $cdna_length = 0;
        foreach my $coordset (@coordsets) {
            my ($lend, $rend) = split (/-/, $coordset);
            my $segment = new CDNA::Alignment_segment($lend, $rend);
            push (@segments, $segment);
            $cdna_length += abs ($rend - $lend) + 1;
        }
        my $alignment = new CDNA::CDNA_alignment($cdna_length, \@segments);
        if ($strand =~ /^[\+\-]$/) {
            $alignment->force_spliced_validation($strand);
        } 
        elsif ($strand eq "?") {
            $alignment->set_spliced_orientation($strand);
        }
        
        
        $alignment->set_acc($acc);
        push (@alignments, $alignment);
    }
    if (@alignments) {
        my $num_alignments = scalar (@alignments);
        print "HEADER: $header\n" if $header;
        my $start_time = time();
        my $assembler = new CDNA::Splice_graph_assembler();
        $assembler->assemble_alignments(@alignments);
        my $end_time = time();
        my $num_seconds = $end_time - $start_time;
        
        # set orientation for display purposes.
        my @assemblies = $assembler->get_assemblies();
        foreach my $assembly (@assemblies) {
            if ((my $orient = $assembly->get_spliced_orientation()) ne '?') {
                $assembly->force_spliced_validation($orient);
            }
            $assembly->remap_cdna_segment_coords();
        }
        
        print $assembler->toAlignIllustration(60);
        print "\n\nTIME: $num_seconds seconds to assemble $num_alignments alignments\n";
        print "\n\n\n";
    
        my $x=0;
        foreach my $assembly ($assembler->get_assemblies()) {
            $x++;
            print "Assembly($x): " .  $assembly->toToken() . "\t" . $assembly->get_acc() . "\n";
        }
    }
}

exit(0);


