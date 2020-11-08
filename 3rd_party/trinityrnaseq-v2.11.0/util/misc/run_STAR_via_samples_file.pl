#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;
use lib("$FindBin::RealBin/../../PerlLib");
use Pipeliner;
use File::Basename;
use Cwd;
use List::Util qw(min);

use Carp;
use Getopt::Long qw(:config no_ignore_case bundling pass_through);


my $usage = <<__EOUSAGE__;

######################################################################
#
#  Required:
#  --genome <string>           target genome to align to
#  --samples_file  <string>    trinity samples file
#  
#  Optional
#  --gtf <string>              annotations in gtf format
#  --CPU <int>                 number of threads (default: 2)
#  --nameSorted                sort bam by name instead of coordinate
#
#######################################################################


__EOUSAGE__

    ;


my ($genome);
my $samples_file;

my $CPU = 2;

my $help_flag;
my $gtf_file;
my $nameSorted;

&GetOptions( 'h' => \$help_flag,
             'genome=s' => \$genome,
             'samples_file=s' => \$samples_file,
             'CPU=i' => \$CPU,
             'gtf=s' => \$gtf_file,
             'nameSorted' => \$nameSorted,
    );


unless ($genome && $samples_file) {
    die $usage;
}

if ($help_flag) {
    die $usage;
}

if (@ARGV) {
    die "Error, cannot recognize opts: @ARGV";
}


my $star_prog = `which STAR`;
chomp $star_prog;
unless ($star_prog =~ /\w/) {
    die "Error, cannot locate STAR program. Be sure it's in your PATH setting.  ";
}


main: {
	
    ## ensure all full paths
    $genome = &Pipeliner::ensure_full_path($genome);
    $gtf_file = &Pipeliner::ensure_full_path($gtf_file) if $gtf_file;
    

    my $num_contigs = `grep '>' $genome | wc -l`;
    chomp $num_contigs;
    $num_contigs = int($num_contigs);
    unless ($num_contigs > 0) {
        die "Error, couldn't determine the number of contigs in genome: $genome  ... shouldn't happen. ";
    }

    my $genomeChrBinNbits = min(18, int(log((-s $genome) / $num_contigs) / log(2) + 0.5) );
    
    my @read_sets = &parse_samples_file($samples_file);    

    my $pipeliner = new Pipeliner(-verbose => 1);
    my $star_index = "$genome.star.idx";
    my $star_index_chkpt = "$star_index/build.ok";
    ## build star index
    unless (-d $star_index) {
        mkdir($star_index) or die "Error, cannot mkdir $star_index";
    }
    
    my $cmd = "$star_prog --runThreadN $CPU --runMode genomeGenerate --genomeDir $star_index "
        . " --genomeFastaFiles $genome "
        . " --genomeChrBinNbits $genomeChrBinNbits "
        . " --limitGenomeGenerateRAM 40419136213 ";
        
    if ($gtf_file) {

        $cmd .= " --sjdbGTFfile $gtf_file "
            .  " --sjdbOverhang 150 ";
    }            
        
    $pipeliner->add_commands( new Command($cmd, $star_index_chkpt));

    $pipeliner->run();
    
    my $checkpoint_dir = "star_aln_chkpts." . basename($genome);
    unless (-d $checkpoint_dir) {
        mkdir($checkpoint_dir) or die "Error, cannot mkdir $checkpoint_dir";
    }


    my $sort_opt = "SortedByCoordinate";
    my $sort_token = "c";
    my $bam_outfile = "Aligned.sortedByCoord.out.bam";
    if ($nameSorted) {
        $sort_opt = "Unsorted";
        $sort_token = "n";
        $bam_outfile = "Aligned.out.bam";
    }
    
    
    foreach my $read_set_aref (@read_sets) {
        my ($sample_id, $left_fq, $right_fq) = @$read_set_aref;
            
        my $cmd = "$star_prog "
            . " --runThreadN $CPU "
            . " --genomeDir $star_index "
            . " --outSAMtype BAM $sort_opt "
            . " --runMode alignReads "
            . " --readFilesIn $left_fq $right_fq "
            . " --twopassMode Basic "
            . " --alignSJDBoverhangMin 10 "
            . " --outSAMstrandField intronMotif "
            . " --outSAMunmapped Within "
            . " --limitBAMsortRAM 20000000000";
        
            
        if ($left_fq =~ /\.gz$/) {
            $cmd .= " --readFilesCommand 'gunzip -c' ";
        }
        
        $pipeliner->add_commands( new Command($cmd, "$checkpoint_dir/star_align.$sample_id." . basename($genome) . ".ok") );
        
        my $renamed_bam_outfile = "$sample_id.${sort_token}Sorted.star." . basename($genome) . ".bam";
        $pipeliner->add_commands( new Command("mv $bam_outfile $renamed_bam_outfile", "$checkpoint_dir/$renamed_bam_outfile.ok") );

        unless ($nameSorted) {
            $pipeliner->add_commands( new Command("samtools index $renamed_bam_outfile", "$checkpoint_dir/$renamed_bam_outfile.bai.ok") );
        }
        
        
        $pipeliner->run();
    }
    
    
	exit(0);
}



####
sub process_cmd {
	my ($cmd) = @_;
	
	print STDERR "CMD: $cmd\n";
	#return;

	my $ret = system($cmd);
	if ($ret) {
		die "Error, cmd: $cmd died with ret ($ret)";
	}

	return;
}

####
sub parse_samples_file {
    my ($samples_file) = @_;

    my @samples;

    open(my $fh, $samples_file) or die "Error, cannot open file $samples_file";
    while (<$fh>) {
        unless (/\w/) { next; }
        chomp;
        my $line = $_;
        my @x = split(/\t/);
        my ($cond, $rep, $fq_a, $fq_b) = @x;

        unless ($fq_a) {
            confess "Error, line in samples file: $samples_file,  line: [$line] not formatted as expected (sample(tab)replicate(tab)left_fq(tab)right_fq)";
        }
        
        if (! defined $fq_b) {
            $fq_b = "";
        }
                
        $fq_a = &Pipeliner::ensure_full_path($fq_a);
        $fq_b = &Pipeliner::ensure_full_path($fq_b) if $fq_b;
        
        push (@samples, [$rep, $fq_a, $fq_b]);
    }
    close $fh;

    return (@samples);
}

