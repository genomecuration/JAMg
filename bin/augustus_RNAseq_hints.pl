#!/usr/bin/env perl

=pod

=head1 NAME

 augustus_RNAseq_hints.pl

=head1 USAGE

Create hint files for Augustus using RNASeq/EST. One is junction reads (excellent for introns), the other is RNASeq/EST coverage

Will use up to 5Gb of RAM

Mandatory options:

 -bam|in           s{1,}  	The input BAM file(s) (co-ordinate sorted).
 -genome|fasta     s  		The genome assembly FASTA file.

Other options:

 -strandness       i		If RNAseq is directional, provide direction: 0 for unknown (default); or 1 for + strand; -1 for - strand
 -min_score        i  		Minimum score for parsing (defaults to 20)
 -window           i  		Window size for coverage graph (defaults to 50)
 -background_fold  i  		Background (defaults to 4), see perldoc
 -no_hints            		Don't create hints file for Augustus, just process junction reads

=head1 DESCRIPTION

Background: The problem of getting the intron boundary correct is that rnaseq doesn't go to 0 at the intron, but continues at a background level.
For that reason, stop if it is -background_fold times lower than a previous 'good' value


=head1 AUTHORS

 Alexie Papanicolaou

        CSIRO Ecosystem Sciences
        alexie@butterflybase.org

=head1 DISCLAIMER & LICENSE

Copyright 2012-2014 the Commonwealth Scientific and Industrial Research Organization. 
See LICENSE file for license info
It is provided "as is" without warranty of any kind.

=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use List::Util qw(sum);
use Pod::Usage;
use File::Basename;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";

my ( $samtools_exec, $bedtools_exec, $bed_to_aug_script ) = &check_program( 'samtools', 'bedtools','bed12_to_augustus_junction_hints.pl' );



#Options
my ( @bamfiles, $genome, $help,$no_hints );
my $cpus = 4;
my $window           = 50;
my $min_score        = 20;
my $strandness       = int(0);
my $background_level = 4;
pod2usage $! unless &GetOptions(
            'help'              => \$help,
            'bam|in:s{,}'       => \@bamfiles,
            'genome|fasta:s'    => \$genome,
            'min_score:i'       => \$min_score,
            'strandness:i'      => \$strandness,
            'window:i'          => \$window,
            'background_fold:i' => \$background_level,
	    'nohints|no_hints'  => \$no_hints,
	    'cpus:i'            => \$cpus
);

pod2usage if $help;

pod2usage "Cannot find the BAM or genome FASTA file\n"
  unless $bamfiles[0]
   && -s $bamfiles[0]
   && $genome
   && ( -s $genome || -s $genome . '.fai' );

my $strand;
if ( !$strandness || $strandness == 0 ) {
 $strand = '.';
}
elsif ( $strandness > 0 ) {
 $strand = '+';
}
elsif ( $strandness < 1 ) {
 $strand = '-';
}
else {
 die;
}

my $master_bamfile;
if (scalar(@bamfiles == 1)){
  $master_bamfile = $bamfiles[0];
}else{
	foreach my $bamfile (@bamfiles){
		die "Cannot find $bamfile\n" unless -s $bamfile;
	}
	$master_bamfile = 'master_bamfile.bam';
	&process_cmd("$samtools_exec merge -@ $cpus -r $master_bamfile ".join(" ",@bamfiles)." && $samtools_exec index $master_bamfile") unless -s $master_bamfile;
}

&process_cmd("$samtools_exec faidx $genome") unless -s $genome . '.fai';
die "Cannot index genome $genome\n" unless -s $genome . '.fai';
unless (-e "$master_bamfile.junctions.completed"){
 &process_cmd("$samtools_exec rmdup -S $master_bamfile - | $bedtools_exec bamtobed -bed12 | $bed_to_aug_script -prio 7 -out $master_bamfile.junctions.bed"
 ."|sort -S 2G -n -k 4,4 | sort -S 2G -s -n -k 5,5 | sort -S 2G -s -n -k 3,3 | sort -S 2G -s -k 1,1 -o $master_bamfile.junctions.hints" ) unless -s "$master_bamfile.junctions.hints";
 unless ($no_hints){
	 # For Augustus
	 &only_keep_intronic("$master_bamfile.junctions.hints");
	 # don't merge before getting intronic.
	 &merge_hints("$master_bamfile.junctions.hints");
 }
 # For JBrowse
 &process_cmd("$bedtools_exec bedtobam -bed12 -g $genome.fai -i $master_bamfile.junctions.bed| $samtools_exec sort -m 4G -@ 4 -o $master_bamfile.junctions.bam -") unless -s "$master_bamfile.junctions.bam";
 &process_cmd("$samtools_exec index $master_bamfile.junctions.bam");
 &touch("$master_bamfile.junctions.completed");
}

unless (-e "$master_bamfile.coverage.bg.completed"){
 # For JBrowse
 &process_cmd("$bedtools_exec genomecov -split -bg -g $genome.fai -ibam $master_bamfile| sort -S 4G -k1,1 -k2,2n -o $master_bamfile.coverage.bg");
 &process_cmd("bedGraphToBigWig $master_bamfile.coverage.bg $genome.fai $master_bamfile.coverage.bw") if `which bedGraphToBigWig`; 
 &touch("$master_bamfile.coverage.bg.completed");
}

unless (-e "$master_bamfile.coverage.hints.completed" && !$no_hints){
 &bg2hints("$master_bamfile.coverage.bg") ;
 &merge_hints("$master_bamfile.coverage.hints");
 &touch("$master_bamfile.coverage.hints.completed");
}

if (    -e "$master_bamfile.junctions.completed"
     && -e "$master_bamfile.coverage.hints.completed" ){
 unless (-e "$master_bamfile.rnaseq.completed"){
   &process_cmd("cat $master_bamfile.junctions.hints $master_bamfile.coverage.hints"
	."|sort -S 2G -n -k 4,4 | sort -S 2G -s -n -k 5,5 | sort -S 2G -s -n -k 3,3 | sort -S 2G -s -k 1,1 -o $master_bamfile.rnaseq.hints" );
   &merge_hints("$master_bamfile.rnaseq.hints");
   &touch("$master_bamfile.rnaseq.completed");
 }
 print "Done!\n";
}
elsif (!$no_hints) {
 die "Something went wrong....\n";
}else{
  print "Done, no hints were processed as requested\n";
}
###
sub check_program() {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  pod2usage "Error, path to a required program ($prog) cannot be found\n\n"
    unless $path =~ /^\//;
  chomp($path);
  $path = readlink($path) if -l $path;
  push( @paths, $path );
 }
 return @paths;
}
###
sub process_cmd {
 my ($cmd) = @_;
 print "CMD: $cmd\n";
 my $ret = system($cmd);
 if ( $ret && $ret != 256 ) {
  die "Error, cmd died with ret $ret\n";
 }
 return $ret;
}

sub bg2hints() {
 my $bg      = shift;
 my $outfile = $bg;
 $outfile =~ s/.bg$/.hints/;
 open( IN, $bg );
 my ( @array, %area );
 while ( my $ln = <IN> ) {
  chomp($ln);
  my @data = split( "\t", $ln );
  next unless $data[3] >= $min_score;
  # store data in an array
  for ( my $i = $data[1] ; $i <= $data[2] ; $i++ ) {
   # co-ords in bg are 0-based; hints/gff is 1-based
   $area{ $data[0] }{$i+1} = $data[3];
  }
 }

 # print final area
 open( OUT, ">$outfile" );
 foreach my $ref ( sort { $a cmp $b } keys %area ) {
  my @coords = sort { $a <=> $b } ( keys %{ $area{$ref} } );
  for ( my $i = $coords[0] ; $i < @coords ; $i++ ) {
   next if ( !$area{$ref}{$i} );
   my $k = $i + $window;
   $k-- while ( !$area{$ref}{$k} );
   next if $k == $i;
   my @splice;
   
   for ( my $v = $i ; $v <= $k ; $v++ ) {
    my $level = $area{$ref}{$v};
    my $previous_level = $v eq $i ? int(0) : $area{$ref}{$v-1};
    my $next_level = $v eq $k ? 1e6 : $area{$ref}{$v+1};
    # the problem of getting the intron boundary correct is that
    # rnaseq doesn't go to 0 at the intron, but continues at a
    # background level. stop if it is 4 times lower than a previous 'good' value
    if (
    !$level ||
     ( $previous_level  && ( $previous_level > ( $level * $background_level ) ))
     || $next_level && ($level > ( $next_level * $background_level ))
     )
    {
     $k = $v - 1;
     last;
    }
    push( @splice, $level );
   }
#   next if scalar(@splice) < ( $window / 2 );
   my $median = &median( \@splice );
   $median = $splice[0] if !$median;
   next unless $median && $median >= $min_score;
   print OUT $ref
     . "\tRNASeq\texonpart\t"
     . $i . "\t"
     . $k . "\t"
     . $median
     . "\t$strand\t.\tsrc=R;pri=4\n";
   $i +=  $window ;
  }
 }

 close OUT;
 close IN;
 &process_cmd("sort -S 2G -n -k 4,4 $outfile| sort -S 2G -s -n -k 5,5 | sort -S 2G -s -n -k 3,3 | sort -S 2G -s -k 1,1 -o $outfile.");
 rename("$outfile.",$outfile);
 return $outfile;
}

sub only_keep_intronic(){
 my $file = shift;
 my %hash;

 # get introns (we will postprocess them later)
 open (IN,$file);
 open (OUT1,">".$file.".intron.only");
 while (my $ln=<IN>){
   print OUT1 $ln if ($ln=~/\tintron\t/);
 }
 close OUT1;
 close IN;
 # group introns and remove solitary ones
 &merge_hints("$file.intron.only");
 &convert_mult_to_score("$file.intron.only");

 # find exons matching intron 
 open (IN,$file.".intron.only");
 while (my $ln=<IN>){
  my @data = split("\t",$ln);
  if ($data[8]=~/grp=([^;]+)/){
   $hash{$1}=$data[5];
  }
 }
 close IN;

 open (IN,$file);
 open (OUT2,">".$file.".intronic");
 while (my $ln=<IN>){
  my @data = split("\t",$ln);
  if ($data[2] eq 'intron'){
   print OUT2 $ln;
  }elsif ($data[8]=~/grp=([^;]+)/){
   my $grp =$1;
   if ($hash{$grp}){
	$data[5]= $hash{$grp};
	print OUT2 join("\t",@data);
   }
  }
 }
 close IN;
 close OUT2;
 rename($file.".intronic",$file);
}

sub merge_hints(){
 my $file = shift;
 open (IN,$file);
 open (OUT,">$file.merged");
 my (@current_line,@previous_line);
 while (my $ln=<IN>) {
    @current_line = split ("\t",$ln);
    if (!@previous_line){
        @previous_line = @current_line;
    }elsif(($current_line[0] eq $previous_line[0]) && ($current_line[2] eq $previous_line[2]) && 
    (($current_line[3] >= $previous_line[3]) && ($current_line[4] <= $previous_line[4]))
      && ($current_line[6] eq $previous_line[6])){
     # update previous_line by adding current to it
        chomp($previous_line[8]);
        $previous_line[8] =~ s/(grp=[^;]*);*//;
        my $grp = $1;
        $grp .= ';' if $grp;
        $grp = '' if !$grp;
        my ($lm,$m)=(1,1);
        if ($previous_line[8] =~ /mult=(\d+);/){
            $lm = $1;
            $previous_line[8] =~ s/mult=\d+;//;
        }
        if ($current_line[8] =~ /mult=(\d+);/){
            $m = $1;
        }
        $previous_line[8] = "mult=" . ($lm+$m) . ";$grp" . $previous_line[8]."\n";
     
    }elsif (
    !(($current_line[0] eq $previous_line[0]) && ($current_line[2] eq $previous_line[2]) && ($current_line[3] == $previous_line[3]) && ($current_line[4] == $previous_line[4])  && ($current_line[6] eq $previous_line[6]))
    ){
        print OUT join("\t",@previous_line);
        @previous_line = @current_line;
    }
    
     else {
        # update previous_line by adding current to it
        chomp($previous_line[8]);
        $previous_line[8] =~ s/(grp=[^;]*);*//;
        my $grp = $1;
        $grp .= ';' if $grp;
        $grp = '' if !$grp;
        my ($lm,$m)=(1,1);
        if ($previous_line[8] =~ /mult=(\d+);/){
            $lm = $1;
            $previous_line[8] =~ s/mult=\d+;//;
        }
        if ($current_line[8] =~ /mult=(\d+);/){
            $m = $1;
        }
        $previous_line[8] = "mult=" . ($lm+$m) . ";$grp" . $previous_line[8]."\n";
    }
 }
  print OUT join("\t",@previous_line) if (@previous_line);
  close IN;
  close OUT;
  unlink($file);
  rename($file.'.merged',$file);
}

sub touch() {
 my $file = shift;
 system("touch $file");
}

sub mean() {
 return sum(@_) / @_;
}

sub median() {
 my $array_ref = shift;
 my @sorted = sort { $a <=> $b } @{$array_ref};
 return $sorted[ int( @sorted / 2 ) ];
}


sub convert_mult_to_score(){
	my $file = shift;
	return unless $file && -s $file;
	open (IN,$file);
	open (OUT,">$file.scored");
	while (my $ln=<IN>){
		my @data = split("\t",$ln);
		if ($data[8] && $data[8]=~/mult=(\d+)/){
			$data[5] = $1;
			print OUT join("\t",@data);
		}
		
	}
	close IN;
	close OUT;
	rename("$file.scored",$file);
}
