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
 -min_coverage     i	  	Minimum coverage for parsing coverage (defaults to 20)
 -window           i  		Window size for coverage graph (defaults to 50)
 -background_fold  i  		Background (defaults to 4), see perldoc
 -no_hints            		Don't create hints file for Augustus, just process junction reads
 -cpus             i		Number of CPUs for sorting (defaults to 4)
 -memory           s            Amount of memory for sorting. Use K/M/G for kilo/mega/giga-bytes (defaults to 5G)
 -min_jr_length    i		Minimum read length to use for junctions (def to 75)
 -min_jr_score     i		Minimum MAQ score to use for junctions (def to 34)
 -min_jr_reads     i            Minimum number of introns to use for print out "known" splice sites database (Augustus and GMAP). Not used for hints. Defaults to 30
 -min_intron       i            Minimum intron length (def 30)

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
use Fasta_reader;

#Options
my ( @bamfiles, $genome, $help,$no_hints, $master_bamfile );

my $cpus = 4;
my $sort_buffer = '5G';
my $tmpdir = $ENV{'TMP'};
$tmpdir = $ENV{'TMPDIR'} if !$tmpdir;
$tmpdir = '/tmp' if !$tmpdir;
my $window           = 50;
my $min_coverage        = 20;
my $min_jr_score        = 34;
my $min_jr_length       = 75;
my $strandness       = int(0);
my $background_level = 4;
my $intron_coverage_cutoff = 30;
my $min_intron_length = 30;

pod2usage $! unless &GetOptions(
            'help'              => \$help,
            'bam|in:s{,}'       => \@bamfiles,
            'genome|fasta:s'    => \$genome,
            'min_coverage:i'    => \$min_coverage,
            'strandness:i'      => \$strandness,
            'window:i'          => \$window,
            'background_fold:i' => \$background_level,
	    'nohints|no_hints'  => \$no_hints,
	    'cpus:i'            => \$cpus,
	    'memory:s'          => \$sort_buffer,
	    'min_jr_length:i'   => \$min_jr_length,
	    'min_jr_score:i'    => \$min_jr_score,
	    'min_intron:i'      => \$min_intron_length,
	    'min_jr_reads:i'    => \$intron_coverage_cutoff,
	    'tmp:s'             => \$tmpdir
);

my ( $samtools_exec, $bedtools_exec, $bed_to_aug_script ) = &check_program( 'samtools', 'bedtools','bed12_to_augustus_junction_hints.pl' );
my $sort_exec = &check_sort_version;

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
 die "Weird strand option given $strandness";
}


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
 my $junction_bam = &grab_intronic_bam($master_bamfile);
 if (-f $junction_bam && (-s $junction_bam) > 10000){
 	&process_cmd("$bedtools_exec bamtobed -bed12 < $junction_bam | $bed_to_aug_script -prio 7 -out $master_bamfile.junctions.bed"
 	."|$sort_exec -n -k 4,4 | $sort_exec -s -n -k 5,5 | $sort_exec -s -n -k 3,3 | $sort_exec -s -k 1,1"
 	." -o $master_bamfile.junctions.all.hints" ) unless (-s "$master_bamfile.junctions.all.hints");
	# i put it here as previous cmd takes ages
	&process_cmd("fg 2>/dev/null;$samtools_exec index $junction_bam.sorted") unless -s "$junction_bam.sorted.bai";

	 unless ($no_hints){
		 # For Augustus
		 &intron_driven_fixes("$master_bamfile.junctions.all.hints");
		 # don't merge before getting intronic.
		 &merge_hints("$master_bamfile.junctions.hints");
	  }
 }else{	
	warn "No junction reads found!";
 }
 &touch("$master_bamfile.junctions.completed");
}

unless (-e "$master_bamfile.coverage.bg.completed"){
 # For JBrowse
 &process_cmd("$bedtools_exec genomecov -split -bg -g $genome.fai -ibam $master_bamfile| $sort_exec -k1,1 -k2,2n -o $master_bamfile.coverage.bg");
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
	."|$sort_exec -n -k 4,4 | $sort_exec -s -n -k 5,5 | $sort_exec -s -n -k 3,3 | $sort_exec -s -k 1,1 -o $master_bamfile.rnaseq.hints" );
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
  next unless $data[3] >= $min_coverage;
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
   next unless $median && $median >= $min_coverage;
   print OUT $ref
     . "\tRNASeq\texonpart\t"
     . $i . "\t"
     . $k . "\t"
     . $median
     . "\t$strand\t.\tsrc=RCOV;pri=4\n";
   $i +=  $window ;
  }
 }

 close OUT;
 close IN;
 &process_cmd("$sort_exec -n -k 4,4 $outfile| $sort_exec -s -n -k 5,5 | $sort_exec -s -n -k 3,3 | $sort_exec -s -k 1,1 -o $outfile.");
 rename("$outfile.",$outfile);
 return $outfile;
}

sub intron_driven_fixes(){
 my $file = shift;
 my $outfile = $file;
 $outfile =~ s/junctions.all.hints/junctions.hints/;

 my %hash;

 # get introns (we will postprocess them later)
 # fix strand of intron based on sequence
 # keep non-canonical ones (singletons purged later)
 my $intronic_file = &get_intron_orient($file);

 # group introns based on co-ords
 &merge_hints($intronic_file);
 # set score to equal number of introns grouped 
 # and remove singletons/solitary ones
 &convert_mult_to_score($intronic_file);

 # find 1 exon pair matching intron that has been grouped
 # and give it the same score AND strand
 open (IN,$intronic_file);
 while (my $ln=<IN>){
  my @data = split("\t",$ln);
  if ($data[8] && $data[8]=~/grp=([^;]+)/){
   $hash{$1} = {'cov' => $data[5],'strand' => $data[6]};
  }
 }
 close IN;

 open (IN,$file);
 open (OUT2,">".$outfile);
 while (my $ln=<IN>){
  my @data = split("\t",$ln);
  # only print those with groups, ie. no solitary ones
  # and introns that met criteria of convert_mult_to_score
  if ($data[8] && $data[8]=~/grp=([^;]+)/){
     my $group = $1;
     next unless exists($hash{$group});
     $data[5] = $hash{$group}{'cov'} if $data[2] ne 'intron';
     $data[6] = $hash{$group}{'strand'};
     print OUT2 join("\t",@data);
  }
 }
 close IN;
 close OUT2;
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
 my $median = $sorted[ int( @sorted / 2 ) ];
 return $array_ref->[0] if !$median;
 return $median;
}


sub convert_mult_to_score(){
	my $file = shift;
	return unless $file && -s $file;
	open (IN,$file);
	open (OUT,">$file.scored");
	while (my $ln=<IN>){
		my @data = split("\t",$ln);
		next unless abs($data[4] - $data[3]) >= $min_intron_length; 
		# discard singletons	
		if ($data[8] && $data[8]=~/mult=(\d+)/){
			$data[5] = $1;
			print OUT join("\t",@data);
		}	
	}
	close IN;
	close OUT;
	rename("$file.scored",$file);
}

sub get_intron_orient(){
	my $hint_file = shift;
	print "Reading genome file $genome\n";
	my $fasta_obj = new Fasta_reader($genome);
	my %fasta_data = $fasta_obj->retrieve_all_seqs_hash();
	# use the introns to set orientation for both introns and exons
	# (this would require parsing the file a couple of time)
	print "Parsing HINT file $hint_file\n";
	open (IN,$hint_file) || die;
	# # spit out splice sites while at it
	open (SPLICE1,">$hint_file.splice.augustus");
	open (SPLICE2,">$hint_file.splice.gmap");
	open (GFF,">$hint_file.intron.only");
	my $printed_counter=int(1);
	my (%intron_count,%intron_printed);
	while (my $ln=<IN>){
		my @data = split("\t",$ln);
		next unless $data[6] && $data[2] eq 'intron';
		chomp($data[8]);
		my ($ref,$type,$start,$end,$strand) = ($data[0],$data[2],$data[3],$data[4],$data[6]);
		die "Cannot find sequence for $ref in genome FASTA file\n" unless $fasta_data{$ref};
		die "Unexpected GFF with start higher than the end\n$ln" if $start > $end;
		# acceptor and donor sites (3' and 5' on + strand)
		my $site1 = uc(substr($fasta_data{$ref},($start-1),2));
		my $site2 = uc(substr($fasta_data{$ref},($end-1-1),2));
		my ($donor,$acceptor,$type1,$type2);
		if (
		(($site1 eq 'GT' || $site1 eq 'GC') && $site2 eq 'AG') ||
		($site1 eq 'AT' && $site2 eq 'AC')
		 ){
			$strand = '+';
			($donor,$acceptor,$type1,$type2) = ($site1,$site2,'dss','ass');
		}elsif (
		(($site2 eq 'AC' || $site2 eq 'GC') && $site1 eq 'CT') ||
		($site2 eq 'AT' && $site1 eq 'GT')
		 ){
			$strand = '-';
			($donor,$acceptor,$type1,$type2) = ($site2,$site1,'ass','dss');
		}else{
			# leave them unchanged
			$data[8] .='noncanonical=true;';
			print GFF join("\t",@data)."\n";
			next;
		}
		$intron_count{"$ref:$start..$end"}++;

		#only if supported by at least intron_coverage_cutoff introns, then save in database
		if (!$intron_printed{"$ref:$start..$end"} && $intron_count{"$ref:$start..$end"} && $intron_count{"$ref:$start..$end"} >= $intron_coverage_cutoff){
			# sequence 40 bp up and 40 bp upstream
			my $site1_seq = lc(substr($fasta_data{$ref},($start-1-40),82));
			my $site2_seq = lc(substr($fasta_data{$ref},($end-1-1-40),82));
			$site1_seq =~ tr/n/-/;$site2_seq =~ tr/n/-/;
			if ($strand eq '-'){
				$site1_seq = &revcomp($site1_seq);
				$site2_seq = &revcomp($site2_seq);
			}
			#augustus
			print SPLICE1 "$type1 $site1_seq\n" if $site1_seq;
			print SPLICE1 "$type2 $site2_seq\n" if $site2_seq;
			#gmap
			my $intron_gsnap_txt = ($strand eq '-') ? ">intron.$printed_counter $ref:$end..$start\n" : ">intron.$printed_counter $ref:$start..$end\n";
			print SPLICE2 $intron_gsnap_txt;
			$intron_printed{"$ref:$start..$end"}++;
			$printed_counter++;
		}
		# print ALL of them
		$data[6] = $strand;
                $data[8].="donor=$donor;acceptor=$acceptor;noncanonical=false;";
		print GFF join("\t",@data)."\n";
	}
	close IN;
	close GFF;
	close SPLICE1;
	close SPLICE2;

	#uniquefy
	system("$sort_exec -u -o $hint_file.splice.augustus. $hint_file.splice.augustus");
	rename("$hint_file.splice.augustus.","$hint_file.splice.augustus");

	system("$sort_exec -u -k2,2 -o $hint_file.splice.gmap. $hint_file.splice.gmap");
	rename("$hint_file.splice.gmap.","$hint_file.splice.gmap");

	return "$hint_file.intron.only";
}

sub revcomp {
  my $dna     = shift;
  return if !$dna;
  my $revcomp = reverse($dna);
  $revcomp =~ tr/ACGTacgt/TGCAtgca/;
  return $revcomp;
}

sub grab_intronic_bam(){
	my $file = shift;
	my $outfile = $file.".junctions.bam";
	return $outfile if -s $outfile;
	# at least 10 bp each side and intron at least 10bp.
	# without \t at end, allows for double introns (cf bed12_to_augustus_junction_hints.pl )
	&process_cmd("$samtools_exec view -@ $cpus -q $min_jr_score -m $min_jr_length $file |grep -P '[0-9]{2,}M[0-9]{2,}N[0-9]{2,}M'|"
		." samtools view -T $genome -@ $cpus -b -o $outfile -");
	&process_cmd("$samtools_exec rmdup -S $outfile - |$samtools_exec sort -m $sort_buffer -o $outfile.sorted - &");
	return $outfile;
}


=pod

=head1 CIGAR Specs

 M 0 alignment match (can be a sequence match or mismatch)
 I 1 insertion to the reference
 D 2 deletion from the reference
 N 3 skipped region from the reference
 S 4 soft clipping (clipped sequences present in SEQ)
 H 5 hard clipping (clipped sequences NOT present in SEQ)
 P 6 padding (silent deletion from padded reference)
 = 7 sequence match
 X 8 sequence mismatch

 H can only be present as the first and/or last operation.
 S may only have H operations between them and the ends of the CIGAR string.
 For mRNA-to-genome alignment, an N operation represents an intron. For other types of alignments,
the interpretation of N is not defined.
 Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ.

=cut

sub check_sort_version(){
	my ($sort_exec) = &check_program('sort');
	my @v=`$sort_exec --version`;

	if ($v[0] && $v[0]=~/(\d+)\.(\d+)\s*$/){
		my $major = $1;
		my $minor = $2;
		if ($major >= 8 && $minor >= 6){
			return "$sort_exec -T $tmpdir --parallel $cpus -S $sort_buffer";
		}else{
			return "$sort_exec -S $sort_buffer -T $tmpdir";
		}
	}else{
		die "Sort of coreutils not found!";
	}
}
