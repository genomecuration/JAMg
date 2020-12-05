#!/usr/bin/perl -w


=head1 NAME

=head1 USAGE

	*  -i|fa|fasta   =s 	=> FASTA file to split
  	   -size         =i     => How many megabases per file (def 0.5)
           -overlap      =i     => How many bp for overlaps (def 1e4)

=cut

$|=1;		# autoflush
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib ("$RealBin/../PerlLib");
$ENV{PATH} .= ":$RealBin:$RealBin/../3rd_party/bin/";

my ($split_fasta,%ref_data,@hintfiles);
my ($tmpdir,$cpus,$sort_buffer) = ("/tmp",6,'4G');
my ($file2split,$splitfiledir);
my $size_mbp=0.5;
my $overlap_bp = 1e4;
my $suffix="seq";

pod2usage $! unless &GetOptions(
	'i|fa|fasta=s'  => \$file2split,
	'size=i'	=> \$size_mbp,
	'overlap=i'     => \$overlap_bp,
);

$size_mbp = $size_mbp * 1e6;

if (!$file2split || !-s $file2split){pod2usage;}
if (!$splitfiledir){$splitfiledir="$file2split"."_split";}
if ($suffix && $suffix!~/^\./){$suffix=".".$suffix;}
my $label=$file2split;
my ($flag);
my $filecount=0;
my $seqcount=0;
mkdir $splitfiledir;			 # to allow multiple runs of this program in wrapper
my $orig_sep = $/;
$/ =  ">";
$|=1;
open (FILE, "$file2split") || die $!;
print "Processing $file2split\n";


#my %OK = ('Herato2001'=>1,'Herato1301'=>1,'Herato1202'=>1,'Herato0101'=>1) if $debug;


while (my $record=<FILE>){
	chomp($record);
	next unless $record;
	my @lines = split("\n",$record);
	my $id = shift (@lines);
	next unless $id;
	my $label = $id;
	$label =~ s/^(\S+).*$/$1/;
	$label =~ s/\|/_/g;
	$label =~ s/_+$//;
#next unless $OK{$label} if $debug;
	$seqcount++;
	my $seq = join("",@lines);
	$seq=~s/\s+//g;
	if (length($seq) <= $size_mbp){
		$filecount++;
		my $outfile=$label.$suffix;
		open (OUT, ">$splitfiledir\/$outfile");
		print OUT ">".$label."\n".$seq."\n";
		close OUT;		
	}else{
		my @seq_array = split("",$seq);
		for (my $split_start=0;$split_start<scalar(@seq_array);$split_start+=$size_mbp-$overlap_bp){
			my $split_end = int($split_start + $size_mbp) >= scalar(@seq_array) ? int(scalar(@seq_array)-1) : int($split_start + $size_mbp);
			if ($split_start > (scalar(@seq_array)-$overlap_bp)){last;}
			my $split_seq = join('',@seq_array[$split_start..$split_end]);
			my $split_label = $label.'_'.$split_start.'-'.$split_end;
			$filecount++;
			my $outfile=$split_label.$suffix;
			open (OUT, ">$splitfiledir\/$outfile");
			print OUT ">".$split_label."\n".$split_seq."\n";
			close OUT;
		}

	}
	print "$filecount seqs\r";
}
close (FILE);
$/=$orig_sep;
$|=0;
print "\nProcessed $filecount files.\n";


###############3
sub combine_results(){

my $total_count=int(0);
open (FASTA,$split_fasta)||die;
while (my $record=<FASTA>){
	chomp($record);next unless $record;
	my @lines = split("\n",$record);
	my $id = shift (@lines);
	if ($id=~/^(\S+)_(\d+)-(\d+)$/){
		my ($ref,$start,$stop) = ($1,$2,$3);# always + strand
		$ref_data{$ref}{$start}{$stop}++;
		$total_count++;
	}

}
close FASTA;
$/ = $orig_sep;
print "Found $total_count values\n\n";
# ==> master_bamfile2.bam.rnaseq.hints <==
#Herato0101      RNASeq  exonpart        436     469     21      .       .       src=RCOV;pri=4
$|=1;


open (OUT,">$split_fasta.hints");


foreach my $hint_file (@hintfiles){
	print "Parsing hints file $hint_file\n";
	my $line_max = `wc -l < $hint_file`;chomp($line_max);
	my $line_count = int(0);
	open (IN,$hint_file);
	while (my $ln =<IN>){
		$line_count++;
		my @data = split("\t",$ln);
		next unless $data[8];
		if ($data[3] > $data[4]){
			my $t = $data[4];
			$data[4] = $data[3];
			$data[3] = $t;
			$data[6] = '-';
		}
		die "Unexpected start > stop in hints file:\n$ln\n" if $data[3] > $data[4] ;
		my ($ref,$hint_start,$hint_stop) = ($data[0],$data[3],$data[4]);
		# if we have a reference sequence match
		if ($ref_data{$ref}){
			# and the hint is part of our dataset
			foreach my $gene_start (keys %{$ref_data{$ref}}){
				next if $hint_start < $gene_start;
				foreach my $gene_stop (keys %{$ref_data{$ref}{$gene_start}}){
					last if $hint_stop > $gene_stop;
					if ($hint_start >= $gene_start && $hint_stop <= $gene_stop){
						my $opt_hint_start = $hint_start - $gene_start +1;
						my $opt_hint_stop  = $hint_stop - $gene_start +1;
						my @opt_data = @data;
						$opt_data[0] = $ref."_$gene_start"."-$gene_stop";
						$opt_data[3] = $opt_hint_start;
						$opt_data[4] = $opt_hint_stop;
						print OUT join("\t",@opt_data);
					}
				}
			}
		}
		my $prog = int(($line_count / $line_max) * 100);
		print "$prog %  \r" if $line_count % 10000 == 0;
	}
	close IN;
	print "100 % \n";
}
$|=0;

close OUT;
}
######################################
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

sub check_program() {
 my @paths;
 foreach my $prog (@_) {
  my $path = `which $prog`;
  pod2usage "Error, path to a required program ($prog) cannot be found\n\n"
    unless $path =~ /^\//;
  chomp($path);
  #$path = readlink($path) if -l $path;
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

