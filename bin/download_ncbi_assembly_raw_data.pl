#!/usr/bin/env perl

=pod

=head1 NAME

download_ncbi_assembly_raw_data.pl

=head1 USAGE

 -file|in|assembly_txt	:s{,} 	=> 1+ assembly report files from NCBI

 One or more of:
  -do_genome			=> Download and process genome data
  -do_dna			=> Download and process DNASeq data
  -do_rna			=> Download and process RNASeq data
  -do_align			N/A
  -do_otherdna                  N/A

 Optional
  -gmap_db_dir		:s	=> Base directory for GMAP databases. Def. to ~/databases/gmap/
  -just_download		=> Just download SRA NCBI files to ~/ncbi/public/sra
  -no_cleanup			=> Don't delete downloaded SRA files

 
=head1 DESCRIPTION

 NB: you need the following software: ncbi-sra-toolkit, blue, JAMg
     a correct Aspera Connect NCBI installation (https://www.ncbi.nlm.nih.gov/books/NBK242625/)


Assembly report files from NCBI File looks like this:

 # Assembly name:  oyster_v9
 # Organism name:  Crassostrea gigas (Pacific oyster)
 # Infraspecific name:  strain=05x7-T-G4-1.051#20
 # Taxid:          29159
 # BioSample:      SAMN00690713
 # BioProject:     PRJNA276446
 # Submitter:      BGI-Shenzhen
 ...

 You can get it from the FTP site
=cut


use strict;
use warnings;
use Carp;
use Data::Dumper;
our $VERSION = '0.1';
use File::Basename;
#use Time::localtime;
use File::stat;
use Pod::Usage;
use Getopt::Long;
use HTML::Strip;
use Cwd qw/abs_path getcwd/;

#SRA
my ($epost_exec, $elink_exec, $efetch_exec,$xtract_exec,$prefetch_exec,$fastq_dump_exec)
 =  &check_programs('epost','elink','efetch','xtract','prefetch','fastq-dump');

#JAMg
die "JAMG_PATH environmental variable missing\n" if !$ENV{'JAMG_PATH'};
$ENV{'PATH'} .= ":" . $ENV{'JAMG_PATH'}."/bin:".$ENV{'JAMG_PATH'}."/3rd_party/bin:".$ENV{'JAMG_PATH'}."/3rd_party/justpreprocessmyreads/";;
my ($trim_fasta_all_exec,$align_dnaseq_gsnap_exec,$align_rnaseq_gsnap_exec,$cleanup_iupac_exec,$N50stats_exec,$gff3_splicesites_exec,$iit_store_exec)
 = &check_programs('trim_fasta_all.pl','align_dnaseq_gsnap.pl','align_rnaseq_gsnap.pl','cleanup_iupac.pl','N50stats.pl','gff3_splicesites','iit_store');
my ($rename_SRA_trimmed_fastq_exec,$auto_sra_preprocess_reads_exec,$create_features_exec,$gt_exec) 
= &check_programs('rename_SRA_orig_fastq.pl','preprocess_illumina_automatic.pl','create_features_from_gff3.pl','gt');

#blue
die "BLUE_PATH environmental variable missing\n" if !$ENV{'BLUE_PATH'};
die "\$BLUE_PATH/Tessel not found\n" unless $ENV{BLUE_PATH}."/Tessel.exe";
my ($mono_exec) = &check_programs('mono');
my $tessel_exec = "$mono_exec ". $ENV{BLUE_PATH}."/Tessel.exe";

#DEW incompatible samtools for time being
#die "DEW_PATH environmental variable missing\n" if !$ENV{'DEW_PATH'};
#$ENV{'PATH'} .= ":" . $ENV{'DEW_PATH'};
my ($salmon_exec);

my (@assembly_txt_files,$just_download,$no_cleanup,$do_rna,$do_genome,$do_dna,$do_align,$do_dna_other,$do_kmers,$do_diffexpr);
my ($just_fetch_metadata,$fetch_ftps,$do_metadata_from_taxid);
# set by NCBI SRA
my $ncbi_download_path = $ENV{'HOME'}."/ncbi/public/sra/";
my $load_modules;
my $cwd = getcwd().'/';
my $gmap_db_dir = $ENV{'HOME'}."/databases/gmap/";
my $salmon_db_dir = $ENV{'HOME'}."/databases/salmon/";
my $ploidy = 2;
GetOptions(
	'file|in|assembly_txt_files:s{,}' => \@assembly_txt_files,
	'just_download' => \$just_download,
	'just_fetch_metadata' => \$just_fetch_metadata,
	'no_cleanup' => \$no_cleanup,
	'gmap_db_dir:s' => \$gmap_db_dir,
	'salmon_db_dir:s' => \$salmon_db_dir,
	'dorna|do_rna|rnaseq' => \$do_rna,
	'dodna|do_dna|dnaseq' => \$do_dna,
	'dootherdna|do_otherdna|otherdnaseq' => \$do_dna_other,
	'dogenome|do_genome|genome' => \$do_genome,
	'doalign|do_align|align' => \$do_align,
	'dodiffexpr|do_diffexpr|diffexpr' => \$do_diffexpr,
	'dokmers|do_kmers|kmers' => \$do_kmers,
	'ploidy:i'		=> \$ploidy,
	'do_taxid|taxid_fetch_metadata:i' => \$do_metadata_from_taxid,
	'fetch_ftps'            => \$fetch_ftps
);

mkdir($gmap_db_dir) if !-d $gmap_db_dir;
mkdir($salmon_db_dir) if !-d $salmon_db_dir;


for (my $i=0;$i< (@assembly_txt_files);$i++){
	unless ($assembly_txt_files[$i] && -f $assembly_txt_files[$i] && -s $assembly_txt_files[$i]){
		warn "SKIP: can't find ".$assembly_txt_files[$i]."\n";
		delete($assembly_txt_files[$i]);
		next;
	}
	$assembly_txt_files[$i] = abs_path($assembly_txt_files[$i]);
}
@assembly_txt_files = sort grep { $_ ne '' } @assembly_txt_files;

pod2usage unless scalar(@assembly_txt_files) && ($do_dna || $do_rna || $do_genome || $do_kmers || $do_align || $do_dna_other || $do_metadata_from_taxid);

if ($do_metadata_from_taxid){
	&get_metadata_from_taxid($do_metadata_from_taxid);
}

if ($do_genome || $do_align){
	($salmon_exec) = &check_programs('salmon');
	foreach my $assembly_txt_file (@assembly_txt_files){
		my ($assembly_path,$genome_fasta,$genome_dbname,$genome_gff,$data_hash_ref) = &prelim_processing($assembly_txt_file);
		unless ($genome_fasta && $genome_dbname && $genome_gff && -s $genome_fasta && -s $genome_gff){
			warn "SKIP: Weird but the FASTA or GFF file are missing. Skipping\n";
		}
		die "Cannot find Genome FASTA ($genome_fasta)\n" unless -s $genome_fasta;
		&process_genome($genome_fasta,$genome_dbname,$genome_gff,$data_hash_ref->{'Organism name'},$data_hash_ref->{'Assembly name'});
	}
}

if ($do_dna){
	foreach my $assembly_txt_file (@assembly_txt_files){
		my ($assembly_path,$genome_fasta,$genome_dbname,$genome_gff,$data_hash_ref) = &prelim_processing($assembly_txt_file);

		my $biosample = $data_hash_ref->{'BioSample'};
		unless ($biosample){
			warn "SKIP: Weird but there is No BioSample for this assembly! Skipping\n";
			next;
		}
		#we assume every assembly has only one biosample ID (which may be a mixture of the other biosamples)
		#this changes into subdir of a biosample ID. the list is inside there
		mkdir('gDNA') unless -d 'gDNA';
		chdir('gDNA');
		my $accession_sra_list_file = &fetch_sra_accessions_biosample($biosample,'GENOMIC','WGS');
		next if $just_fetch_metadata || !$accession_sra_list_file || !-s $accession_sra_list_file;
		mkdir("$assembly_path/gDNA/$biosample/") unless -d "$assembly_path/gDNA/$biosample";
		chdir("$assembly_path/gDNA/$biosample/");
		&download_sras($accession_sra_list_file);
		next if $just_download;
		my $check = &convert_sras($accession_sra_list_file);
		next if $check && $check eq 'FAIL';
		$check = &process_sras($genome_fasta,$genome_dbname);
		next if $check && $check eq 'FAIL';
	}
}

if ($do_rna){
	foreach my $assembly_txt_file (@assembly_txt_files){
		my ($assembly_path,$genome_fasta,$genome_dbname,$genome_gff,$data_hash_ref) = &prelim_processing($assembly_txt_file);
	
		my $taxid = $data_hash_ref->{'Taxid'};
		if (!$taxid){
			warn "SKIP: WEIRD but there is no TaxID for this file. Skipping\n";
			next;
		}
		my @biosample_list = &fetch_biosample_taxid($taxid);
		mkdir('other') unless -d 'other';
		chdir('other');
		mkdir('RNA-Seq') unless -d 'RNA-Seq';
		foreach my $biosample (@biosample_list){
			chdir($assembly_path.'/other/');
			my $accession_sra_list_file = &fetch_sra_accessions_biosample($biosample,'TRANSCRIPTOMIC','RNA-Seq','Illumina');
			next if $just_fetch_metadata || !$accession_sra_list_file || !-s $accession_sra_list_file;
			mkdir("$assembly_path/other/RNA-Seq/$biosample") unless -d "$assembly_path/other/RNA-Seq/$biosample";
			chdir("$assembly_path/other/RNA-Seq/$biosample");
			&download_sras($accession_sra_list_file);
			next if $just_download;
			my $check = &convert_sras($accession_sra_list_file);
			next if $check && $check eq 'FAIL';
			$check = &process_sras();
			next if $check && $check eq 'FAIL';
		}
	}
}

if ($do_dna_other){
	foreach my $assembly_txt_file (@assembly_txt_files){
		my ($assembly_path,$genome_fasta,$genome_dbname,$genome_gff,$data_hash_ref) = &prelim_processing($assembly_txt_file);
		my $taxid = $data_hash_ref->{'Taxid'};
		if (!$taxid){
			warn "SKIP: WEIRD but there is no TaxID for this file. Skipping\n";
                        next;
	        }
		my @biosample_list = &fetch_biosample_taxid($taxid);
		mkdir('other') unless -d 'other';
		chdir('other');
		mkdir('gDNA') unless -d 'gDNA';
		foreach my $biosample (@biosample_list){
			next if $biosample eq $data_hash_ref->{'BioSample'};
			chdir($assembly_path.'/other/');
			my $accession_sra_list_file = &fetch_sra_accessions_biosample($biosample,'GENOMIC','WGS','Illumina');
			next if $just_fetch_metadata || !$accession_sra_list_file || !-s $accession_sra_list_file;
			mkdir("$assembly_path/other/gDNA/$biosample") unless -d "$assembly_path/other/gDNA/$biosample";
			chdir("$assembly_path/other/gDNA/$biosample");
			&download_sras($accession_sra_list_file);
			my $check = &convert_sras($accession_sra_list_file);
			next if $check && $check eq 'FAIL';
			$check = &process_sras($genome_fasta,$genome_dbname);
			next if $check && $check eq 'FAIL';
		}
	}
}

if ($do_align){
#not tested
	# what do we want to do here?
	# we want to align d/rna reads for jbrowse => we use gsnap or gmap for that
		# coverage
		# intron splice
	foreach my $assembly_txt_file (@assembly_txt_files){
		my ($assembly_path,$genome_fasta,$genome_dbname,$genome_gff,$data_hash_ref) = &prelim_processing($assembly_txt_file);

		if (-d $assembly_path.'gDNA'){
			chdir($assembly_path.'gDNA');
			my @dirs = glob("./S*");
			mkdir('aligns') if !-d 'aligns';
			chdir('aligns');
			foreach my $dir (@dirs){
				$dir = '../'.$dir;
				next unless $dir && -d $dir;
				&perform_dna_align($genome_fasta,$genome_dbname,$dir);
			}
		}
		if (-d $assembly_path.'/other/gDNA'){
			chdir($assembly_path.'/other/gDNA');
			my @dirs = glob("./S*");
			mkdir('aligns') if !-d 'aligns';
			chdir('aligns');
			foreach my $dir (@dirs){
				$dir = '../'.$dir;
				next unless $dir && -d $dir;
				&perform_dna_align($genome_fasta,$genome_dbname,$dir);
			}
		}
		if (-d $assembly_path.'/other/RNA-Seq'){		
			chdir($assembly_path.'/other/RNA-Seq');
			my @dirs = glob("./S*");
			mkdir('aligns') if !-d 'aligns';
			chdir('aligns');
			foreach my $dir (@dirs){
				$dir = '../'.$dir;
				next unless $dir && -d $dir;
				&perform_rna_align($genome_fasta,$genome_dbname,$dir);
			}
		}
	}
}

if ($do_diffexpr){
	# we want to align rna reads for DE estimation => we use salmon for that
#not done yet
	($salmon_exec) = &check_programs('salmon');
	foreach my $assembly_txt_file (@assembly_txt_files){
		my ($assembly_path,$genome_fasta,$genome_dbname,$genome_gff,$data_hash_ref) = &prelim_processing($assembly_txt_file);
		next unless -s "$genome_gff.tidy.mRNA" && -d "$salmon_db_dir/$genome_dbname";
		if (-d $assembly_path.'/other/RNA-Seq'){
			chdir($assembly_path.'/other/RNA-Seq');
			my @dirs = glob("./S*");
			mkdir('diffexpr') if !-d 'diffexpr';
			chdir('diffexpr');
			foreach my $dir (@dirs){
				$dir = '../'.$dir;
				next unless $dir && -d $dir;
				&perform_diffexpr_salmon('quant',$genome_gff.'.tidy.mRNA',"$salmon_db_dir/$genome_dbname",$dir);
			}
		}
	}
}

if ($do_kmers){
	foreach my $assembly_txt_file (@assembly_txt_files){
		my ($assembly_path,$genome_fasta,$genome_dbname,$genome_gff,$data_hash_ref) = &prelim_processing($assembly_txt_file);
		my $biosample = $data_hash_ref->{'BioSample'};
		unless (-d $assembly_path.'/gDNA/'.$biosample){
			warn "SKIP: Biosample gDNA/$biosample has not been processed\n";
			next;
		}
		chdir($assembly_path.'/gDNA/'.$biosample);
		&process_kmers($genome_fasta,$ploidy,25);

		if (-d $assembly_path.'/other/gDNA/'){
			chdir($assembly_path.'/other/gDNA/');
			my @dirs = glob("./S*");
			foreach my $dir (@dirs){
				next unless -d $dir;
				chdir($assembly_path.'/other/gDNA/'.$dir);
				&process_kmers($genome_fasta,$ploidy,25);
			}
		}
	}
}


############################################################################
sub process_kmers(){
        my ($fasta,$local_ploidy,$kmer_size) = @_;
	#can change local ploidy if needed
        return if !$fasta;
	my @fastqs = glob("*trimmomatic");
	return unless $fastqs[0] || -s 'kmers_'.$kmer_size.'_histo.txt'; 
	my $dir = getcwd();
	my $sample_dir = basename($dir);
	$local_ploidy = 2 if !$local_ploidy || $local_ploidy < 1;
	print "Kmer-izing data for $sample_dir (".scalar(@fastqs)." files at $dir)\n";
        my $kmer_table_size = 2 * int(-s $fasta); #overkill but otherwise it can crash
	my ($reads_number,$read_tiles);
        &process_cmd("$tessel_exec -k $kmer_size -g $kmer_table_size -t 6 -tmp ".$ENV{'TMP'}." -m 4 -canonical kmers *.trimmomatic") unless -s 'kmers_'.$kmer_size.'_histo.txt';
	
        open (IN,'kmers_'.$kmer_size.'_histo.txt') || die "Cannot produce histogram file ".'kmers_'.$kmer_size.'_histo.txt'."\n";
       	open (OUT,'>kmers_'.$kmer_size.'_histo.csv') if (!-s 'kmers_'.$kmer_size.'_histo.csv');

	my ($check_curve,$check_prev_point);
        while (my $ln=<IN>){
		if ($ln=~/(\d+)\s+reads/){
			$reads_number = $1;
		}elsif($ln=~/(\d+)\s+mers tiled from reads/){
			$read_tiles = $1;
			
		}else{
	               next if $ln=~/^\s*$/ || $ln=~/[^0-9\s\.]/;
       		       print OUT $ln if (!-s 'kmers_'.$kmer_size.'_histo.csv');
			my @data = split("\t",$ln);
			if ($data[0] < 100 and $data[2] > 1){
				$check_curve++ if $check_prev_point && $data[1] && $check_prev_point < $data[1];
				$check_prev_point = $data[1];
			}
		}
        }
       	close OUT if (!-s 'kmers_'.$kmer_size.'_histo.csv');
        close IN;
	my $av_read_length = int($read_tiles /$reads_number) + 1;
	print "\t average read length: $av_read_length\n";
	if (!$check_curve){
		warn "SKIP: histogram does not seem suitable for kmer-based genome size estimation\n";
		return;
	}
	&estimate_genome_coverage('kmers_'.$kmer_size.'_histo.csv',$local_ploidy,$av_read_length,$kmer_size);
}

sub estimate_genome_coverage(){
	# from allpathslg KmerSpectra.cc
	my ($kmer_file,$local_ploidy,$av_read_length,$kmer_size) = @_;
	my $kf_min1_arg = 2;	 #user provided in allpaths
	my $max ;
	$av_read_length = 150 if !$av_read_length || $av_read_length < 1;
	my (@ndk,@cndk,@nk,@cnk);
	# set element 0 so index is matches coverage
	push(@ndk,int(0));
	push(@cndk,int(0));
	push(@nk,int(0));
	push(@cnk,int(0));
	open (IN,$kmer_file) || return;
	while (my $ln=<IN>){
		chomp($ln);
		my @data = split("\t",$ln);
		next unless $data[2];
		# cov	ukmers	cov*ukmers	prop
		# 1       794452826       794452826       9.09
		# 2       65227991        130455982       1.49
		last if $max && $data[0] > $max;
		push(@cndk,$cndk[-1]+$data[1]); #cumulative
		push(@ndk,$data[1]);

		push(@cnk,$cnk[-1]+$data[2]);	#cumulative
		push(@nk,$data[2]);
	}
	close IN;
	my $kf_ceil = scalar(@ndk)-1;
	my $kf_min1 = $kf_min1_arg;
	if ($kf_min1 > $kf_ceil){
		warn "SKIP: error with kmerisation.\n";
		return;
	}
	while ($kf_min1 - 1 >= 2 && $nk[$kf_min1 - 1] < $nk[$kf_min1]){
	    $kf_min1--;
	}
	while ($kf_min1 <= $kf_ceil && $nk[$kf_min1 + 1] < $nk[$kf_min1]){
	    $kf_min1++;
	}
	if ($kf_min1 > $kf_ceil) {
		warn "SKIP: Could not find kmer spectrum local minimum.\n";
		return;
	}
	my $kf_max2 = $kf_min1;
	  for (my $kf = $kf_min1 + 1; $kf < 0.8 * $kf_ceil ; $kf++){
	   if ($nk[$kf] > $nk[$kf_max2]){
	      $kf_max2 = $kf;
	   }
          }

	  if ($kf_max2 < $kf_min1_arg) {
	    warn "SKIP: Coverage too low. Can't estimate genome size.\n";
	    return;
	  }
	 if ($local_ploidy == 2) {
    		my $ndk_half   = $ndk[$kf_max2 / 2];
		my $ndk_double = $ndk[$kf_max2 * 2];
    		if ($ndk_double > $ndk_half){ $kf_max2 *= 2;}
	 }
	my $kf_max1 = $kf_max2 / 2;
  	my $kf_min2 = int($kf_max1 * (2 * $ndk[$kf_max1] + $ndk[$kf_max2]) / ($ndk[$kf_max1] + $ndk[$kf_max2]));
#  	my $kf_min2 = int($kf_max2 / ln2);
			

	for (my $kf = $kf_min1 + 1; $kf < $kf_max1; $kf++){
		if ($nk[$kf] < $nk[$kf_min1]){
	      		$kf_min1 = $kf;
		}
	}

	if ($local_ploidy == 1) {
		$kf_min2 = $kf_min1;
		$kf_max1 = $kf_min2;
	}

	my $kf_min3 = $kf_max2 * 3 / 2;
	if ($kf_min3 > $kf_ceil) {
	    warn "SKIP:  Can't estimate genome size.\n";
 	    return;
	}
	#allpaths:
	#my $kf_hi = (2 * $kf_max2 < scalar(@ndk)-1) ? $kf_max2 * sqrt(4 * $ndk[2 * $kf_max2] * $kf_max2) : $kf_max2 * sqrt(4 * $ndk[-1] * $kf_max2);
	my $kf_hi = (2 * $kf_max2 < scalar(@ndk)-1) ? sqrt(4 * $ndk[2 * $kf_max2] * $kf_max2) : sqrt(4 * $ndk[-1] * $kf_max2);
	$kf_hi = int($kf_hi);
	$kf_hi = $kf_ceil if ($kf_hi > $kf_ceil);

	# ---- number of read kmers in the various intervals
	my $nk_total       = $cnk[-1];
	my $nk_bad_low_kf  = $cnk[$kf_min1];
	my $nk_good_snp    = $cnk[$kf_min2] - $cnk[$kf_min1];
	my $nk_good_uniq   = $cnk[$kf_min3] - $cnk[$kf_min2];
	my $nk_good_rep    = $cnk[$kf_hi] - $cnk[$kf_min3];
	my $nk_bad_high_kf = $nk_total  - $cnk[$kf_hi];

	# ---- number of distinct kmers in the various intervals
	my $ndk_total       = $cndk[-1];
	my $ndk_bad_low_kf  = $cndk[$kf_min1];
	my $ndk_good_snp    = $cndk[$kf_min2] - $cndk[$kf_min1];
	my $ndk_good_uniq   = $cndk[$kf_min3] - $cndk[$kf_min2];
	my $ndk_good_rep    = $cndk[$kf_hi]   - $cndk[$kf_min3];
	my $ndk_bad_high_kf = $ndk_total     - $cndk[$kf_hi];

	my $kf_ave_uniq = sprintf("%.2f",($nk_good_uniq / $ndk_good_uniq));

	my $genome_size            = ($nk_total - $nk_bad_low_kf - $nk_bad_high_kf) / $kf_ave_uniq;
	my $genome_size_unique     = $ndk_good_uniq + $ndk_good_snp / 2;
	my $genome_size_repetitive = $genome_size - $genome_size_unique;
	my $coverage = ($genome_size ? $nk_total / $genome_size : 0);
	my $d_SNP = 0;
	my $stddev_bias = 0;

	if ($coverage < 20) {
	    $genome_size = 0;
	    $genome_size_unique = 0;
	    $genome_size_repetitive = 0;
	    $d_SNP = 0;
	    $stddev_bias = 0;
	}else {
		if ($local_ploidy == 2) {
     			$d_SNP = ($ndk_good_snp > 0) ? 1.0 / (1.0 - ((1.0 - 0.5 * $ndk_good_snp / $genome_size) ** (1.0 / $kmer_size))) : 1000000;
			$d_SNP = int($d_SNP);
		}

		my $sig2_SNP = 0;
		if ($local_ploidy == 2) {
			my $alpha = sprintf("%.2f",($ndk_good_snp / $ndk_good_uniq));
		        my $tmp = sprintf("%.2f",(1.0 / ($alpha + 2.0)));
		        $sig2_SNP =  sprintf("%.2f",($alpha * $tmp * $tmp));
		}

		my $sum = int(0);
		my $mu = int(0);
		my $sig2 = int(0);
		my $sig = int(0);
		$kf_min3 = $kf_max2 + ($kf_max2 - $kf_min1);
		for (my $kf = $kf_min1; $kf != $kf_min3; $kf++){ $sum += $ndk[$kf];}
		for (my $kf = $kf_min1; $kf != $kf_min3; $kf++){ $mu += $kf * $ndk[$kf];}
		$mu /= $sum;
		for (my $kf = $kf_min1; $kf != $kf_min3; $kf++){ $sig2 += ($kf - $mu) * ($kf - $mu) * $ndk[$kf];}
		$sig2 /= $sum;
		$sig2 = sprintf("%.2f",abs($sig2));
		$sig = sprintf("%.2f",(sqrt($sig2)));
     		my $sig2_bias = sprintf("%.2f",( ($sig2 - $mu)/($mu * $mu) - $sig2_SNP  ));
      		$stddev_bias = sprintf("%.2f",(sqrt($sig2_bias) ));
	}
	if (!$genome_size || $genome_size < 10){
		warn "SKIP: Genome kmerization failed\n";
		return;
	}
	$genome_size = &thousands(int($genome_size));
	$genome_size_unique = &thousands(int($genome_size_unique));
	$genome_size_repetitive = &thousands(int($genome_size_repetitive));
	print "INFO:\tGenome size is $genome_size bp\n\tGenome size (unique) is $genome_size_unique bp\n\tGenome size (repeats) is $genome_size_repetitive bp\n";
	print "\tSNPs are every $d_SNP bp on average (stddev bias: $stddev_bias)\n";
	print "\tKmer cutoffs were:\n\t\tlower,upper cutoff for good heterozygous kmers: $kf_min1,$kf_min2\n"
		."\t\tupper cuttoff for non-repetitive kmers: $kf_min3\n"
		."\t\tabsolute heterozygous maximum: $kf_max1\n"
		."\t\tabsolute homozygous maximum: $kf_max2\n"
		."\t\tignored repeats over: $kf_hi\n\n";
}

sub perform_diffexpr_salmon(){
	my ($salmon_type,$mrna_file,$index_dir,$input_dir) = @_;
	#3'UTR quantseq needs --noLengthCorrection and no --seqBias --gcBia

	if (!-s "$mrna_file.map"){
		open (IN,$mrna_file);
		open (OUT,">$mrna_file.map");
		while (my $ln=<IN>){
			if ($ln=~/^>(\S+)\s.+gene:(\S+)/){
				print OUT "$1\t$2\n";
			}
		}
		close IN;
		close OUT;
	}
	#output should be SAMN not SRR RNA-Seq/SAMN07333176/SRR
	my $common_cmd = "$salmon_exec --no-version-check quant --libType A --seqBias --gcBias --fldMax 800 --fldMean 400 --geneMap $mrna_file.map --auxDir salmon.info ";
	if ($salmon_type=~/ali/i){
#		&process_cmd($common_cmd." --alignments in.BAMs -t $mrna_file --threads 6 --sampleOut out.bam --output directory ");
	}else{
#		&process_cmd($common_cmd." --index $index_dir --threads 10 --mates1 file(s) --mates2 file(s) --output file ");
	}

}

sub perform_dna_align(){
#TOSDO: preprocess the ../samplename.text
#GENOMIC WGS     PacBio RS II

	#TODO mate pair detection?
	# max length is 300 for gsnap: now automatically done by align helper script
	my ($fasta,$dbname,$dir) = @_;
	if (!-d "$gmap_db_dir/$dbname" || !-s "$fasta.trim.x"){
		warn "SKIP: I cannot perform alignments if the genome is not preprocessed as .trim.x and indexed at $gmap_db_dir/$dbname.\n";
		return;
	}
	my $cmd = "$align_dnaseq_gsnap_exec -fasta $fasta.trim.x -dbname $dbname -gmap_dir $gmap_db_dir -suffix -input_dir $dir ";
	$cmd .= " -large_genome " if -s $fasta > 3000000000;
	$cmd .= " -filetype trimmomatic -nofail -pattern1 _1_fastq "; 
	if (()=glob("$dir/*_2_fastq*trimmomatic")){
		$cmd .= " -pattern2 _2_fastq ";
	}else{
		$cmd .= " -notpaired ";
	}
	&process_cmd($cmd);
}

sub perform_rna_align(){
#TOSDO: preprocess the ../samplename.text
#TRANSCRIPTOMIC  miRNA-Seq       Illumina Genome Analyzer IIx

	# max length is 300 for gsnap: now automatically done by align helper script
	my ($fasta,$dbname,$dir) = @_;
	if (!-d "$gmap_db_dir/$dbname" || !-s "$fasta.trim.x"){
		warn "SKIP: I cannot perform alignments if the genome is not preprocessed as .trim.x and indexed at $gmap_db_dir/$dbname.\n";
		return;
	}
	my $cmd = "$align_rnaseq_gsnap_exec -fasta $fasta.trim.x -dbname $dbname -gmap_dir $gmap_db_dir -input_dir $dir ";
	$cmd .= " -large_genome " if -s $fasta > 3000000000;
	$cmd .= " -filetype trimmomatic -nofail -pattern1 _1_fastq "; 
	if (()=glob("$dir/*_2_fastq*trimmomatic")){
		$cmd .= " -pattern2 _2_fastq ";
	}else{
		$cmd .= " -notpaired ";
	}
	&process_cmd($cmd);
}

sub prelim_processing(){
	my $assembly_txt_file = shift;
	chdir($cwd);
	my $assembly_path = dirname($assembly_txt_file);
	chdir($assembly_path);
	print "Processing assembly: $assembly_txt_file at $assembly_path\n";
	my $data_hash_ref = &read_ncbi_genome_txt($assembly_txt_file);
	my $genome_fasta = $assembly_txt_file;
	$genome_fasta =~s/assembly_report.txt/genomic.fna/;
	my $genome_dbname = $data_hash_ref->{'Organism name'}.'_'.$data_hash_ref->{'Assembly name'};
	$genome_dbname=~s/\W+/_/g; 
	$genome_dbname=~s/__+/_/g;
	my $genome_gff = $genome_fasta;
	$genome_gff =~s/_genomic.fna$/_genomic.gff/;
	return ($assembly_path,$genome_fasta,$genome_dbname,$genome_gff,$data_hash_ref);
}

sub process_genome(){
	my ($fasta,$dbname,$gff_file,$species_name,$assembly_name) = @_;
	return unless $fasta && -s $fasta;
        my $file_date = &file_mod_month_year($fasta);

	if (!-s "$fasta.trim.x.n50"){
		&process_cmd("$trim_fasta_all_exec  -le 5000 -df -i $fasta")  unless -s "$fasta.trim";
		# convert to UC and remove IUPAC
		&process_cmd("$cleanup_iupac_exec $fasta.trim") unless -s "$fasta.trim.x";
		&process_cmd("$N50stats_exec $fasta.trim.x") ;
	}if (!-d "$gmap_db_dir/$dbname"){
		my $cmd = "$align_dnaseq_gsnap_exec -fasta $fasta.trim.x -dbname $dbname -gmap_dir $gmap_db_dir -suffix -build_only";
		$cmd .= " -large_genome " if -s $fasta > 3000000000;
		&process_cmd($cmd." &");
	}
	# samtools faidx already done in align_dnaseq
	# galaxy loc (1,2 of faidx)
		die "No .fai file found" unless -s "$fasta.trim.x.fai";
		&process_cmd("cut -f 1,2 $fasta.trim.x.fai > $fasta.trim.x.loc") unless -s "$fasta.trim.x.loc";

	# dbkey
    		# $dbname + friendly name + location of loc
	unless (-s "$fasta.trim.x.galaxy.dbkey.loc"){
		open (OUT,">$fasta.trim.x.galaxy.dbkey.loc");
                print OUT "$dbname\t$species_name ($assembly_name from $file_date)\t$fasta.trim.x.loc\n";
		close OUT;
	}

	unless (-s "$fasta.trim.x.galaxy.fasta.loc"){
		open (OUT,">$fasta.trim.x.galaxy.fasta.loc");
                print OUT "$dbname\t$dbname\t$species_name ($assembly_name from $file_date)\t$fasta.trim.x\n";
		close OUT;
	}
	my ($makeblastdb_exec,$bowtie_build_exec,$bwa_exec,$faToTwoBit_exec,$blat_exec,$bowtie2_build) =
		&check_programs('makeblastdb','bowtie-build','bwa','faToTwoBit','blat','bowtie2-build');

	# blast-ncbi
 		&process_cmd("sleep 10 && makeblastdb -dbtype nucl -parse_seqids -hash_index -in $fasta.trim.x &") unless -s "$fasta.trim.x.nin";
	# bowtie1
 		&process_cmd("sleep 10 && bowtie-build $fasta.trim.x $fasta.trim.x &") unless -s "$fasta.trim.x.1.ebwt";

	# bwa
		&process_cmd("sleep 10 && bwa index $fasta.trim.x") unless -s "$fasta.trim.x.amb" || -s "$fasta.trim.x.ann";

	# blat
 		&process_cmd("sleep 10 && faToTwoBit $fasta.trim.x $fasta.trim.x.2bit && blat $fasta.trim.x.2bit /dev/null /dev/null -makeOoc=$fasta.trim.x.11.ooc") unless -s "$fasta.trim.x.2bit"; 

	# bowtie2
 		&process_cmd("sleep 10 && bowtie2-build --threads 4 $fasta.trim.x $fasta.trim.x") unless -s "$fasta.trim.x.1.bt2";

	# annotation specific:
		# rnastar todo
	
	#should always have GFF but not always annotated
	if ($gff_file && -s $gff_file){
		my $annotation_check = `grep "mRNA" $gff_file`; # do we have real annotations?
		if ($annotation_check){
			&process_cmd("$gt_exec gff3 -sort -tidy -retainids $gff_file|uniq|grep -v '^#'|grep -vP '"
				.'\tcDNA_match\t|\tprotein_match\t|\ttRNA\t|\trRNA\t|\tpiRNA\t|\tlnc_RNA\t|\ttranscript\t'
				."' > $gff_file.tidy") unless -s "$gff_file.tidy";
			system("grep -P '"
				.'\ttRNA\t|\trRNA\t|\tpiRNA\t|\tlnc_RNA\t'
				."' $gff_file > $gff_file.tidy.ncRNA") unless -s "$gff_file.tidy.ncRNA";
			system("grep -P '"
				.'\ttranscript\t'
				."' $gff_file > $gff_file.tidy.transcriptRNA") unless -s "$gff_file.tidy.transcriptRNA";
			# we set the -lettername because sometimes the GenBank accession is not unique in the assembly GFF: the Name is used and it can actually link 
			# two scaffolds together
			my $tag = 'tag';
			if ($gff_file=~/^(G[\w_]+)/){$tag = $1;}
			&process_cmd("$create_features_exec -gff $gff_file.tidy -genome $fasta.trim.x -name -lettername -simple -genbank $tag 2> $gff_file.tidy.log ") unless -s "$gff_file.tidy.mRNA";
			if (!-s "$gff_file.tidy.mRNA"){
				#failed, ie no real annotation (often just tRNAs)
				system("rm -f $gff_file.tidy.*");
			}else{
		 		&process_cmd("makeblastdb -dbtype prot -parse_seqids -hash_index -in $gff_file.tidy.pep") unless -s "$gff_file.tidy.pep.pin";
				&process_cmd("$salmon_exec --no-version-check index -k 31 -p 6 -t $gff_file.tidy.mRNA -i $salmon_db_dir/$dbname") unless -d "$salmon_db_dir/$dbname";
				&process_cmd("$gff3_splicesites_exec < $gff_file.tidy.gff3 | $iit_store_exec -o $gmap_db_dir/$dbname/$dbname.maps/ncbi_gff.splicesites >/dev/null 2>/dev/null") unless -s "$gmap_db_dir/$dbname/$dbname.maps/ncbi_gff.splicesites.iit" || !-s "$gff_file.tidy.gff3";
			}
			# even if we have annotations, some genomes (eg fungi) have no introns (weird, huh? bad annotation?)
			if (-s "$gmap_db_dir/$dbname/$dbname.maps/ncbi_gff.splicesites.iit" && -s "$gmap_db_dir/$dbname/$dbname.maps/ncbi_gff.splicesites.iit" < 100){
				unlink("$gmap_db_dir/$dbname/$dbname.maps/ncbi_gff.splicesites.iit");
				system("echo Unannotated > $gff_file.unannotated");
			}
		}
	}


	# kmers for repeat estimation
	return if (-s 'genome_kmers_31_histo.csv');
	my $kmer_table_size = int(-s $fasta)+10;
	&process_cmd("$tessel_exec -k 31 -g $kmer_table_size -t 2 -tmp ".$ENV{'TMP'}." -canonical -f fasta "
		." genome_kmers $fasta.trim.x");
	open (IN,'genome_kmers_31_histo.txt');
	open (OUT,'>genome_kmers_31_histo.csv');
	while (my $ln=<IN>){
		next if $ln=~/^\s*$/ || $ln=~/[^0-9\s\.]/;
		print OUT $ln;
	}
	close OUT;
	close IN;
	unlink("genome_kmers_31.cbt");
}

sub download_sras(){
	my $accession_file = shift;
	open (IN,$accession_file);
	my $check = <IN>;
	close IN;
	return if $check =~/Unavailable/;
	open (IN,$accession_file);
	open (OUT,">$accession_file.t");
	while (my $ln=<IN>){
		chomp($ln);
		print OUT $ln."\n" unless ()=glob("./$ln*");
	}
	close OUT;
	close IN;
	return if !-s "$accession_file.t";
	print "\tDOWNLOADING SRA DATA\n";
	my $max_size = 100 * 1024 * 1024; #100 gb; def 20gb
	&process_web("$prefetch_exec --transport fasp --progress 1 --force no --max-size $max_size --option-file $accession_file.t");
}

sub process_sras(){
	my () = @_;
	
	print "\tPROCESSING FASTQ DATA\n";
	&process_cmd("$auto_sra_preprocess_reads_exec -naming sra -noqc -no_delete_raw") unless ()=glob("./*trimmomatic*");
	my @check = (glob("./*fastq"),glob("./*fastq.bz2"));
	if ( @check ){
		&process_cmd("pbzip2 -dp6k *.trimmomatic.bz2 2>/dev/null");
		my @trimmed = glob("./*trimmomatic");
		unless (@trimmed){
			warn "Processing of FASTQ failed\n";
			return 'FAIL';
		}
		#don't need the untrimmed
		&process_cmd("rm -f *unpaired* *fastqc*");
		foreach my $trim (@trimmed){
			$trim=~s/.trimmomatic$//;
			unlink($trim) if -s $trim;
			unlink($trim.'.bz2') if -s $trim.'.bz2';
		}
	}
}

sub convert_sras(){
	my @check1 = glob("./*trimmomatic*");
	my @check2 = glob("./?RR??????");
	return if scalar(@check1) == scalar(@check2) && scalar(@check1) >1;
	print "\tCONVERTING SRA DATA TO FASTQ\n";
	my $accession_file = shift;
	open (IN,$accession_file);
	my @sras = <IN>;
	chomp(@sras);
	close IN;
	foreach my $sra_accession (@sras){
		next if ()=glob("./$sra_accession*");
		my $sra_file = $ncbi_download_path . $sra_accession.'.sra';
		return 'FAIL' if !-s $sra_file;
		&process_cmd("$fastq_dump_exec --split-spot --split-files --skip-technical -F -Q 33 -W -T -R pass $sra_file");
	}
	&process_cmd("$rename_SRA_trimmed_fastq_exec");
	unless (()=glob("./*fastq*")){
		warn "We couldn't produce FASTQ files from this directory\n";
		return 'FAIL';
	}
	#removing DIRs
	foreach my $sra_accession (@sras){
		system("find $sra_accession -type d -empty -delete \; 2>/dev/null") if -d $sra_accession;
		my $sra_file = $ncbi_download_path . $sra_accession.'.sra';
		unlink($sra_file.'.sra') unless $no_cleanup;
	}
}


sub fetch_sra_accessions_biosample(){
	my ($biosample_id,$library_source,$library_strategy,$seq_platform) = @_;
	return if -s "$biosample_id.unavailable";
	unless (-f "$biosample_id.xml"){
		&process_web("$epost_exec -db biosample -format acc -id $biosample_id| $elink_exec -target sra | $efetch_exec -format xml > $biosample_id.xml");
	}
	&process_cmd("$xtract_exec -input $biosample_id.xml -pattern EXPERIMENT_PACKAGE -element PRIMARY_ID LIBRARY_SOURCE LIBRARY_STRATEGY INSTRUMENT_MODEL "
		. " > $biosample_id.text") if -s "$biosample_id.xml" && !-s "$biosample_id.text";
	return unless -s "$biosample_id.text";
	my $outfile = getcwd()."/$biosample_id.SRA.ids";
	return $outfile if -s $outfile;
	print "\tACQUIRING SRA ACCESSIONS for $biosample_id\n";
	open (IN, "$biosample_id.text");
	open (OUT,">$outfile");
	while (my $ln=<IN>){
		chomp($ln);
		next if !$ln || $ln=~/^\s*$/;
		my @data = split("\t",$ln);
		next unless $data[3]; #min 4 cols
		next if $library_source && $data[-3] !~/$library_source/i;
		next if $library_strategy && $data[-2] !~/$library_strategy/i;
		next if $seq_platform && $data[-1] !~/$seq_platform/i;
		while ($ln=~/\s([ES]RR\d+)\s/g){
			print OUT $1."\n";
		}
	}
	close OUT;
	close IN;
	if (!-s $outfile){
		system("echo UnavailableSRA > $biosample_id.unavailable");
		system("echo UnavailableSRA > $outfile");
		unlink($outfile);
		return;
	}
	return $outfile;
}

sub fetch_biosample_taxid(){
	my ($taxid) = @_;
	my @biosample_list;
	return if -s "$taxid.unavailable";
	my $outfile = "$taxid.biosample.ids";
	if (-s $outfile){
		open (IN,$outfile);
		@biosample_list = <IN>;
		chomp(@biosample_list);
		close IN;
		return @biosample_list;
	}
	print "\tACQUIRING BioSample ACCESSIONS for Taxonomy $taxid\n";
	unless (-s "$taxid.biosamples.xml"){
		&process_web("$epost_exec -db taxonomy -id $taxid | elink -target biosample|efetch -format xml > $taxid.biosamples.xml");
	}
	&process_cmd("$xtract_exec -input $taxid.biosamples.xml -pattern BioSample -first Id > $taxid.biosamples.text") if -s "$taxid.biosamples.xml" && !-s "$taxid.biosamples.text";
	return unless -s "$taxid.biosamples.text";
	open (IN,"$taxid.biosamples.text");
	open (OUT,">$outfile");
	while (my $ln = <IN>){
		chomp($ln);
		my @data = split("\t",$ln);
		next unless $data[0] && $data[0] =~/^SAM/;
		print OUT $data[0]."\n";
		push(@biosample_list,$data[0]);
	}
	close OUT;
	close IN;
	if (!-s $outfile){
		system("echo UnavailableBioSamples > ../$taxid.unavailable");
		system("echo UnavailableBioSamples > $outfile");
		unlink($outfile);
		chdir("../");
		rmdir($taxid);
		return;
	}
	return @biosample_list;
}


sub read_ncbi_genome_txt(){
	my $file = shift;
	print "PARSING GENOME ASSEMBLY INFO FILE $file\n";
	my %out_hash;
	open (IN, $file);
	while (my $line=<IN>){
		last if $line!~/^#/;
		chomp($line);
		my @data = split(":",$line);
		last unless $data[1];
		$data[1] =~s/^\s+//;
		$data[1] =~s/\s+$//;
		$data[0] =~s/\s+$//;
		$data[0] =~s/^#\s+//;
		$out_hash{$data[0]} = $data[1];

	}
	close IN;
	return \%out_hash;
}



sub get_metadata_from_taxid(){
	my $taxid = shift;
	return unless $taxid && $taxid=~/^\d+$/;
	my $html_file = "$taxid.html";
	system("wget -q 'https://www.ncbi.nlm.nih.gov/assembly/organism/browse-organism/$taxid/latest/?p\$l=Ajax&page=1&pageSize=99999&sortCol=1&sortDir=1' -O $html_file >/dev/null");
	my $data_hash_ref = &process_ncbi_genome_metadata($html_file,$taxid);
}

sub process_ncbi_genome_metadata(){
	my ($file, $taxid) = @_;
	my (@urls_to_fetch,%out_hash);
	open (OUT,">$file.out.tsv");

	my $orig_sep = $/;
	$/ = '<tr';
	open (IN,$file);
	my $total_data_line = <IN>;
	my ($total,$counter)=(int(0),int(0));
	if ($total_data_line=~/\s(\d+)\s-->/){
		$total = $1;
	}
	my $hs = HTML::Strip->new();
	while (my $record=<IN>){
		chomp($record);
		next unless $record;
		my $url;
		if ($record=~/(www.ncbi.nlm.nih.gov[^"]+)/){
			$url = 'https://'.$1;
		}
		next unless $url;
		$record = $hs->parse( $record );
		my @lines = split("\n",$record);
		my $skip = shift(@lines);
		my @data = ($taxid);
#TODO: cleanup columns for directories
		foreach my $ln (@lines){
			next if !$ln || $ln=~/^\s+$/;
			$ln=~s/^\s+//;
			$ln=~s/\s+$//;
			push(@data,$ln);
		}
		push(@data,$url);
		push(@data,&process_ncbi_genome_metadata_parse_download($url));
		#$data[-2] is the FTP URL
		push(@urls_to_fetch, $data[-2]);
		$counter++;
		print OUT join("\t",@data)."\n";
		print "Processed: $counter / $total       \r";
	}
	$hs->eof;
	close (IN);
	close (OUT);
	$/ = $orig_sep;
	return \%out_hash;
	&grab_ftp(@urls_to_fetch) if $fetch_ftps;
}

sub grab_ftp(){
	my @urls = @_;
	print "Fetching assembly metadata\n";
	system("wget -c -q --mirror --tries=10 --accept assembly_report.txt,genomic.fna.gz,genomic.gff.gz,rm.out.gz ".join(' ',@urls));

}


sub process_ncbi_genome_metadata_parse_download(){
	my $url = shift;
	my $document = `wget -q $url -O /dev/stdout 2>/dev/null`;
	my ($assembly_ftp,$n50);
	my @lines = split("\n",$document);
	for (my $i=0;$i<scalar(@lines);$i++){
#		last if $i < scalar(@lines);
		my $ln = $lines[$i];
		if (!$assembly_ftp && $ln=~/a href="(ftp:\/\/ftp.ncbi.nlm.nih.gov\/genomes\/all\/[^"]+")/){
			$assembly_ftp = dirname($1);
                }elsif ($ln=~/<td>Scaffold N50<\/td>\D+([\d,]+)<\/td>/){
			$n50 = $1;
		}elsif (!$n50 && $ln=~/<td>Contig N50<\/td>\D+([\d,]+)<\/td>/){
                        $n50 = $1;
                }
	}
	$n50=~s/,//g;
	return ($assembly_ftp,$n50);
}


sub mytime() {
 my @mabbr = qw(January February March April May June July August September October November December);
 my @wabbr = qw(Sunday Monday Tuesday Wednesday Thursday Friday Saturday);
 my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
 #my $sec   = localtime->sec() < 10 ? '0' . localtime->sec() : localtime->sec();
 #my $min   = localtime->min() < 10 ? '0' . localtime->min() : localtime->min();
 #my $hour = localtime->hour() < 10 ? '0' . localtime->hour() : localtime->hour();
 #my $wday = $wabbr[ localtime->wday() ];
 #my $mday = localtime->mday();
 #my $mon  = $mabbr[ localtime->mon ];
 #my $year = localtime->year() + 1900;

 $sec   = $sec < 10 ? '0' . $sec : $sec;
 $min   = $min < 10 ? '0' . $min : $min;
 $hour = $hour < 10 ? '0' . $hour : $hour;
 $wday = $wabbr[ $wday ];
 $mon  = $mabbr[ $mon ];
 $year = $year + 1900;
 return "$wday, $mon $mday, $year: $hour:$min:$sec\t";
}

sub process_cmd {
 my ( $cmd, $tries ) = @_;
 $tries = 3 if !$tries;
 print &mytime . "CMD: $cmd\n";
 my $ret = system($cmd);
 if ( $ret){
	my $my_tries = 1;
	# redo if NCBI query fails
	while (($ret == 65280 || $ret == 768) && $my_tries < $tries){
		sleep(20);
		$ret = system($cmd);
		$my_tries++;
	}
	if ($ret == 65280){
		warn "NCBI says no data for this item\n";
		return;
	}
 	if ( $ret != 256 ) {
	  confess "Error, cmd died with ret $ret\n";
	}
 }
 return;
}

sub process_web {
 my ( $cmd, $tries ) = @_;
 $tries = 10 if !$tries;
 print &mytime . "CMD: $cmd\n";
 my $ret = system($cmd);
 if ( $ret){
	my $my_tries = 1;
	# redo if NCBI query fails
	while ($ret && $ret != 256 && $my_tries < $tries){
		sleep(20);
		$ret = system($cmd);
		$my_tries++;
	}
 }
 return;
}

sub check_programs() {
  my @progs = @_;
  my @paths;
  foreach my $prog (@progs){
	  my $path = `which $prog`;	
	  unless ($path =~ /^\//){
	    $prog = basename($prog);
	    $path = `which $prog`;
	  }
	  confess "Error, path to $prog cannot be found"
	    unless $path =~ /^\//;
	  chomp($path);
	  $path = readlink($path) while -l $path;
	push (@paths,$path);
 }
 return @paths;
}


sub thousands(){
        my $val = shift;
        $val = sprintf("%.0f", $val);
        return $val if length($val)<4;
        1 while $val =~ s/(.*\d)(\d\d\d)/$1,$2/;
        return $val;
}

sub file_mod_month_year(){
 my $filename = shift;
 my $mtime = (stat ($filename)->[9]);
 my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($mtime);
 #my $mon = localtime -> mon($mtime);
 #my $year = localtime -> year($mtime);
 #my $mday = localtime -> mday($mtime);
 my @mabbr = qw(January February March April May June July August September October November December);
 $mon  = $mabbr[ $mon ] || die "Cannot get month from file $filename";
 $year = $year + 1900;
 return "$mday/$mon/$year";
}
