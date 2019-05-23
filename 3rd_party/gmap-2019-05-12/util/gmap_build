#! /usr/bin/perl
# $Id: gmap_build.pl.in 215174 2018-05-12 01:22:21Z twu $

use warnings;	

my $gmapdb = "/mnt/nfs/home/30042108/software/JAMg/3rd_party/gmap-2019-05-12/share";
my $package_version = "2019-05-12";

use IO::File;
use File::Copy;	
#use File::Basename;
use Getopt::Long;


Getopt::Long::Configure(qw(no_auto_abbrev no_ignore_case_always));

# Default values
my $bindir = "/mnt/nfs/home/30042108/software/JAMg/3rd_party/gmap-2019-05-12/bin";   # dirname(__FILE__)
my $sampling = 3;
my $sleeptime = 2;

GetOptions(
    'B=s' => \$bindir,		# binary directory
    'T=s' => \$builddir,	# temporary build directory

    'D|dir=s' => \$destdir,	# destination directory
    'd|db=s' => \$dbname,	# genome name
    'n|names=s' => \$chrnamefile,   # substitute chromosomal names

    'M|mdfile=s' => \$mdfile,	# NCBI MD file
    'C|contigs-are-mapped' => \$contigs_mapped_p, # Each contig contains a chromosome tag in the FASTA header

    'z|compression=s' => \$compression_types, # compression types
    'j|local=s' => \$localsize, # k-mer size for local index (allowed: 8 or less)
    'k|kmer=s' => \$kmersize, # k-mer size for genomic index (allowed: 16 or less)
    'q=s' => \$sampling,	   # sampling interval for genome (default: 3)

    's|sort=s' => \$sorting,	# Sorting
    'g|gunzip' => \$gunzipp,	# gunzip files
    'E|fasta-pipe=s' => \$fasta_pipe,		  # Pipe for providing FASTA files
    'Q|fastq' => \$fastqp, # fastq files
    'R|revcomp' => \$revcompp, # reverse complement all reads

    'w=s' => \$sleeptime, # waits (sleeps) this many seconds between steps.  Useful if there is a delay in the filesystem.

    'c|circular=s' => \$circular,    # Circular chromosomes
    '2|altscaffold=s' => \$altscaffold,  # File with altscaffold info

    'e|nmessages=s' => \$nmessages,  # Max number of warnings or messages to print

    'p|part=s' => \$part	# Build in parts
    );


if (defined($builddir)) {
    print STDERR "Note: the -T flag is no longer necessary, since gmap_build now builds directly in the destination directory\n";
}

if (defined($compression_types)) {
    $compression_flag = "-z $compression_types";
} else {
    $compression_flag = "";
}

if (!defined($kmersize)) {
    print STDERR "-k flag not specified, so building with default 15-mers\n";
    $kmersize = 15;
}

if (!defined($localsize)) {
    print STDERR "-j flag not specified, so building with default 8-mers\n";
    $localsize = 8;
}


if (!defined($dbname)) {
    print_usage();
    die "Must specify genome database name with -d flag.";
} elsif ($dbname =~ /(\S+)\/(\S+)/) {
    $dbdir = $1;
    $dbname = $2;
    if (defined($destdir) && $destdir =~ /\S/) {
	# Note: The -D and -F arguments to gmapindex are different from the -D argument to gmap/gsnap.
	# For gmapindex, we use -D /path/to/dir/dbname -d dbname.  For gmap/gsnap, we use -D /path/to/dir -d dbname.
	$destdir = $destdir . "/" . $dbname;
    } else {
	$destdir = $dbdir;
    }
}

$dbname =~ s/\/$//;	# In case user gives -d argument with trailing slash

if (!defined($destdir) || $destdir !~ /\S/) {
    print STDERR "Destination directory not defined with -D flag, so writing to $gmapdb\n";
    $destdir = $gmapdb;
}

if (defined($sorting)) {
    if ($sorting ne "names") {
	$chr_order_flag = "-s $sorting";
    } elsif (!defined($chrnamefile)) {
	die "For --sort=names option, must also provide a filename to --names";
    } else {
	$chr_order_flag = "-s $sorting -n $chrnamefile";
    }
} else {
    # Default is to order genomes
    print STDERR "Sorting chromosomes in chrom order.  To turn off or sort other ways, use the -s flag.\n";
    $chr_order_flag = "";
}

if (!defined($gunzipp)) {
    $gunzip_flag = "";
} elsif (defined($fasta_pipe)) {
    die "Cannot use the -E (--fasta-pipe) flag with the -g (--gunzip) flag";
} else {
    $gunzip_flag = "-g";
}

if (!defined($fastqp)) {
    $fastq_flag = "";
} else {
    $fastq_flag = "-Q";
}

if (!defined($revcompp)) {
    $revcomp_flag = "";
} else {
    $revcomp_flag = "-R";
}

if (defined($circular)) {
    $circular_flag = "-c $circular";
} else {
    $circular_flag = "";
}

if (defined($altscaffold)) {
    $altscaffold_flag = "-2 $altscaffold";
} else {
    $altscaffold_flag = "";
}

if (defined($nmessages)) {
    $nmessages_flag = "-e $nmessages";
} else {
    $nmessages_flag = "";
}

if (defined($contigs_mapped_p)) {
    $contigs_mapped_flag = "-C";
} else {
    $contigs_mapped_flag = "";
}


#my @quoted = ();
#foreach $fasta (@ARGV) {
#    push @quoted,"\"$fasta\"";
#}
#my $genome_fasta = join(" ",@quoted);

$dbdir = create_db($destdir,$dbname);
$genomecompfile = "$dbdir/$dbname.genomecomp";

my $coordsfile = "$destdir/$dbname.coords";
my $fasta_sources = "$destdir/$dbname.sources";

$FP = new IO::File(">$fasta_sources") or die "Could not create $fasta_sources";
foreach $fasta (@ARGV) {
    print $FP $fasta . "\n";
}
close($FP);


#####################################################################################

check_compiler_assumptions();

if (!defined($part) || $part == 1) {

    create_genome_version($dbdir,$dbname);

    create_coords($mdfile,$fastq_flag,$fasta_pipe,$gunzip_flag,$circular_flag,$altscaffold_flag,$contigs_mapped_flag,$chrnamefile,
		  $bindir,$coordsfile,$fasta_sources);
    if (!(-s "$coordsfile")) {
	die "ERROR: $coordsfile not found";
    } else {
	$gmap_process_pipe = make_gmap_process_pipe($fastq_flag,$fasta_pipe,$gunzip_flag,$bindir,$coordsfile,$fasta_sources);
    }

    make_contig($nmessages_flag,$chr_order_flag,
		$bindir,$dbdir,$dbname,$gmap_process_pipe);

    compress_genome($nmessages_flag,$bindir,$dbdir,$dbname,$gmap_process_pipe);

    unshuffle_genome($bindir,$dbdir,$dbname,$genomecompfile);
}

if (!defined($part) || $part == 2) {
    $index_cmd = "\"$bindir/gmapindex\" -k $kmersize -q $sampling $nmessages_flag -d $dbname -F \"$dbdir\" -D \"$dbdir\"";

    if (($huge_genome_p = count_index_offsets($index_cmd,$genomecompfile)) == 1) {
	$index_cmd .= " -H";
    }

    create_index_offsets($index_cmd,$compression_flag,$genomecompfile);

    create_index_positions($index_cmd,$genomecompfile);
}

if (!defined($part) || $part == 3) {
    if (!defined($huge_genome_p)) {
	$huge_genome_p = count_index_offsets($index_cmd,$genomecompfile);
    }

    if ($huge_genome_p == 1) {
	print STDERR "Currently not creating localdb for large genomes\n";
    } else {
	$index_cmd = "\"$bindir/gmapindex\" -j $localsize $nmessages_flag -d $dbname -F \"$dbdir\" -D \"$dbdir\"";
	create_localdb($index_cmd,$compression_flag,$genomecompfile);
    }
}


#if (!defined($part) || $part == 4) {
#    if ($sarrayp == 1) {
#	make_enhanced_suffix_array($bindir,$dbdir,$dbname);
#    }
#}

if (!defined($part) || $part == 4) {
    system("rm -f \"$fasta_sources\"");
    system("rm -f \"$coordsfile\"");
}

exit;


#####################################################################################

sub check_compiler_assumptions {
    if (system("\"$bindir/gmapindex\" -9") != 0) {
	print STDERR "There is a mismatch between this computer system and the one where gmapindex was compiled.  Exiting.\n";
	exit(9);
    }
}

sub create_db {
    my ($destdir, $dbname) = @_;

    print STDERR "Creating files in directory $destdir/$dbname\n";
    system("mkdir -p \"$destdir\"");
    system("mkdir -p \"$destdir/$dbname\"");
    system("mkdir -p \"$destdir/$dbname/$dbname.maps\"");
    system("chmod 755 \"$destdir/$dbname/$dbname.maps\"");

    return "$destdir/$dbname";
}


sub create_genome_version {
    my ($dbdir, $dbname) = @_;

    open GENOMEVERSIONFILE, ">$dbdir/$dbname.version" or die $!;
    print GENOMEVERSIONFILE "$dbname\n";
    close GENOMEVERSIONFILE or die $!;
    sleep($sleeptime);
    return;
}

sub create_coords {
    my ($mdfile, $fastq_flag, $fasta_pipe, $gunzip_flag, $circular_flag, $altscaffold_flag, $contigs_mapped_flag, $chrnamefile,
	$bindir, $coordsfile, $fasta_sources) = @_;
    my ($cmd, $rc);

    if (defined($mdfile)) {
	# MD file cannot specify that a chromosome is circular or altscaffold
	$cmd = "\"$bindir/md_coords\" -o \"$coordsfile\" $mdfile";
    } else {
	if (defined($fasta_pipe)) {
	    $cmd = "$fasta_pipe | \"$bindir/fa_coords\" $revcomp_flag $fastq_flag $circular_flag $altscaffold_flag $contigs_mapped_flag -o \"$coordsfile\"";
	} else {
	    $cmd = "\"$bindir/fa_coords\" $gunzip_flag $revcomp_flag $fastq_flag $circular_flag $altscaffold_flag $contigs_mapped_flag -o \"$coordsfile\"";
	}
	if (defined($chrnamefile)) {
	    $cmd .= " -n $chrnamefile";
	}
	if (!defined($fasta_pipe)) {
	    $cmd .= " -f \"$fasta_sources\"";
	}
    }
    print STDERR "Running $cmd\n";
    if (($rc = system($cmd)) != 0) {
	die "$cmd failed with return code $rc";
    }
    sleep($sleeptime);
    return;
}

sub make_gmap_process_pipe {
    my ($fastq_flag, $fasta_pipe, $gunzip_flag, $bindir, $coordsfile, $fasta_sources) = @_;

    if (defined($fasta_pipe)) {
	return "$fasta_pipe | \"$bindir/gmap_process\" $fastq_flag -c \"$coordsfile\"";
    } else {
	return "\"$bindir/gmap_process\" $fastq_flag $gunzip_flag -c \"$coordsfile\" -f \"$fasta_sources\"";
    }
}

sub make_contig {
    my ($nmessages_flag, $chr_order_flag,
	$bindir, $dbdir, $dbname, $gmap_process_pipe) = @_;
    my ($cmd, $rc);

    $cmd = "$gmap_process_pipe | \"$bindir/gmapindex\" $nmessages_flag -d $dbname -D \"$dbdir\" -A $chr_order_flag";
    print STDERR "Running $cmd\n";
    if (($rc = system($cmd)) != 0) {
	die "$cmd failed with return code $rc";
    }
    sleep($sleeptime);
    return;
}

sub compress_genome {
    my ($nmessages_flag, $bindir, $dbdir, $dbname, $gmap_process_pipe) = @_;
    my ($cmd, $rc);

    $cmd = "$gmap_process_pipe | \"$bindir/gmapindex\" $nmessages_flag -d $dbname -F \"$dbdir\" -D \"$dbdir\" -G";
    print STDERR "Running $cmd\n";
    if (($rc = system($cmd)) != 0) {
	die "$cmd failed with return code $rc";
    }
    sleep($sleeptime);
    return;
}

sub unshuffle_genome {
    my ($bindir, $dbdir, $dbname, $genomecompfile) = @_;
    my ($cmd, $rc);

    $cmd = "cat \"$genomecompfile\" | \"$bindir/gmapindex\" -d $dbname -U > \"$dbdir/$dbname.genomebits128\"";
    print STDERR "Running $cmd\n";
    if (($rc = system($cmd)) != 0) {
	die "$cmd failed with return code $rc";
    }
    sleep($sleeptime);
    return;
}

# No longer supported
#sub full_ASCII_genome {
#    make_contig();
#	
#    $cmd = "\"$bindir/gmap_process\" $gunzip_flag -c $dbdir/$dbname.coords $genome_fasta | \"$bindir/gmapindex\" $nmessages_flag -d $dbname -F $dbdir -D $dbdir -l -G";
#    print STDERR "Running $cmd\n";
#    if (($rc = system($cmd)) != 0) {
#	die "$cmd failed with return code $rc";
#    }
#    sleep($sleeptime);
#    return;
#}

sub count_index_offsets {
    my ($index_cmd, $genomecompfile) = @_;
    my $huge_genome_p;
    my ($cmd, $noffsets);
    
    $cmd = "cat \"$genomecompfile\" | $index_cmd -N";
    print STDERR "Running $cmd\n";
    $noffsets = `$cmd`;
    chop $noffsets;
    if ($noffsets <= 4294967295) {
	print STDERR "Number of offsets: $noffsets => pages file not required\n";
	$huge_genome_p = 0;
    } else {
	print STDERR "Number of offsets: $noffsets => pages file required\n";
	$huge_genome_p = 1;
    }
    sleep($sleeptime);
    return $huge_genome_p;
}

sub create_index_offsets {
    my ($index_cmd, $compression_flag, $genomecompfile) = @_;
    my ($cmd, $rc);

    $cmd = "$index_cmd -O $compression_flag \"$genomecompfile\"";
    print STDERR "Running $cmd\n";
    if (($rc = system($cmd)) != 0) {
	die "$cmd failed with return code $rc";
    }
    sleep($sleeptime);
    return;
}

sub create_index_positions {
    my ($index_cmd, $genomecompfile) = @_;
    my ($cmd, $rc);

    $cmd = "$index_cmd -P \"$genomecompfile\"";
    print STDERR "Running $cmd\n";
    if (($rc = system($cmd)) != 0) {
	die "$cmd failed with return code $rc";
    }
    sleep($sleeptime);
    return;
}

sub create_localdb {
    my ($index_cmd, $compression_flag, $genomecompfile) = @_;
    my ($cmd, $rc);

    $cmd = "$index_cmd -Q $compression_flag \"$genomecompfile\"";
    print STDERR "Running $cmd\n";
    if (($rc = system($cmd)) != 0) {
	die "$cmd failed with return code $rc";
    }
    sleep($sleeptime);
    return;
}

# No longer supported
sub make_enhanced_suffix_array {
    my ($bindir, $dbdir, $dbname) = @_;
    my ($cmd, $rc);

    # Suffix array
    $cmd = "\"$bindir/gmapindex\" -d $dbname -F \"$dbdir\" -D \"$dbdir\" -S";
    print STDERR "Running $cmd\n";
    if (($rc = system($cmd)) != 0) {
	die "$cmd failed with return code $rc";
    }
    sleep($sleeptime);

    # LCP and child arrays
    $cmd = "\"$bindir/gmapindex\" -d $dbname -F \"$dbdir\" -D \"$dbdir\" -L";
    print STDERR "Running $cmd\n";
    if (($rc = system($cmd)) != 0) {
	die "$cmd failed with return code $rc";
    }
    sleep($sleeptime);

    # Compressed suffix array
    # $cmd = "\"$bindir/gmapindex\" -d $dbname -F \"$dbdir\" -D \"$dbdir\" -C";
    # print STDERR "Running $cmd\n";
    # if (($rc = system($cmd)) != 0) {
    # die "$cmd failed with return code $rc";
    # }
    # sleep($sleeptime);

    return;
}


sub print_usage {
  print <<TEXT1;

gmap_build: Builds a gmap database for a genome to be used by GMAP or GSNAP.
Part of GMAP package, version $package_version.

A simplified alternative to using the program gmap_setup, which creates a Makefile.

Usage: gmap_build [options...] -d <genomename> <fasta_files>

Options:
    -D, --dir=STRING          Destination directory for installation (defaults to gmapdb directory specified at configure time)
    -d, --db=STRING           Genome name

    -n, --names=STRING        Substitute names for chromosomes, provided in a file.  The file should have one line
                                for each chromosome name to be changed, with the original FASTA name in column 1 and
                                the desired chromosome name in column 2.  This provides an easy way to change the
                                names of chromosomes, for example, to add or remove the "chr" prefix.  Column 2 may
                                be blank, which indicates no name change.  This file can also be combined with
                                --sort=names to provide a particular order for the chromosomes in the genome index.

    -M, --mdflag=STRING       Use MD file from NCBI for mapping contigs to chromosomal coordinates
    -C, --contigs-are-mapped  Find a chromosomal region in each FASTA header line.  Useful for contigs that have been mapped
                                to chromosomal coordinates.  Ignored if the --mdflag is provided.

    -k, --kmer=INT            k-mer value for genomic index (allowed: 15 or less, default is 15)
    -q INT                    sampling interval for genomoe (allowed: 1-3, default 3)

    -s, --sort=STRING         Sort chromosomes using given method:
                                none - use chromosomes as found in FASTA file(s)
                                alpha - sort chromosomes alphabetically (chr10 before chr 1)
                                numeric-alpha - chr1, chr1U, chr2, chrM, chrU, chrX, chrY
                                chrom - chr1, chr2, chrM, chrX, chrY, chr1U, chrU
                                names - sort chromosomes based on file provided to --names flag

    -g, --gunzip              Files are gzipped, so need to gunzip each file first
    -E, --fasta-pipe=STRING   Interpret argument as a command, instead of a list of FASTA files
    -Q, --fastq               Files are in FASTQ format
    -R, --revcomp             Reverse complement all contigs
    -w INT                    Wait (sleep) this many seconds after each step (default 2)

    -c, --circular=STRING     Circular chromosomes (either a list of chromosomes separated by a comma, or
                                a filename containing circular chromosomes, one per line).  If you use the
                                --names feature, then you should use the original name of the chromosome,
                                not the substitute name, for this option.                                                        

    -2, --altscaffold=STRING  File with alt scaffold info, listing alternate scaffolds, one per line, tab-delimited,
                                with the following fields: (1) alt_scaf_acc, (2) parent_name, (3) orientation,
                                (4) alt_scaf_start, (5) alt_scaf_stop, (6) parent_start, (7) parent_end.

    -e, --nmessages=INT       Maximum number of messages (warnings, contig reports) to report (default 50)

TEXT1
  return;
}

