Alexie's Illumina preprocessing pipeline; no warranty but we use it...
complaints to alexie@butterflybase.org

Identifies rRNA, contaminants, trims in a standard fashion etc. Maintains read pairs. Built for fire-and-forget high throughput projects (terabytes of data).


# INSTALL

Uses pbzip2

On Ubuntu you can install most of these software:
$ sudo apt-get install pbzip2

Uses some perl libraries so:

$ sudo cpan Data::Dumper threads Getopt::Long Pod::Usage Digest::MD5 BSD::Resource FindBin

For Blue we will need mono:

$ sudo apt-get install mono-runtime mono-mcs

Finally to compile 3rd party software, run make

$ make

# USAGE

$ preprocess_illumina.pl <infile> <infile2> etc 

or for pairs:
$ preprocess_illumina.pl -paired <infile1> <infile2>

see
$ perldoc preprocess_illumina.pl
or code within
$ vim preprocess_illumina.pl


+++++++++++++++++++++++++
Programs used:

Trimmomatic
 http://www.usadellab.org/cms/index.php?page=trimmomatic
 Lohse M, Bolger AM, Nagel A, Fernie AR, Lunn JE, Stitt M, Usadel B. RobiNA: a user-friendly, integrated software solution for RNA-Seq-based transcriptomics. Nucleic Acids Res. 2012 Jul;40(Web Server issue):W622-7.

pbzip2
 Jeff Gilchrist
 http://compression.ca
 Debian/Ubuntu:
 apt-get install pbzip2

fastqc
 simon.andrews@babraham.ac.uk
 http://www.bioinformatics.bbsrc.ac.uk

allpathslg
 http://www.broadinstitute.org/software/allpaths-lg/blog/

blue
 http://www.bioinformatics.csiro.au/blue/
