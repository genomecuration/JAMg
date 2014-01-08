INCLUDED FILES:

 validate_gtf.pl	The GTF validator.

 eval.pl		The Eval GUI.
 evaluate_gtf.pl	The command line interface to the Evaluate function.
 get_general_stats.pl  	The command line interface to the GenStats function.
 filter_gtfs.pl		The command line interface to the Filter function.
 graph_gtfs.pl		The command line interface to the Graph function.
 get_overlap_stats.pl	The command line interface to the Overlap function.
 get_distribution.pl	The command line interface to the Dist function.

 GTF.pm			A Perl library for dealing with GTF files.
 Eval.pm		A Perl library containing all the Eval functions.

 eval-documentation.pdf	The complete Eval documentation including a
                        user guide and code documentation.
 eval-summary.pdf	A summary of the Eval system.

 Also included are three gtf files containing annotation of the NCBI
  build 31 / November 2002 Golden Path human chromosome 22.  A fasta
  file of the genomic sequence is available for download from
  genome.ucsc.edu.

 chr22.twinscan.gtf	
 chr22.genscan.gtf
 chr22.refseq.gtf
 

REQUIREMENTS:

The Eval system requires Perl 5.0 or later and Perl::Tk 8.0 or later.
The Perl::Tk module can be found at www.cpan.org.  To display graphs
from the GUI the gnuplot utility is required.  gnplot can be found at
www.gnuplot.info.


INSTALLATION:
	
To install the Eval package just unpack the files to any directory.
The library files (*.pm files) should be placed in a directory in the
users Perl library path.  For example, if the *.pm files are located
in /usr/lib/eval, then the Perl library path should include /usr/lib/eval:

PERL5LIB=/usr/lib:/usr/lib/eval

Alternatively the files may all be placed in a single directory, which
will allow the programs to be run only from that directory.

Correct installation can be tested using the three provided GTF files.
The command:

eval.pl -v chr22.refseq.gtf chr22.twinscan.gtf chr22.genscan.gtf 

should load the three files in the GUI without any errors if the
libraries have been installed correctly.

USE:

See the documentation for a detailed description of how to use the
programs.  
