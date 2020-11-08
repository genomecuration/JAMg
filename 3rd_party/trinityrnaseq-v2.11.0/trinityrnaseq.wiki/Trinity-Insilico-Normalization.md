# Trinity's In silico Read Normalization 


Large RNA-Seq data sets, such as those exceeding 300M pairs, are best suited for in silico normalization prior to running Trinity.

>Note, all recent versions of Trinity will perform in silico normalization by default.  You can turn it off with Trinity --no_normalize_reads

If for any reason, you'd like to normalize data separately, you can do so as follows:

     %  $TRINITY_HOME/util/insilico_read_normalization.pl 

     ###############################################################################
     #
     # Required:
     #
     #  --seqType <string>      :type of reads: ( 'fq' or 'fa')
     #  --JM <string>            :(Jellyfish Memory) number of GB of system memory to use for 
     #                            k-mer counting by jellyfish  (eg. 10G) *include the 'G' char
     #                     
     #
     #  --max_cov <int>         :targeted maximum coverage for reads.
     #
     #
     #  If paired reads:
     #      --left  <string>    :left reads
     #      --right <string>    :right reads
     #
     #  Or, if unpaired reads:
     #      --single <string>   :single reads
     #
     #  Or, if you have read collections in different files you can use 'list' files, where each line in a list
     #  file is the full path to an input file.  This saves you the time of combining them just so you can pass
     #  a single file for each direction.
     #      --left_list  <string> :left reads, one file path per line
     #      --right_list <string> :right reads, one file path per line
     #
     ####################################
     ##  Misc:  #########################
     #
     #  --pairs_together                :process paired reads by averaging stats between pairs and retaining linking info.
     #
     #  --SS_lib_type <string>          :Strand-specific RNA-Seq read orientation.
     #                                   if paired: RF or FR,
     #                                   if single: F or R.   (dUTP method = RF)
     #                                   See web documentation.
     #  --output <string>               :name of directory for output (will be
     #                                   created if it doesn't already exist)
     #                                    
     #
     #  --CPU <int>                     :number of threads to use (default: = 2)
     #  --PARALLEL_STATS                :generate read stats in parallel for paired reads
     #
     #  --KMER_SIZE <int>               :default 25
     #
     #  --max_pct_stdev <int>           :maximum pct of mean for stdev of kmer coverage across read (default: 200)
     #
     #  --no_cleanup                    :leave intermediate files                      
     #  --tmp_dir_name <string>         default("tmp_normalized_reads");
     #
     ###############################################################################



This should be run on a machine that has a suitably large amount of RAM (typically hundreds of GB of RAM). 
The command-line options are quite similar to those used by Trinity itself.

There are several potential invocations of the normalization process, depending on your interests.  

Ideally, the process would be run as follows:

     % $TRINITY_HOME/util/insilico_read_normalization.pl --seqType fq \
               --JM 100G --max_cov 30 --left left.fq --right right.fq \
               --pairs_together --PARALLEL_STATS --CPU 10 

To roughly halve memory requirements but take almost twice as long to run, remove the '--PARALLEL_STATS' parameter.  In this case, the left.fq and right.fq files will be examined separately instead of in parallel.

## Results 

Normalized versions of the reads will be written, with file names corresponding to the original files but with an added extension to indicate the targeted maximum coverage levels.

##Sample data

You can run the normalization process on the sample data provided like so:

     % cd $TRINITY_HOME/sample_data/test_InSilicoReadNormalization/

     %  ./test_PE_normalization.sh

## References

The algorithm used for in silico read normalization is described in the http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3875132/bin/NIHMS537313-supplement-supplementary_text.pdf [supplement] to our http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3875132/ [Nature Protocol paper (Nat Protoc. Aug 2013; 8(8): 10.1038/nprot.2013.084)].

Our method is based on that earlier described by C. Titus Brown and implemented in http://ged.msu.edu/papers/2012-diginorm/ [diginorm].  Titus has a thoughtful [blog post](http://ivory.idyll.org/blog/trinity-in-silico-normalize.html) contrasting our approaches.
