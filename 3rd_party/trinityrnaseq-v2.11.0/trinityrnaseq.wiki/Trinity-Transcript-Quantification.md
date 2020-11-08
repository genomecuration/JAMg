# Trinity Transcript Quantification

There are now several methods available for estimating transcript abundance in a genome-free manner, and these include alignment-based methods (aligning reads to the transcript assembly) and alignment-free methods (typically examining k-mer abundances in the reads and in the resulting assemblies).

In Trinity, we provide direct support for running the alignment-based quantification methods [RSEM](http://deweylab.biostat.wisc.edu/rsem/) and [eXpress](http://bio.math.berkeley.edu/eXpress/), as well as the ultra-fast alignment-free method [kallisto](http://pachterlab.github.io/kallisto/) and 'wicked-fast' [salmon](http://salmon.readthedocs.org/en/latest/salmon.html).

The Trinity software does not come pre-packaged with any of these software tools, so be sure to download and install any that you wish to use. The tools should be available via your PATH setting (so, typing 'which kallisto' or 'which express' on the linux command line returns the path to where the tool is installed on your system).

>If you have multiple RNA-Seq data sets that you want to compare (eg. different tissues sampled from a single organism), be sure to **generate a single Trinity assembly** and to then run the abundance estimation separately for each of your samples.

# Table of Contents

* [Estimating transcript abundance](#estimating-transcript-abundance)
  * [Alignment based abundance estimation methods](#alignment-based-abundance-estimation-methods)
    * [RSEM output](#rsem-output)
    * [eXpress output](#express-output)
  * [Alignment-free abundance estimation methods](#alignment-free-abundance-estimation)
    * [kallisto output](#kallisto-output)
    * [salmon output](#salmon-output)
* [Building transcript and gene expression matrices](#building-expression-matrices)
* [Counting expressed transcripts and genes](#counting-expressed-genes)
* [Filtering Transcripts Based on Expression Values](#filtering-transcripts)

<a name='#estimating-transcript-abundance'></a>
## Estimating Transcript Abundance 

The Trinity toolkit comes with a script to facilitate running your choice of the above tools to quantitate transcript abundance:

     % $TRINITY_HOME/util/align_and_estimate_abundance.pl 

    #########################################################################
    #
    #  --transcripts <string>           transcript fasta file
    #  --seqType <string>               fq|fa
    # 
    #  If Paired-end:
    #
    #     --left <string>
    #     --right <string>
    #  
    #   or Single-end:
    #
    #      --single <string>
    #   or
    #      --samples_file <string>    tab-delimited text file indicating biological replicate relationships.
    #                                   ex.
    #                                        cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
    #                                        cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
    #                                        cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
    #                                        cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
    #
    #                      # if single-end instead of paired-end, then leave the 4th column above empty. 
    #
    #
    #
    #  --est_method <string>           abundance estimation method.
    #                                        alignment_based:  RSEM|eXpress       
    #                                        alignment_free: kallisto|salmon
    #  
    # --output_dir <string>            write all files to output directory
    #  
    #
    #  if alignment_based est_method:
    #       --aln_method <string>            bowtie|bowtie2|(path to bam file) alignment method.  (note: RSEM requires bowtie)
    #                                       (if you already have a bam file, you can use it here instead of rerunning bowtie)
    #
    # Optional:
    #  
    # --SS_lib_type <string>           strand-specific library type:  paired('RF' or 'FR'), single('F' or 'R').
    #                                         (note, no strand-specific mode for kallisto)
    #
    # --thread_count                   number of threads to use (default = 4)
    #
    # --debug                          retain intermediate files
    #
    #  --gene_trans_map <string>        file containing 'gene(tab)transcript' identifiers per line.
    #     or  
    #  --trinity_mode                   Setting --trinity_mode will automatically generate the gene_trans_map and use it.
    #
    #
    #  --prep_reference                 prep reference (builds target index)
    #
    #
    ########################################
    #
    #  Parameters for single-end reads:
    #
    #  --fragment_length <int>         specify RNA-Seq fragment length (default: 200) 
    #  --fragment_std <int>            fragment length standard deviation (defalt: 80)
    #
    ########################################
    #  
    #   bowtie-related parameters: (note, tool-specific settings are further below)
    #
    #  --max_ins_size <int>             maximum insert size (bowtie -X parameter, default: 800)
    #  --coordsort_bam                  provide coord-sorted bam in addition to the default (unsorted) bam.
    #
    ########################################
    #  RSEM opts:
    #
    #  --bowtie_RSEM <string>          if using 'bowtie', default: "--all --best --strata -m 300 --chunkmbs 512"
    #  --bowtie2_RSEM <string>         if using 'bowtie2', default: "--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 "
    #  --include_rsem_bam              provide the RSEM enhanced bam file including posterior probabilities of read assignments.
    #  --rsem_add_opts <string>        additional parameters to pass on to rsem-calculate-expression
    #
    ##########################################################################
    #  eXpress opts:
    #
    #  --bowtie_eXpress <string>  default: "--all --best --strata -m 300 --chunkmbs 512"
    #  --bowtie2_eXpress <string> default: "--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200 "
    #  --eXpress_add_opts <string>  default: ""
    #
    ##########################################################################
    #  kallisto opts:
    #
    #  --kallisto_add_opts <string>  default:   
    #
    ##########################################################################
    #
    #  salmon opts:
    #
    #  --salmon_idx_type <string>    quasi|fmd (defalt: quasi)
    #  --salmon_add_opts <string>    default: 
    #
    #
    #  Example usage
    #
    #   ## Just prepare the reference for alignment and abundance estimation
    #
    #    /home/unix/bhaas/GITHUB/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference
    #
    #   ## Run the alignment and abundance estimation (assumes reference has already been prepped, errors-out if prepped reference not located.)
    #
    #    /home/unix/bhaas/GITHUB/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left reads_1.fq --right reads_2.fq --est_method RSEM --aln_method bowtie --trinity_mode --output_dir rsem_outdir
    #
    ##  ## prep the reference and run the alignment/estimation
    #
    #    /home/unix/bhaas/GITHUB/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --seqType fq --left reads_1.fq --right reads_2.fq --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference --output_dir rsem_outdir
    #
    #########################################################################

If you have strand-specific data, be sure to include the '--SS_lib_type' parameter.

Before running the above, please consider the following:

>Please use the --samples_file parameter with the abundance estimation utility.  This will organize your outputs so that each replicate will be organized in its own output directory named according to the corresponding replicate.

>It is useful to first run 'align_and_estimate_abundance.pl' to only prep your reference database for alignment, using '--prep_reference', and then subsequently running it on each of your sets of reads in parallel to obtain sample-specific abundance estimates.

>If you quality-trimmed your reads using the --trimmomatic parameter in Trinity, you should consider using the corresponding quality-trimmed reads for the abundance estimation process outlined here. You'll find the quality-trimmed reads in the trinity_out_dir/ with a 'P.qtrim.gz' extension.



Any of the abundance estimation methods will provide transcript-level estimates of the count of RNA-Seq fragments that were derived from each transcript, in addition to a normalized measure of transcript expression that takes into account the transcript length, the number of reads mapped to the transcript, and the the total number of reads that mapped to any transcript. Normalized expression metrics may be reported as 'fragments per kilobase transcript length per million fragments mapped' (FPKM) or 'transcripts per million transcripts' (TPM).  The TPM metric is generally preferred to FPKM, given the property that all values will always sum up to 1 million (FPKM values will tend to not sum up to the same value across samples). For info on TPM vs. FPKM, see [Wagner et al. 2012. Theory Biosci](http://lynchlab.uchicago.edu/publications/Wagner,%20Kin,%20and%20Lynch%20(2012).pdf), and [Li and Dewey, BMC Bioinf. 2011](http://www.biomedcentral.com/1471-2105/12/323).


>Be sure to use the --gene_trans_map or --trinity_mode parameters in order to get a *gene counts matrix* in addition to the isoform counts matrix.

The estimated fragment counts are generally needed by many differential expression analysis tools that use count-based statistical models, and the normalized expression values (FPKM or TPM) are used almost everywhere else, such as plotting in heatmaps.

<a name='alignment-based-abundance-estimation-methods'></a>
## Alignment-based abundance estimation methods 
The alignment step generates the file 'bowtie.bam', which is then fed directly into either RSEM or eXpress.  Note, if parameter '--coordsort_bam ' is set, the process also generates a 'bowtie.csorted.bam' file, which is a coordinate-sorted bam file that can be used for visualization using IGV.


<a name='rsem-output'></a>
### RSEM output 
The RSEM computation generates two primary output files containing the abundance estimation information:

  RSEM.isoforms.results  : EM read counts per Trinity transcript
  RSEM.genes.results     : EM read counts on a per-Trinity-gene, 'gene' used loosely here.


The output for the isoforms file looks like so:

|transcript_id|   gene_id| length|  effective_length|        expected_count|  TPM|     FPKM|    IsoPct|
|-------------:|---------:|-----:|-----------------:|---------------------:|-----:|--------:|--------:|
|TRINITY_DN100_c0_g1_i1|  TRINITY_DN100_c0_g1|     443|     181.37|  4.84|    1670.06| 12311.85|        100.00|
|TRINITY_DN101_c0_g1_i1|  TRINITY_DN101_c0_g1|     251|     19.37|   1.00|    3231.22| 23820.87|        100.00|
|TRINITY_DN103_c0_g1_i1|  TRINITY_DN103_c0_g1|     1219|    957.37|  0.00|    0.00|    0.00|    100.00|
|TRINITY_DN103_c0_g2_i1|  TRINITY_DN103_c0_g2|     414|     152.41|  0.00|    0.00|    0.00|    100.00|
|TRINITY_DN104_c0_g1_i1|  TRINITY_DN104_c0_g1|     408|     146.44|  0.00|    0.00|    0.00|    0.00|
|TRINITY_DN106_c0_g1_i1|  TRINITY_DN106_c0_g1|     1111|    849.37|  1.00|    73.70|   543.30|  100.00|
|TRINITY_DN106_c1_g1_i1|  TRINITY_DN106_c1_g1|     339|     81.68|   0.00|    0.00|    0.00|    0.00|


<a name='express-output'></a>
### eXpress output 

The eXpress runner also generates two files:

     results.xprs  : transcript abundance estimates (generated by eXpress)
     results.xprs.genes : gene abundance estimates (provided here by summing up transcript values per gene)

Each includes the estimated counts, FPKM, and TPM measures, in addition to a large number of other metrics - see the [eXpress documentation](http://bio.math.berkeley.edu/eXpress/manual.html) for details.

<a name='alignment-free-abundance-estimation'></a>
## Alignment-free abundance estimation methods

<a name='kallisto-output'></a>
### kallisto 

The kallisto runner similarly generates two files:

     abundance.tsv  : transcript abundance estimates (generated by kallisto)
     abundance.tsv.genes : gene abundance estimates (provided here by summing up transcript values per gene)

The format of the output is short and sweet, providing the key essential metrics:

|target_id|       length|  eff_length|      est_counts|      tpm|
|:---------|-----------:|-----------:|---------------:|--------:|
|TRINITY_DN10_c0_g1_i1|   334|     100.489| 13|      4186.62|
|TRINITY_DN11_c0_g1_i1|   319|     87.9968| 0|       0|
|TRINITY_DN12_c0_g1_i1|   244|     38.2208| 2|       1693.43|
|TRINITY_DN17_c0_g1_i1|   229|     30.2382| 5|       5351.21|
|TRINITY_DN18_c0_g1_i1|   633|     384.493| 19|      1599.2|
|TRINITY_DN18_c1_g1_i1|   289|     65.795|  1|       491.864|
|TRINITY_DN19_c0_g1_i1|   283|     61.0618| 10|      5299.91|


<a name='salmon-output'></a>
### Salmon 

Likewise, after running salmon you'll find output files:

    quant.sf : transcript abundance estimates (generated by salmon)
    quant.sf.genes : gene-level abundance estimates (generated here by summing transcript values)

and the output format is also short and sweet:

|Name|    Length|  EffectiveLength| TPM|     NumReads|
|:---|---------:|----------------:|---:|------------:|
|TRINITY_DN10_c0_g1_i1|   334|     135.084| 0|       0|
|TRINITY_DN11_c0_g1_i1|   319|     120.926| 6158.51| 13|
|TRINITY_DN12_c0_g1_i1|   244|     57.9931| 0|       0|
|TRINITY_DN17_c0_g1_i1|   229|     47.8216| 2395.84| 2|
|TRINITY_DN18_c0_g1_i1|   633|     432.567| 662.169| 5|
|TRINITY_DN18_c1_g1_i1|   289|     93.8233| 3663.47| 6|
|TRINITY_DN19_c0_g1_i1|   283|     88.6594| 646.141| 1|
|TRINITY_DN21_c0_g1_i1|   242|     56.5791| 1012.5|  1|

<a name='building-expression-matrices'></a>
## Build Transcript and Gene Expression Matrices

Using the transcript and gene-level abundance estimates for each of your samples, construct a matrix of counts and a matrix of normalized expression values using the following script:

    $TRINITY_HOME/util/abundance_estimates_to_matrix.pl 

    ####################################################################################
    #
    # Usage:  trinityrnaseq/util/abundance_estimates_to_matrix.pl --est_method <method>  sample1.results sample2.results ...
    #
    #      or  trinityrnaseq/util/abundance_estimates_to_matrix.pl --est_method <method> --quant_files file.listing_target_files.txt
    #
    #      Note, if only a single input file is given, it's expected to contain the paths to all the target abundance estimation files.
    #
    # Required:
    #            
    #  --est_method <string>           RSEM|eXpress|kallisto|salmon  (needs to know what format to expect)
    #
    #  --gene_trans_map <string>           the gene-to-transcript mapping file. (if you don't want gene estimates, indicate 'none'.
    #
    #
    # Options:
    #
    #  --cross_sample_norm <string>         TMM|UpperQuartile|none   (default: TMM)
    #
    #  --name_sample_by_basedir             name sample column by dirname instead of filename
    #      --basedir_index <int>            default(-2)
    #
    #  --out_prefix <string>                default: value for --est_method
    #
    #  --quant_files <string>              file containing a list of all the target files.
    #
    ######################################################################################




For example, suppose you have two samples (sampleA and sampleB), and you ran kallisto to estimate transcript abundances, you might generate matrices like so:

    $TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method kallisto \
        --gene_trans_map Trinity.fasta.gene_trans_map \
        --out_prefix kallisto \
        --name_sample_by_basedir \
         sampleA/abundance.tsv \
         sampleB/abundance.tsv 

which would generate the following three files:

      kallisto.isoform.counts.matrix  : the estimated RNA-Seq fragment counts (raw counts)
      kallisto.isoform.TPM.not_cross_norm  : a matrix of TPM expression values (not cross-sample normalized)
      kallisto.isoform.TMM.EXPR.matrix : a matrix of TMM-normalized expression values

>I recommend using the --name_sample_by_basedir parameter here, which requires that you have all your expression results organized into separate directories, each named according to the corresponding replicate name.   If you ran the abundance estimation script above using the --samples_file parameter, it will automatically organize the data accordingly.

The 'counts.matrix' file is used for downstream analyses of differential expression.  The TMM.EXPR.matrix file is used as the gene expression matrix in most other analyses.  For information on the importance of TMM (or cross-sample normalization in general), see [Robinson & Oshlack, Genome Biology 2010](http://www.genomebiology.com/2010/11/3/R25) and [Dillies et al., Brief Bioinf, 2012](http://bib.oxfordjournals.org/content/14/6/671.long).

>When you include the --gene_trans_map file above, it will automatically generate the gene-level count and expression matrices, using the 'scaledTPM' method as described in [txImport](http://bioconductor.org/packages/release/bioc/html/tximport.html) but implemented here directly in the Trinity script.  This 'scaledTPM' method for estimating gene counts accounts for differences in isoform lengths that could otherwise lead to false gene DE reporting under situations where it is differential transcript usage (DTU) as opposed to differential gene expression (DGE) occurring. See [Soneson et al., F1000 Research, 2016](https://f1000research.com/articles/4-1521/v2) for details.


<a name='counting-expressed-genes'></a>
## Counting Numbers of Expressed Transcripts or Genes

Presumably, a transcript is expressed if it has been assembled from RNA-Seq data, but as we know, transcription can be quite pervasive, and many transcripts, particularly the very lowly expressed ones, have questionable biological significance.  Note that some transcripts may have artificially low (or zero) expression values simply because they are incompletely assembled and do not recruit both pairs of PE reads in order to be properly accounted for during abundance estimation.  If we assume that most biologically relevant transcripts are reasonably well assembled and well quantified by the abundance estimation method used, we might infer the approximate number of expressed genes or transcripts as the number that are expressed above some minimum expression threshold.

Given a matrix of TPM values (ideally, in this case, **not** the TMM normalized version), you can plot the number of genes (or transcripts) that are expressed above a minimum TPM expression threshold in any sample like so.


    $TRINITY_HOME/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
              genes_matrix.TPM.not_cross_norm | tee genes_matrix.TPM.not_cross_norm.counts_by_min_TPM

    and

     $TRINITY_HOME/util/misc/count_matrix_features_given_MIN_TPM_threshold.pl \
              trans_matrix.TPM.not_cross_norm | tee trans_matrix.TPM.not_cross_norm.counts_by_min_TPM


and, looking at the output for gene counts as a function of minimum TPM value we see:

|neg_min_tpm|     num_features|
|-----------:|---------------:|
|-177001| 1|
|-167039| 2|
|-163407| 3|
|-162288| 4|
|-115688| 5|
|-106147| 6|
|-94130|  7|
|-77788|  8|
|-75498|  9|
|...|....|
|-8|	35526|
|-7|	39776|
|-6|	46264|
|-5|	58324|
|-4|	84328|
|-3|	147918|
|-2|	311108|
|-1|	847297|
|0|	1388798|


The above table indicates that we have 847,297 'genes' that are expressed by at least 1 TPM in any one of the many samples in this expression matrix.  No, there are probably not so many of what we would call biologically relevant 'genes' in this data set, but instead, due to the sensitivity of RNA-Seq and our de novo transcriptome assembly, we were able to reconstruct contigs that represent that many features with evidence of being expressed at that minimum threshold.  If we increase our stringency to a minimum of 5 TPM, we report only 58,324 'genes', which many would consider a more reasonable estimate - even if still a probable exaggeration.  

Plotting the number of 'genes' (or 'transcripts') as a function of minimum TPM threshold, we can see that the vast majority of all expressed features have very little expression support.  Using R (or your own favorite data analysis package), we might extrapolate the number of expressed 'genes' based on the trend prior to the massive influx of lowly expressed transcripts:

     % R
     > data = read.table("genes_matrix.TPM.not_cross_norm.counts_by_min_TPM", header=T)
     > plot(data, xlim=c(-100,0), ylim=c(0,100000), t='b')

<img src="https://raw.githubusercontent.com/wiki/trinityrnaseq/trinityrnaseq/images/gene_count_vs_min_TPM.png" width=400 >

Another approach for exploring this is to estimate the [E90 transcript count](Transcriptome-Contig-Nx-and-ExN50-stats).

<a name='filtering-transcripts'></a>
## Filtering Transcripts Based on Expression Values

Most downstream analyses should be applied to the entire set of assembled transcripts, including functional annotation and differential expression analysis.  

If you do decide that you want to filter transcripts to exclude those that are lowly expressed, you can use the following script:

    %  $TRINITY_HOME/util/filter_low_expr_transcripts.pl 

    ##########################################################################################
    #
    #  --matrix|m <string>            expression matrix (TPM or FPKM, *not* raw counts)
    #
    #  --transcripts|t <string>       transcripts fasta file (eg. Trinity.fasta)
    #
    #
    #  # expression level filter:
    #
    #     --min_expr_any <float>      minimum expression level required across any sample (default: 0)
    #
    #  # Isoform-level filtering
    #
    #     --min_pct_dom_iso <int>         minimum percent of dominant isoform expression (default: 0)
    #          or  
    #     --highest_iso_only          only retain the most highly expressed isoform per gene (default: off)
    #                                 (mutually exclusive with --min_pct_iso param)
    #
    #     # requires gene-to-transcript mappings
    #
    #     --trinity_mode              targets are Trinity-assembled transcripts
    #         or
    #     --gene_to_trans_map <string>   file containing gene-to-transcript mappings
    #                                    (format is:   gene(tab)transcript )
    #
    #########################################################################################



The input to the script is the matrix of transcript expression values (this would ideally be your TPM matrix - or TMM-normalized TPM matrix), and your assembled transcripts fasta file.   Depending on the parameters you choose, you can filter based on minimum transcript expression levels, require some minimum percent expression of the dominant isoform, or simply retain the single most highly expressed isoform for each gene.  

>Be cautious in filtering transcripts solely based on expression values, as you can easily discard biologically relevant transcripts from your data.  Ideally, any filtering would consider a combination of expression values and functional annotation data, and filtering is currently more of an art than a science, and again, simply not needed in most circumstances unless you have a very clear objective in doing so.