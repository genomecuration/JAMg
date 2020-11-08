# Trinity: Frequently Asked Questions

- [There are too many transcripts! What do I do?](#ques_why_so_many_transcripts)
- [What computing resources are required?](#ques_comp_resources_required)
- [How long should this take?](#ques_how_long)
- [How do I use reads I downloaded from SRA?](#ques_sra_fq_conversion)
- [How can I run this in parallel on a computing grid?](#ques_computing_grid)
- [How do I identify the specific reads that were incorporated into the transcript assemblies?](#ques_reads_in_assembly)
- [How how do I combine paired and single reads?](#ques_mult_seq_library_types)
- [How do I combine strand-specific and non-strand-specific reads?](#ques_ss_and_not_ss)
- [Trinity process died due to 'std::bad_alloc'](#ques_bad_alloc)
- [Butterfly fails with java Error: Cannot create GC thread. Out of system resources.](#ques_butterfly_GC_thread_fail)
- [Last resort - sharing your data privately with Trinity developers for debugging](#share_data_trin_devel)

<a name="ques_why_so_many_transcripts"></a>
## There are too many transcripts!  What do I do?

Welcome to RNA-Seq de novo assembly!  :)  There are several aspects of de novo assembled transcripts to be aware of:

-  *Lots* of transcripts is the rule rather than the exception.  

-  Most of the transcripts are very lowly expressed, and the deeper you sequence and the more complex your genome, the larger the number of lowly expressed transcripts you will be able to assemble.  Biological relevance of the lowly expressed transcripts could be questionable - some are bound to be very relevant.

-  There's really no good reason to immediately filter them out.  They can be 'passengers' throughout all of your data analyses, and if any of them are important, they'll ideally surface in the relevant study.   You can put them all through Trinotate.github.io for annotation/analysis, and you can put them through DE studies just fine (those with insufficent reads will get directly filtered out during the DE analysis protocols to avoid problems associated with multiple hypothesis testing);  If the read counts are few or lacking, they simply won't surface as significant DE entries, but if there's protein homology or other interesting features, you'll want to continue to capture this info - hence don't feel the need to immediately filter!

-  If you want to guestimate 'how many expressed genes/transcripts do I really have?', then examine [counts of transcripts or genes vs. min expression thresholds](Trinity-Transcript-Quantification). This will enable you to count the number that reflect the majority of the reads, excluding counting those that have little read support.  The entries with >= ~1 fpkm or tpm tend to be heavily enriched for transcripts that correspond to what we typically think of as 'genes' in the pre-RNA-Seq era, and what typically are awarded annotation status in genome annotation. 

- Related to the above point, consider computing the [contig Ex90N50 value](Transcriptome-Contig-Nx-and-ExN50-stats) (contig N50 computed based on transcript contigs representing the top 90% of expressed transcripts), in addition to the Ex90 count of transcripts (the number of transcripts that represent this top 90% of expression data).  You'll find these numbers to be more palatable as compared to the hundreds of thousands (or millions) of transcript contigs being reconstructed - which we well know are biased towards any lowly expressed transcripts that can be reconstructed.

-  If for some reason you do want to filter out those transcripts that have little read support, you can [filter them based on the abundance estimates](Trinity-Transcript-Quantification#filtering-transcripts).

-  [CD-HIT](http://weizhongli-lab.org/cd-hit/) can be used to cluster highly similar sequences and retain a single representative sequence per group. Stringent clustering can be done like so:   'cd-hit-est -o cdhit -c 0.98 -i Trinity.fasta -p 1 -d 0 -b 3 -T 10'

-  If you're concerned about having too many isoforms for a given gene, you can use the filtering method above to remove isoforms that are weakly expressed compared to more dominant isoforms.  Also, you can simply perform DE analysis at the 'gene' level in addition to the 'isoform' level, so as to increase your power for DE analysis.  Other methods such as [Corset](http://genomebiology.com/2014/15/7/410) can also be used to regroup relevant transcripts into 'gene' clusters. Since the 'gene' aggregates the expression data from all of the isoforms, problems with having too many isoforms quickly dissipate.


-  If you have compelling reasons for filtering, please share them with us on our [Trinity google forum](https://groups.google.com/forum/#!forum/trinityrnaseq-users), and we'll update our documentation accordingly!

<a name="ques_comp_resources_required"></a>
## What computing resources are required?

Ideally, you will have access to a large-memory server, roughly having ~1G of RAM per 1M reads to be assembled.  The memory usage mostly depends on the complexity of the RNA-Seq data set, specifically on the number of unique k-mers.  If you do not have access to a high-memory server, [other freely available options are available](Accessing-Trinity-on-Publicly-Available-Compute-Resources).

<a name="ques_how_long"></a>
## How long should this take?

Trinity run-time depends on a number of factors including the number of reads to be assembled and the complexity of the transcript graphs.  The assembly from start to finish can take anywhere from ~1/2 hour to 1 hour per million reads (your mileage may vary).  If you have a large data set, be sure to use the --normalize_reads parameter to greatly improve run times.



<a name="ques_sra_fq_conversion"></a>
## How do I use reads I downloaded from SRA?

RNA-Seq data downloaded from SRA tends to exist in a .sra file that needs to be converted to fastq file format.  This can be done using the [SRA Toolkit](http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software) like so:

      SRA_TOOLKIT/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files file.sra

>The '@' and '$' in the above should be provided exactly as those characters.  They're not symbols for you to personally interpolate here.

Instead of file.sra, just give it your SRR accession, and it'll retrieve the data for you.  This requires a modern version of the fastq-dump utility, which you can get from here <https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software>.

For example:
```
  SRA_TOOLKIT/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files SRR390728
```


### Other example strategies for downloading directly from SRA via fastq-dump:

The command below will generate an 'interleaved' fastq file, where record 1 is followed immediately by record 2, and we'll extract the top 1M read pairs (which = 8M top lines due to the interleaving).

>below, we retrieve SRA accession: SRR390728 and use the -X 5 parameter to simply limit the number of reads in this small example.  When you run it with your sample of interest, use your SRR-value and do not use the -X 5 parameter.

    % fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files -X 5 -Z  SRR390728 | \
           head -n8000000 | gzip > SRR390728.interleaved.fastq.gz

Now, to de-interleave and generate the two separate fastq files for the 'left' and 'right' read mates, we can do the following:

    % gunzip -c SRR390728.interleaved.fastq.gz | \
      paste - - - - - - - - | \
      tee >(cut -f 1-4 | tr '\t' '\n' | gzip > SRR390728_1.fastq.gz) | \
      cut -f 5-8 | tr '\t' '\n' | gzip -c > SRR390728_2.fastq.gz

>above is adapted from: https://biowize.wordpress.com/2015/03/26/the-fastest-darn-fastq-decoupling-procedure-i-ever-done-seen/

>If the read names contain _forward or _reverse in the names, this must be removed somehow (such as via a sed command).

<a name="ques_computing_grid"></a>
## How can I run this in parallel on a computing grid?

The initial Inchworm and Chrysalis steps currently need to be run on a single server, however parts of Chrysalis and all of Butterfly can be run as naive parallel processes on a computing grid. Both SGE and LSF are currently supported. Configuring Trinity for leveraging these grid-computing platforms is described [here](Running-Trinity#grid_conf).


<a name="ques_reads_in_assembly"></a>
## How do I identify the specific reads that were incorporated into the transcript assemblies?

Currently, the mappings of reads to transcripts are not reported.  To obtain this information, we recommend realigning the reads to the assembled transcripts [Bowtie](using http://bowtie-bio.sourceforge.net/index.shtml). Capturing both properly paired read alignments and single unpaired alignments can be captured and quantified as described [here](RNA-Seq-Read-Representation-by-Trinity-Assembly).


<a name="ques_mult_seq_library_types"></a>
## How do I combine paired-end and single-end reads?

Our recommendation for doing this has changed over time and across different Trinity software versions.  Trinity doesn't have great support for using mixed single-end and paired-end data, and the setup process for shoe-horning the data into Trinity involves a few manual steps.

1. If you want to quality trim your data, run [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (or other) directly.

2. Run [in silico normalization](Trinity-Insilico-Normalization) separately for each of your libraries.

3. Combine all your reads into one fastq file.  If your data are strand-specific, orient your reads to all the sense transcript orientation.

4. Run Trinity like so:

    %  Trinity --single combined_reads.fastq --no_normalize_reads  --run_as_paired  ...

>if your data are strand-specific and were oriented to the forward strand orientation, be sure to include '--SS_lib_type F'

<a name="ques_ss_and_not_ss"></a>
## How do I combine strand-specific and non-strand-specific reads?

There is no good way to combine strand-specific data with non-strand-specific data, unless you decide to treat the entire data set as non-strand-specific.


<a name="ques_bad_alloc"></a>
## Trinity process died due to 'std::bad_alloc'

This is an indicator that the process ran out of available RAM. If you have more RAM resources to make available to Trinity, then simply rerun your original Trinity command with the altered resources allocated and it should resume approximately where it left off.  

If you are resource limited, please consider running Trinity [here](Accessing-Trinity-on-Publicly-Available-Compute-Resources).  If you want to continue to try to run Trinity given your available resources, you can reduce the total RAM requirements by running Trinity with parameter '--min_kmer_cov 2'. Although the assembly should still be of high quality and require less RAM, lowly expressed transcripts may be more highly fragmented in the assembly.


<a name="ques_butterfly_GC_thread_fail"></a>
## Butterfly fails with java Error: Cannot create GC thread. Out of system resources.

There are a couple reasons why this error message might creep up.

1.  *all memory has been consumed on the machine*.  Each butterfly process wants to reserve 10G of maximum heap space.  If there's less than 10G of free memory on the machine per butterfly (--CPU setting), then java may not be able to initialize (depends on your OS configuration).  Try reducing the --CPU setting and rerunning your Trinity command. It should resume where it left off.

2.  Your server might be configured to allow only limited numbers of processes per user.  Check your hardware configuration like so:

```
    %  cat /etc/security/limits.conf 

     soft    memlock         unlimited
     hard    memlock         unlimited
     soft    stack           unlimited
     hard    stack           unlimited
        -    nofile          131072
     soft    nproc           unlimited
     hard    nproc           unlimited
```

There are various ways to deal with restricted settings (Google can help).


3.  *NUMA architecture*:  one of our users found that the java invocation required: -XX:ParallelGCThreads=<Numerical Thread Count>, otherwise it would try to use too many threads.

<a name="share_data_trin_devel"></a>
## Last resort - sharing your data privately with Trinity developers for debugging

A great way to share your data with the Trinity developers is through:
<https://mega.nz/>
as, a command-line interface can be used to easily download data from there.

Alternatively, [google drive](https://www.google.com/drive/) or [dropbox](http://www.dropbox.com) can be used.

Get the Trinity developers direct email address, rather than posting links to your data on the google forum for privacy concerns.  Trinity developers will not share your data, so no worries!

