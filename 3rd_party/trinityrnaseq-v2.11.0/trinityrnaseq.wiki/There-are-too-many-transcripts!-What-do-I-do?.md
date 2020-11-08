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
