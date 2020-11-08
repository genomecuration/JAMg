# Post- Transcriptome Assembly Downstream Analyses

Once your assembly is complete, there are several analyses you will likely want to pursue to explore aspects of the biology of your organism based on your assembled transcripts and the input RNA-Seq data. Such studies might include:

* [Quantifying transcript and gene abundance](Trinity-Transcript-Quantification).  This is a prerequisite to many other analyses such as examining differentially expressed transcripts among your samples.

* [Quality checking your samples and biological replicates](QC-Samples-and-Replicates). Make sure the relationships among your samples and biological replicates make sense.  If there are any confounding aspects to the data, such as outliers or batch effects, you'll want to catch them as early as possible and account for them in your further data explorations.

* [Perform differential expression analysis](Trinity-Differential-Expression). Trinity provides direct support for several DE analysis methods, including edgeR, DESeq2, Limma/Voom, and ROTS.

* Extract coding regions using [TransDecoder](Coding-Region-Identification-in-Trinity-Assemblies) and functionally annotate the transcripts using [Trinotate](Functional-Annotation-of-Transcripts).

* If your organism has an assembled genome, consider using the Trinity transcriptome assemblies for gene structure annotation using [PASA](Genome-Annotation-Using-Trinity--and-PASA).