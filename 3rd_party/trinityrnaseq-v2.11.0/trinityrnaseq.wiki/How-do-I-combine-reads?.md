## How do I combine paired-end and single-end reads?

Our recommendation for doing this has changed over time and across different Trinity software versions.  Trinity doesn't have great support for using mixed single-end and paired-end data, and the setup process for shoe-horning the data into Trinity involves a few manual steps.

1. If you want to quality trim your data, run [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic) (or other) directly.

2. Run [in silico normalization](Trinity-Insilico-Normalization) separately for each of your libraries.

3. Combine all your reads into one fastq file.  If your data are strand-specific, orient your reads to all the sense transcript orientation.

4. Run Trinity like so:

    %  Trinity --single combined_reads.fastq --no_normalize_reads  --run_as_paired  ...

>if your data are strand-specific and were oriented to the forward strand orientation, be sure to include '--SS_lib_type F'


## How do I combine strand-specific and non-strand-specific reads?

There is no good way to combine strand-specific data with non-strand-specific data, unless you decide to treat the entire data set as non-strand-specific.
