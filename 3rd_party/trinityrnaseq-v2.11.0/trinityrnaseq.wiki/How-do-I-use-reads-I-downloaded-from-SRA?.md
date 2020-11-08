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
