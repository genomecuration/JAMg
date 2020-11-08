#Trinity Hardware and Configuration Requirements

The Inchworm and Chrysalis steps can be memory intensive.  A basic recommendation is to have ~1G of RAM per ~1M pairs of Illumina reads. Simpler transcriptomes (lower eukaryotes) require less memory than more complex transcriptomes such as from vertebrates.  

If you are able to run the entire Trinity process on a single high-memory multi-core server, indicate the number of butterfly processes to run in parallel by the --CPU parameter. 

Our experience is that the entire process can require ~1/2 hour to one hour per million pairs of reads in the current implementation (see link:[FAQ](Trinity-FAQ)).  We're striving to improve upon both memory and time requirements.


If you do not have direct access to a high memory machine (typically having 256G or 512G of RAM), consider [running Trinity on one of the externally available resources](Accessing-Trinity-on-Publicly-Available-Compute-Resources).
