## What computing resources are required?

Ideally, you will have access to a large-memory server, roughly having ~1G of RAM per 1M reads to be assembled.  The memory usage mostly depends on the complexity of the RNA-Seq data set, specifically on the number of unique k-mers.  If you do not have access to a high-memory server, [other freely available options are available](Accessing-Trinity-on-Publicly-Available-Compute-Resources).


## How long should this take?

Trinity run-time depends on a number of factors including the number of reads to be assembled and the complexity of the transcript graphs.  The assembly from start to finish can take anywhere from ~1/2 hour to 1 hour per million reads (your mileage may vary).  If you have a large data set, be sure to use the --normalize_reads parameter to greatly improve run times.
