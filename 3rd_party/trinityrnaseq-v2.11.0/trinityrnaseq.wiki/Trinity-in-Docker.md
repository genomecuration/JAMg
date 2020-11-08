# Run Trinity Using Docker

<img src="https://s3.amazonaws.com/media-p.slid.es/uploads/602799/images/3238125/moby.svg" width=500 />

If you have [Docker](https://www.docker.com/) installed, you can pull [our image from DockerHub](https://hub.docker.com/r/trinityrnaseq/trinityrnaseq/), which contains Trinity and all software used for downstream analyses supported within the larger Trinity framework.

Pull the latest Docker image for Trinity like so:

    % docker pull trinityrnaseq/trinityrnaseq


Run Trinity like so (eg. as shown where with a very small test data set):

    % docker run --rm -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq Trinity \
          --seqType fq \
          --left `pwd`/reads_1.fq.gz \
          --right `pwd`/reads_2.fq.gz
          --max_memory 1G --CPU 4 --output `pwd`/trinity_out_dir


## Downstream analyses using Dockerized Trinity:

Trinity is installed in the Docker container at '/usr/local/bin/trinityrnaseq'.

Just use that path to access all tools installed within Trinity.

eg.

     %  docker run -v`pwd`:`pwd` trinityrnaseq/trinityrnaseq \
          /usr/local/bin/trinityrnaseq/util/align_and_estimate_abundance.pl


>With the above, just be sure to specify full paths to inputs and outputs.

<a name='trinity_singularity'></a>
# Running Trinity Using Singularity

<img src="https://sylabs.io/assets/images/logos/singularity.png" width=200 >

Singularity is easier and safer to use than Docker, and is our preferred method for running Trinity. All modern releases of Trinity have a Singularity image (.simg) offered for download from our [Trinity Singularity Image Archive](https://data.broadinstitute.org/Trinity/TRINITY_SINGULARITY/).  If you have [Singularity](https://sylabs.io/docs/) installed and the .simg file downloaded, you can run Trinity like so:

```
    %  singularity exec -e Trinity.simg  Trinity \
          --seqType fq \
          --left `pwd`/reads_1.fq.gz  \
          --right `pwd`/reads_2.fq.gz \
          --max_memory 1G --CPU 4 \
          --output `pwd`/trinity_out_dir
```

All downstream analyses can be accessed similarly to the Docker instructions above, but using the Singularity execution syntax.

