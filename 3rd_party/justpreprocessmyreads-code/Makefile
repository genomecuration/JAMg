SHELL := /bin/bash

all:
	cd 3rd_party && unzip fastqc_v*.zip && unzip Trimmomatic-*.zip && unzip allpaths.zip && chmod +x FastQC/fastqc allpaths/* 
	cd 3rd_party && unzip blue_113.zip && cd blue_113 && find . -name "*sh" -exec /bin/bash '{}' \; && chmod +x *.exe
	ln -fs `pwd`/3rd_party/FastQC/fastqc .
clean:
	cd 3rd_party && rm -rf FastQC Trimmomatic-0.33 allpaths && cd ../test && ./cleanme.sh

prep_ecoli:
	cd dbs && if [[ ! -s ecoli_pseudomonas.fsa.masked.nr.1.bt2 ]]; then echo "Downloading contam database"; wget "http://sourceforge.net/projects/justpreprocessmyreads/files/ecoli_pseudomonas.fsa.masked.nr.bz2/download" -O ecoli_pseudomonas.fsa.masked.nr.bz2; bunzip2 ecoli_pseudomonas.fsa.masked.nr.bz2; bowtie2-build ecoli_pseudomonas.fsa.masked.nr ecoli_pseudomonas.fsa.masked.nr; fi

prep_human:
	cd dbs && if [[ ! -s human_genome.fasta.1.bt2 ]]; then echo "Downloading human genome"; wget "http://sourceforge.net/projects/justpreprocessmyreads/files/human_genome.fasta.bz2/download" -O human_genome.fasta.bz2; bunzip2 human_genome.fasta.bz2; bowtie2-build human_genome.fasta human_genome.fasta; fi

test:
	cd test && ./runme.sh
