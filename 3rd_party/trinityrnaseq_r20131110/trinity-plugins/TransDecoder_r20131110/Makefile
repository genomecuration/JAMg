SHELL := /bin/bash

parafly:
	cd 3rd_party/parafly && ./configure --prefix=`pwd`/../../util && $(MAKE) install

all: prep_pfam ffindex parafly

prep_pfam:
	if [[ -e pfam/Pfam-AB.hmm.bin && ! -e pfam/Pfam-AB.hmm.bin.h3i ]]; then echo "Preparing PFAM database" && hmmpress pfam/Pfam-AB.hmm.bin; fi

ffindex:
	cd 3rd_party/ffindex && $(MAKE) HAVE_MPI=1 && $(MAKE) install INSTALL_DIR=../../../util HAVE_MPI=1

clean:
	cd 3rd_party/parafly && $(MAKE) clean
	cd 3rd_party/ffindex && $(MAKE) clean
	rm -f pfam/Pfam-AB.hmm.bin.*
	cd sample_data && ./cleanme.pl


test:
	cd sample_data/ && ./runMe.sh

test_pfam:
	cd sample_data/ && ./__runMe_use_pfam.sh

test_ffindex:
	cd sample_data && ./__runMe_use_pfam_MPI_ffindex.sh
