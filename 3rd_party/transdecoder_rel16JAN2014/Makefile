SHELL := /bin/bash

all: prep_pfam ffindex parafly cdhit

nopfam: ffindex parafly cdhit

parafly:
	cd 3rd_party/parafly && ./configure --prefix=`pwd`/../../util && $(MAKE) install

prep_pfam:
	if [[ ! -s pfam/Pfam-AB.hmm.bin ]]; then echo "Downloading PFAM database" && wget 'http://sourceforge.net/projects/transdecoder/files/Pfam-AB.hmm.bin/download' -O pfam/Pfam-AB.hmm.bin; fi
	if [[ ! -e pfam/Pfam-AB.hmm.bin.h3i ]]; then echo "Preparing PFAM database" && hmmpress pfam/Pfam-AB.hmm.bin; fi

ffindex:
	cd 3rd_party/ffindex && $(MAKE) HAVE_MPI=1 && $(MAKE) install INSTALL_DIR=../../../util HAVE_MPI=1
	cd 3rd_party/ffindex && $(MAKE) install INSTALL_DIR=../../../util HAVE_MPI=1

cdhit:
	cd 3rd_party/cd-hit && $(MAKE) openmp=yes && $(MAKE) install PREFIX=../../util/bin

clean:
	cd 3rd_party/parafly && $(MAKE) clean
	cd 3rd_party/ffindex && $(MAKE) clean
	cd 3rd_party/cd-hit && $(MAKE) clean
	rm -f util/bin/* util/lib64/* util/lib/* util/include/*
#	rm -f pfam/Pfam-AB.hmm.bin.*
	cd sample_data && ./cleanme.pl

test:
	cd sample_data/ && ./runMe.sh

test_pfam:
	cd sample_data/ && ./__runMe_use_pfam.sh

test_ffindex:
	cd sample_data && ./__runMe_use_pfam_MPI_ffindex.sh
