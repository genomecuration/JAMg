SHELL := /bin/bash

all:
	cd 3rd_party/parafly && ./configure --prefix=`pwd`/../../util && $(MAKE) install
	cd 3rd_party/ffindex && $(MAKE) HAVE_MPI=1 && $(MAKE) install INSTALL_DIR=../../../util HAVE_MPI=1
	if [[ -e pfam/Pfam-AB.hmm.bin && ! -e pfam/Pfam-AB.hmm.bin.h3i ]]; then echo "Preparing PFAM database" && hmmpress pfam/Pfam-AB.hmm.bin; fi

clean:
	cd 3rd_party/parafly && $(MAKE) clean
	cd 3rd_party/ffindex && $(MAKE) clean
	rm -f pfam/Pfam-AB.hmm.bin.*
