SHELL := /bin/bash

all: 
	if [ ! -d 3rd_party/bin ]; then mkdir 3rd_party/bin; fi
	cd 3rd_party/cdbtools/cdbfasta && $(MAKE) && cp cdbfasta ../../bin/ && cp cdbyank ../../bin/ && $(MAKE) clean
	cd 3rd_party/parafly && ./configure --prefix=`pwd`/../ && if [ ! -d bin ]; then mkdir bin; fi && $(MAKE) install
	cd 3rd_party/aatpackage && ./configure --prefix=`pwd`/../ && $(MAKE) install && cp -r bin/* matrices ../bin/ && $(MAKE) clean
	cd 3rd_party/preprocess_reads && $(MAKE) 3rd_party
	cd 3rd_party/PASA && if [ ! -d bin ]; then mkdir bin; fi && $(MAKE) && cd seqclean && find . -name "*.tar.gz" -exec tar xzf '{}' \; && cd mdust && $(MAKE) && cd ../psx && $(MAKE) && cd ../trimpoly && $(MAKE) && cd .. && find . -type f -executable -exec cp -u '{}' ../../bin/ \;
	cd 3rd_party/transdecoder && $(MAKE) nopfam
	cd 3rd_party/trinityrnaseq && $(MAKE)
	cd 3rd_party/RepeatMasker && if [ -e Libraries/Dfam.hmm.bz2 ]; then bunzip2 -vd Libraries/Dfam.hmm.bz2; fi && if [ -e taxonomy.dat.bz2 ]; then bunzip2 -vd taxonomy.dat.bz2 ; fi && echo "I will download TRF from http://tandem.bu.edu/trf. You must read and accept their (v. short) license (see http://tandem.bu.edu/trf/trf.license.html) otherwise exit now:" && sleep 5 && echo "The author of this software grants to any individual or organization the right to use and to make an unlimited number of copies of this software. You may not de-compile, disassemble, reverse engineer, or modify the software. This software cannot be sold, incorporated into commercial software or redistributed. The author of this software accepts no responsibility for damages resulting from the use of this software and makes no warranty or representation, either express or implied, including but not limited to, any implied warranty of merchantability or fitness for a particular purpose. This software is provided as is, and the user assumes all risks when using it." && sleep 10 && wget -r http://tandem.bu.edu/trf/downloads/trf407b.linux64 -O trf && chmod +x trf && echo "I will configure RepeatMasker. Use the RMBlast NCBI option when asked 'Add a Search Engine', the files are in `pwd`/ncbi-blast"
	cd 3rd_party/augustus && $(MAKE) && cp bin/* ../bin/
	cd 3rd_party/geneid && $(MAKE)
	cd 3rd_party/cegma && $(MAKE)
	cd 3rd_party/GlimmerHMM/sources && $(MAKE) && cd ../train && $(MAKE)
	cd 3rd_party/snap && $(MAKE) && cp fathom ../bin/
	cd 3rd_party/gmap && ./configure --prefix=`pwd`/../ --with-gmapdb=`pwd`/../../databases/gmap/ && $(MAKE) && $(MAKE) check && $(MAKE) install
	cd 3rd_party/samtools && $(MAKE) && cp samtools ../bin/ && $(MAKE) clean
	cd 3rd_party/bedtools && if [ ! -d bin ]; then mkdir bin; fi && $(MAKE) && cp bin/* ../bin/ && $(MAKE) clean
	cd databases/hhblits && echo "Uncompressing databases, this may take a while..." find . -name "*tar.bz2" -exec tar -xjf '{}' \;
	chmod a+rx bin/* 3rd_party/bin/*
clean:
	rm -f bin/cdb* && cd 3rd_party/cdbtools/cdbfasta && $(MAKE) clean
	rm -f bin/ParaFly && cd 3rd_party/parafly && $(MAKE) clean
	cd 3rd_party/snap && $(MAKE) clean
	cd 3rd_party/transdecoder && $(MAKE) clean
	cd 3rd_party/RepeatMasker && if [ -e Libraries/Dfam.hmm ] ; then bzip2 -v Libraries/Dfam.hmm; fi && if [ -e taxonomy.dat ] ; then bzip2 -v taxonomy.dat; fi && rm -f trf
	cd 3rd_party/augustus && $(MAKE) clean
	cd 3rd_party/geneid && $(MAKE) clean
	cd 3rd_party/cegma && if [ -e bin ] ; then $(MAKE) clean; fi
	cd 3rd_party/snap && $(MAKE) clean
	cd 3rd_party/gmap && $(MAKE) clean
	cd 3rd_party/samtools && $(MAKE) clean
	cd 3rd_party/bedtools && $(MAKE) clean
	cd 3rd_party/GlimmerHMM/sources && rm -f *.o
	cd 3rd_party/PASA && $(MAKE) clean && cd seqclean && rm -rf mdust psx seqclean trimpoly
	cd 3rd_party/aatpackage && $(MAKE) clean
	cd 3rd_party/bin && rm -fr *
	cd test_suite && bash cleanme.sh

test:
	cd test_suite && bash runme.sh && bash cleanme.sh

###################################################################


