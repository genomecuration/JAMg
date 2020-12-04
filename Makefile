SHELL := /bin/bash

all: 
	if [ ! -d 3rd_party/bin ]; then mkdir 3rd_party/bin; fi
	if [ ! -d 3rd_party/lib ]; then mkdir 3rd_party/lib; fi
	if [ ! -d 3rd_party/share ]; then mkdir 3rd_party/share; fi
	if [ ! -d 3rd_party/include ]; then mkdir 3rd_party/include; fi
	if [ ! -d 3rd_party/lib64 ]; then mkdir 3rd_party/lib64; fi
	export JAMG_PATH=$(PWD)
	export LD_LIBRARY_PATH="$(LD_LIBRARY_PATH):$(PWD)/3rd_party/lib:$(PWD)/3rd_party/lib64:$(PWD)/3rd_party/mysql/lib"
	export LDFLAGS="$(LDFLAGS) -L$(PWD)/3rd_party/lib -L$(PWD)/3rd_party/lib64 -L$(PWD)/3rd_party/mysql/lib"
	export CPPFLAGS="$(CPPFLAGS) -I$(PWD)/3rd_party/include -I$(PWD)/3rd_party/mysql/include -fPIC"
	export CFLAGS="$(CFLAGS) -I$(PWD)/3rd_party/include -I$(PWD)/3rd_party/mysql/include -fPIC"
	cd 3rd_party && tar xzf zlib-1.2.8.tar.gz && cd zlib-1.2.8 && ./configure --prefix=`pwd`/../ --64 && $(MAKE) && $(MAKE) install
	cd 3rd_party && tar xzf libpng-1.6.19.tar.gz && cd libpng-1.6.19 && ./configure --prefix=`pwd`/../ --with-zlib-prefix=`pwd`/../ && $(MAKE) && $(MAKE) install
	cd 3rd_party/mysql && find . -name "*bz2" -exec bunzip2 '{}' \; 
	cd 3rd_party/EMBOSS && if [[ ! -e emboss/data/TAXONOMY/division.dmp ]]; then bunzip2 -k emboss/data/TAXONOMY/*bz2; fi && ./configure --prefix=`pwd`/../ --without-mysql --without-x --without-java --without-hpdf --without-auth --without-axis2c --without-postgresql --without-pngdriver && $(MAKE) -j3 && $(MAKE) install
	if [ ! -d 3rd_party/hhsuite ]; then cd 3rd_party && git clone --depth 1 https://github.com/soedinglab/hh-suite.git hhsuite && cd hhsuite && git submodule init; else cd 3rd_party/hhsuite && git pull; fi 
	cd 3rd_party/hhsuite && git submodule update && rm -rf build && mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -G "Unix Makefiles" -DCMAKE_INSTALL_PREFIX=../ .. && $(MAKE) && $(MAKE) install
	cd 3rd_party/hhsuite && cp -fp lib64/* ../lib64/ && cp -fp include/* ../include/ && cp -fp bin/* ../bin/
	cd 3rd_party/ncbi-blast/bin && ln -fs `pwd`/* ../../bin/
	cd 3rd_party/cdbtools/cdbfasta && $(MAKE) && cp -pf cdbfasta ../../bin/ && cp -pf cdbyank ../../bin/ && $(MAKE) clean
	cd 3rd_party && tar -xjf samtools-1.3.tar.bz2 && ln -sf samtools-1.3 samtools && cd samtools && ./configure --without-curses --prefix=`pwd`/../ && $(MAKE) && $(MAKE) install
	cd 3rd_party/bedtools && if [ ! -d bin ]; then mkdir bin; fi && $(MAKE) -j 3 && cp -pf bin/* ../bin/ && $(MAKE) clean
	cd 3rd_party/blat && ln -fs `pwd`/* ../bin/
	cd 3rd_party/parafly && ./configure --prefix=`pwd`/../ && if [ ! -d bin ]; then mkdir bin; fi && $(MAKE) install
	cd 3rd_party/aatpackage && ./configure --prefix=`pwd`/../ && $(MAKE) install && cp -fpr bin/* matrices ../bin/ && $(MAKE) clean
	cd 3rd_party/justpreprocessmyreads && $(MAKE)
	cd 3rd_party/PASA && if [ ! -d bin ]; then mkdir bin; fi && $(MAKE) && cd seqclean && wget -qc ftp://ftp.ncbi.nih.gov/pub/UniVec/UniVec && find . -name "*.tar.gz" -exec tar xzf '{}' \; && rm -rf seqclean/bin && cd mdust && $(MAKE) && cd ../psx && $(MAKE) && cd ../trimpoly && $(MAKE) && cd .. && find . -type f -executable -exec cp -pf '{}' ../../bin/ \;
	cd 3rd_party/transdecoder && $(MAKE)
	cd 3rd_party/trinityrnaseq && $(MAKE)
	cd 3rd_party && tar -xjf mira_4.9.5_2_linux-gnu_x86_64_static.tar.bz2 && ln -sf mira_4.9.5_2_linux-gnu_x86_64_static mira
	cd 3rd_party/augustus && $(MAKE) && cp -pf bin/* ../bin/ && cp -fp scripts/gff2gbSmallDNA.pl scripts/filterGenes.pl ../bin/ 
	cd 3rd_party/geneid && $(MAKE)
	cd 3rd_party/GlimmerHMM/sources && $(MAKE) && cd ../train && $(MAKE)
	cd 3rd_party/snap && $(MAKE) && cp -fp fathom ../bin/
	cd 3rd_party/gmap && ./configure --prefix=`pwd`/../ --with-gmapdb=`pwd`/../../databases/gmap/ && $(MAKE) -j3 && $(MAKE) check && $(MAKE) install
	cd databases/hhblits && printf "\n\nUncompressing databases, this may take a while...\n\n" find . -name "*tar.bz2" -exec tar -xjf '{}' \;
	cd 3rd_party/fasta && cp -fp bin/* ../bin/ && cd ../bin/ && ln -sf fasta36 fasta
	cd 3rd_party/exonerate/bin && ln -fs `pwd`/* ../../bin/
	cd 3rd_party/tRNAscan-SE && $(MAKE) && $(MAKE) install
	cd 3rd_party/aragorn && gcc -O3 -o aragorn aragorn1.2.36.c && cp -fp aragorn ../bin/
	cd 3rd_party/minimap2/ && make && cp -f minimap2 ../bin/
	cd 3rd_party/RECON/src && $(MAKE) && $(MAKE) install && cp -fp ../bin/* ../../bin/
	cd 3rd_party/RepeatScout && $(MAKE)
	cd 3rd_party/RepeatMasker && \
  if [[ ! -e Libraries/Dfam.hmm ]]; \
  then wget -c -O Libraries/Dfam.hmm.gz http://dfam.org/web_download/Current_Release/Dfam.hmm.gz && gunzip Libraries/Dfam.hmm.gz; fi && \
  if [[ -e Libraries/taxonomy.dat.bz2 && ! -e Libraries/taxonomy.dat ]]; then bunzip2 -kd Libraries/taxonomy.dat.bz2 ; fi && \
  if [[ ! -e trf ]]; then printf "\n\nI will download TRF from http://tandem.bu.edu/trf. You must read and accept their (v. short) license (see http://tandem.bu.edu/trf/trf.license.html) otherwise exit now:\n" && sleep 5 && printf "\n\nThe author of this software grants to any individual or organization the right to use and to make an unlimited number of copies of this software. You may not de-compile, disassemble, reverse engineer, or modify the software. This software cannot be sold, incorporated into commercial software or redistributed. The author of this software accepts no responsibility for damages resulting from the use of this software and makes no warranty or representation, either express or implied, including but not limited to, any implied warranty of merchantability or fitness for a particular purpose. This software is provided as is, and the user assumes all risks when using it.\n" && sleep 10 && wget -c -r http://tandem.bu.edu/trf/downloads/trf407b.linux64 -O trf; fi && chmod +x trf && \
  printf "\n\nYou have to now configure RepeatMasker. Answer the questions. DEFAULTS SHOULD WORK. Can configure RMBLAST AND HMMER. Use the RMBlast NCBI option as DEFAULT.\n" && sleep 10 && ./configure
	cd 3rd_party/RepeatModeler && printf "\n\nYou have to now configure RepeatModeler.  Answer the questions. DEFAULTS SHOULD WORK. Use RMBLAST NCBI option as DEFAULT\n" && sleep 5 && ./configure ;
	chmod -R a+rx bin 3rd_party/bin
	printf "\n\nInstallation complete.\n\n"
clean:
	cd 3rd_party/cdbtools/cdbfasta && $(MAKE) clean
	cd 3rd_party/parafly && $(MAKE) clean
	cd 3rd_party/snap && $(MAKE) clean
	cd 3rd_party/RepeatMasker && rm -f trf Libraries/taxonomy.dat Libraries/Dfam.hmm 
	cd 3rd_party/augustus && $(MAKE) clean
	cd 3rd_party/geneid && $(MAKE) clean
	cd 3rd_party/snap && $(MAKE) clean
	cd 3rd_party/gmap && $(MAKE) clean
	cd 3rd_party/samtools && $(MAKE) clean
	cd 3rd_party/bedtools && $(MAKE) clean
	cd 3rd_party/GlimmerHMM/sources && rm -f *.o
	cd 3rd_party/PASA && $(MAKE) clean && cd seqclean && rm -rf mdust psx seqclean trimpoly
	cd 3rd_party/aatpackage && $(MAKE) clean
	rm -fr 3rd_party/bin 
	cd test_suite && bash cleanme.sh

test:
	cd test_suite && bash runme.sh && bash cleanme.sh

###################################################################


