SHELL := /bin/bash

all: 
	cd 3rd_party/cdbtools/cdbfasta && $(MAKE) && cp cdbfasta ../../../bin/ && cp cdbyank ../../../bin/
	cd 3rd_party/parafly && ./configure --prefix=`pwd`/../../ && $(MAKE) install
	cd 3rd_party/aatpackage && ./configure --prefix=`pwd`/../../ && $(MAKE) install && cp bin/* ../../bin/
	chmod a+rx bin/*
	cd 3rd_party/transdecoder && $(MAKE)
	cd 3rd_party/augustus && $(MAKE)
	cd 3rd_party/geneid && $(MAKE)
	cd 3rd_party/snap && $(MAKE) && cp fathom ../../bin/
	cd 3rd_party/gmap && ./configure --prefix=`pwd`/../../ --with-gmapdb=`pwd`/../../databases/ && $(MAKE) && $(MAKE) check && $(MAKE) install
clean:
	rm -f bin/cdb* && cd 3rd_party/cdbtools/cdbfasta && $(MAKE) clean
	rm -f bin/ParaFly && cd 3rd_party/parafly && $(MAKE) clean
	cd bin && rm -f fathom dds dps ext filter gap2 nap show AAT.pl AAT_btab_to_alignment_gff3.pl alignment_gff3_to_bed_format.pl extCollapse.pl atoiindex cmetindex dbsnp_iit fa_coords get-genome gff3_* gmap* gsnap* gtf_* iit_* gvf_iit md_coords psl_* snpindex uniqscan uniqscanl vcf_iit && cd ../3rd_party/aatpackage && $(MAKE) clean
	cd 3rd_party/snap && $(MAKE) clean
	cd 3rd_party/transdecoder && $(MAKE) clean
	cd 3rd_party/augustus && $(MAKE) clean
	cd 3rd_party/geneid && $(MAKE) clean
	cd 3rd_party/snap && $(MAKE) clean
	cd 3rd_party/gmap && $(MAKE) clean && rm -fr autom4te.cache config.log
	cd test_suite && bash cleanme.sh

test:
	cd test_suite && bash runme.sh && bash cleanme.sh

###################################################################


