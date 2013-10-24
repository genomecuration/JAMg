
all: 
	cd 3rd_party/cdbtools/cdbfasta && $(MAKE) && cp cdbfasta ../../../bin/ && cp cdbyank ../../../bin/
	cd 3rd_party/parafly && ./configure --prefix=`pwd`/../../ && $(MAKE) install
	cd 3rd_party/aatpackage && ./configure --prefix=`pwd`/../../ && $(MAKE) install

clean:
	rm -f bin/cdb* && cd 3rd_party/cdbtools/cdbfasta && $(MAKE) clean
	rm -f bin/ParaFly && cd 3rd_party/parafly && $(MAKE) clean
	cd bin && rm -f dds dps ext filter gap2 nap show && cd ../3rd_party/aatpackage && $(MAKE) clean


###################################################################


