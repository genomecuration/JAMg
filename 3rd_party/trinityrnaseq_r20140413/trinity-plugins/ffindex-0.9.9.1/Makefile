OS:= $(shell uname)

ifeq ($(OS), Darwin)
MFILE=Makefile.osx
else
MFILE=Makefile
endif

all:
	$(MAKE) -C src -f $(MFILE) $@

%:
	$(MAKE) -C src -f $(MFILE) $@

clean:
	$(MAKE) -C src -f $(MFILE) clean && rm -f bin/*

release:
	hg archive -t tgz ffindex-`cat VERSION`.tar.gz
