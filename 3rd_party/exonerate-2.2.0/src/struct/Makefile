# Makefile.in generated automatically by automake 1.4-p6 from Makefile.am

# Copyright (C) 1994, 1995-8, 1999, 2001 Free Software Foundation, Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.


SHELL = /bin/bash

srcdir = .
top_srcdir = ../..

prefix = /loginhome/jzucker/src/jamg-code/3rd_party/exonerate/../
exec_prefix = ${prefix}

bindir = ${exec_prefix}/bin
sbindir = ${exec_prefix}/sbin
libexecdir = ${exec_prefix}/libexec
datadir = ${datarootdir}
sysconfdir = ${prefix}/etc
sharedstatedir = ${prefix}/com
localstatedir = ${prefix}/var
libdir = ${exec_prefix}/lib
infodir = ${datarootdir}/info
mandir = ${datarootdir}/man
includedir = ${prefix}/include
oldincludedir = /usr/include

DESTDIR =

pkgdatadir = $(datadir)/exonerate
pkglibdir = $(libdir)/exonerate
pkgincludedir = $(includedir)/exonerate

top_builddir = ../..

ACLOCAL = /loginhome/jzucker/src/jamg-code/3rd_party/exonerate/missing aclocal-1.4
AUTOCONF = autoconf
AUTOMAKE = /loginhome/jzucker/src/jamg-code/3rd_party/exonerate/missing automake-1.4
AUTOHEADER = autoheader

INSTALL = /usr/bin/install -c
INSTALL_PROGRAM = ${INSTALL} $(AM_INSTALL_PROGRAM_FLAGS)
INSTALL_DATA = ${INSTALL} -m 644
INSTALL_SCRIPT = ${INSTALL}
transform = s,x,x,

NORMAL_INSTALL = :
PRE_INSTALL = :
POST_INSTALL = :
NORMAL_UNINSTALL = :
PRE_UNINSTALL = :
POST_UNINSTALL = :
host_alias = 
host_triplet = x86_64-unknown-linux-gnu
CC = gcc
GLIB_CONFIG = 
HAVE_LIB = @HAVE_LIB@
LIB = @LIB@
LTLIB = @LTLIB@
MAKEINFO = /loginhome/jzucker/src/jamg-code/3rd_party/exonerate/missing makeinfo
PACKAGE = exonerate
PKG_CONFIG = /usr/bin/pkg-config
VERSION = 2.2.0
codegen_extra_ldadd = viterbi.o scheduler.o c4_model_archive.a
codegen_extra_sources = c4_model_archive.h
datarootdir = ${prefix}/share
glib_cflags = -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include  
glib_libs = -lglib-2.0  
host = x86_64-unknown-linux-gnu
installed_util_list =                                                fastaclean fastaclip fastachecksum fastacomposition            fastadiff fastaexplode fastafetch fastahardmask                fastaindex fastalength fastanrdb fastaoverlap                  fastareformat fastaremove fastarevcomp fastasoftmask           fastasort fastasplit fastasubseq fastatranslate                fastavalidcds fasta2esd fastaannotatecdna esd2esi
source_root_dir = /loginhome/jzucker/src/jamg-code/3rd_party/exonerate
use_pthreads = 

TESTS = fsm.test matrix.test pqueue.test slist.test rangetree.test         vfsm.test recyclebin.test sparsecache.test dejavu.test             bitarray.test splaytree.test noitree.test

noinst_PROGRAMS = $(TESTS)

noinst_HEADERS = fsm.h matrix.h pqueue.h slist.h rangetree.h vfsm.h                  recyclebin.h sparsecache.h dejavu.h bitarray.h                      splaytree.h noitree.h


fsm_test_SOURCES = fsm.test.c fsm.c
fsm_test_LDADD = recyclebin.o

matrix_test_SOURCES = matrix.test.c matrix.c
pqueue_test_SOURCES = pqueue.test.c pqueue.c recyclebin.c
slist_test_SOURCES = slist.test.c slist.c recyclebin.c
rangetree_test_SOURCES = rangetree.test.c rangetree.c
rangetree_test_LDADD = recyclebin.o

vfsm_test_SOURCES = vfsm.test.c vfsm.c
vfsm_test_LDADD = -lm

recyclebin_test_SOURCES = recyclebin.test.c recyclebin.c
sparsecache_test_SOURCES = sparsecache.test.c sparsecache.c
dejavu_test_SOURCES = dejavu.test.c dejavu.c
bitarray_test_SOURCES = bitarray.test.c bitarray.c

splaytree_test_SOURCES = splaytree.test.c splaytree.c
splaytree_test_LDADD = recyclebin.o

noitree_test_SOURCES = noitree.test.c noitree.c
noitree_test_LDADD = recyclebin.o splaytree.o

# Files to clear away

MAINTAINERCLEANFILES = Makefile.in
mkinstalldirs = $(SHELL) $(top_srcdir)/mkinstalldirs
CONFIG_CLEAN_FILES = 
PROGRAMS =  $(noinst_PROGRAMS)


DEFS = -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE=\"exonerate\" -DVERSION=\"2.2.0\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIMITS_H=1 -DHAVE_ERRNO_H=1 -DHAVE_MATH_H=1 -DHAVE_STRINGS_H=1 -DHAVE_CTYPE_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_TIME_H=1 -DHAVE_UNISTD_H=1 -I. -I$(srcdir) 
CPPFLAGS = 
LDFLAGS =  -lpthread
LIBS =  -lglib-2.0  
fsm_test_OBJECTS =  fsm.test.o fsm.o
fsm_test_DEPENDENCIES =  recyclebin.o
fsm_test_LDFLAGS = 
matrix_test_OBJECTS =  matrix.test.o matrix.o
matrix_test_LDADD = $(LDADD)
matrix_test_DEPENDENCIES = 
matrix_test_LDFLAGS = 
pqueue_test_OBJECTS =  pqueue.test.o pqueue.o recyclebin.o
pqueue_test_LDADD = $(LDADD)
pqueue_test_DEPENDENCIES = 
pqueue_test_LDFLAGS = 
slist_test_OBJECTS =  slist.test.o slist.o recyclebin.o
slist_test_LDADD = $(LDADD)
slist_test_DEPENDENCIES = 
slist_test_LDFLAGS = 
rangetree_test_OBJECTS =  rangetree.test.o rangetree.o
rangetree_test_DEPENDENCIES =  recyclebin.o
rangetree_test_LDFLAGS = 
vfsm_test_OBJECTS =  vfsm.test.o vfsm.o
vfsm_test_DEPENDENCIES = 
vfsm_test_LDFLAGS = 
recyclebin_test_OBJECTS =  recyclebin.test.o recyclebin.o
recyclebin_test_LDADD = $(LDADD)
recyclebin_test_DEPENDENCIES = 
recyclebin_test_LDFLAGS = 
sparsecache_test_OBJECTS =  sparsecache.test.o sparsecache.o
sparsecache_test_LDADD = $(LDADD)
sparsecache_test_DEPENDENCIES = 
sparsecache_test_LDFLAGS = 
dejavu_test_OBJECTS =  dejavu.test.o dejavu.o
dejavu_test_LDADD = $(LDADD)
dejavu_test_DEPENDENCIES = 
dejavu_test_LDFLAGS = 
bitarray_test_OBJECTS =  bitarray.test.o bitarray.o
bitarray_test_LDADD = $(LDADD)
bitarray_test_DEPENDENCIES = 
bitarray_test_LDFLAGS = 
splaytree_test_OBJECTS =  splaytree.test.o splaytree.o
splaytree_test_DEPENDENCIES =  recyclebin.o
splaytree_test_LDFLAGS = 
noitree_test_OBJECTS =  noitree.test.o noitree.o
noitree_test_DEPENDENCIES =  recyclebin.o splaytree.o
noitree_test_LDFLAGS = 
CFLAGS = -g -O2 -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include   -DG_DISABLE_ASSERT -O3 -finline-functions -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DUSE_PTHREADS
COMPILE = $(CC) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)
CCLD = $(CC)
LINK = $(CCLD) $(AM_CFLAGS) $(CFLAGS) $(LDFLAGS) -o $@
HEADERS =  $(noinst_HEADERS)

DIST_COMMON =  Makefile.am Makefile.in


DISTFILES = $(DIST_COMMON) $(SOURCES) $(HEADERS) $(TEXINFOS) $(EXTRA_DIST)

TAR = tar
GZIP_ENV = --best
SOURCES = $(fsm_test_SOURCES) $(matrix_test_SOURCES) $(pqueue_test_SOURCES) $(slist_test_SOURCES) $(rangetree_test_SOURCES) $(vfsm_test_SOURCES) $(recyclebin_test_SOURCES) $(sparsecache_test_SOURCES) $(dejavu_test_SOURCES) $(bitarray_test_SOURCES) $(splaytree_test_SOURCES) $(noitree_test_SOURCES)
OBJECTS = $(fsm_test_OBJECTS) $(matrix_test_OBJECTS) $(pqueue_test_OBJECTS) $(slist_test_OBJECTS) $(rangetree_test_OBJECTS) $(vfsm_test_OBJECTS) $(recyclebin_test_OBJECTS) $(sparsecache_test_OBJECTS) $(dejavu_test_OBJECTS) $(bitarray_test_OBJECTS) $(splaytree_test_OBJECTS) $(noitree_test_OBJECTS)

all: all-redirect
.SUFFIXES:
.SUFFIXES: .S .c .o .s
$(srcdir)/Makefile.in: Makefile.am $(top_srcdir)/configure.in $(ACLOCAL_M4) 
	cd $(top_srcdir) && $(AUTOMAKE) --gnu --include-deps src/struct/Makefile

Makefile: $(srcdir)/Makefile.in  $(top_builddir)/config.status
	cd $(top_builddir) \
	  && CONFIG_FILES=$(subdir)/$@ CONFIG_HEADERS= $(SHELL) ./config.status


mostlyclean-noinstPROGRAMS:

clean-noinstPROGRAMS:
	-test -z "$(noinst_PROGRAMS)" || rm -f $(noinst_PROGRAMS)

distclean-noinstPROGRAMS:

maintainer-clean-noinstPROGRAMS:

.c.o:
	$(COMPILE) -c $<

.s.o:
	$(COMPILE) -c $<

.S.o:
	$(COMPILE) -c $<

mostlyclean-compile:
	-rm -f *.o core *.core

clean-compile:

distclean-compile:
	-rm -f *.tab.c

maintainer-clean-compile:

fsm.test: $(fsm_test_OBJECTS) $(fsm_test_DEPENDENCIES)
	@rm -f fsm.test
	$(LINK) $(fsm_test_LDFLAGS) $(fsm_test_OBJECTS) $(fsm_test_LDADD) $(LIBS)

matrix.test: $(matrix_test_OBJECTS) $(matrix_test_DEPENDENCIES)
	@rm -f matrix.test
	$(LINK) $(matrix_test_LDFLAGS) $(matrix_test_OBJECTS) $(matrix_test_LDADD) $(LIBS)

pqueue.test: $(pqueue_test_OBJECTS) $(pqueue_test_DEPENDENCIES)
	@rm -f pqueue.test
	$(LINK) $(pqueue_test_LDFLAGS) $(pqueue_test_OBJECTS) $(pqueue_test_LDADD) $(LIBS)

slist.test: $(slist_test_OBJECTS) $(slist_test_DEPENDENCIES)
	@rm -f slist.test
	$(LINK) $(slist_test_LDFLAGS) $(slist_test_OBJECTS) $(slist_test_LDADD) $(LIBS)

rangetree.test: $(rangetree_test_OBJECTS) $(rangetree_test_DEPENDENCIES)
	@rm -f rangetree.test
	$(LINK) $(rangetree_test_LDFLAGS) $(rangetree_test_OBJECTS) $(rangetree_test_LDADD) $(LIBS)

vfsm.test: $(vfsm_test_OBJECTS) $(vfsm_test_DEPENDENCIES)
	@rm -f vfsm.test
	$(LINK) $(vfsm_test_LDFLAGS) $(vfsm_test_OBJECTS) $(vfsm_test_LDADD) $(LIBS)

recyclebin.test: $(recyclebin_test_OBJECTS) $(recyclebin_test_DEPENDENCIES)
	@rm -f recyclebin.test
	$(LINK) $(recyclebin_test_LDFLAGS) $(recyclebin_test_OBJECTS) $(recyclebin_test_LDADD) $(LIBS)

sparsecache.test: $(sparsecache_test_OBJECTS) $(sparsecache_test_DEPENDENCIES)
	@rm -f sparsecache.test
	$(LINK) $(sparsecache_test_LDFLAGS) $(sparsecache_test_OBJECTS) $(sparsecache_test_LDADD) $(LIBS)

dejavu.test: $(dejavu_test_OBJECTS) $(dejavu_test_DEPENDENCIES)
	@rm -f dejavu.test
	$(LINK) $(dejavu_test_LDFLAGS) $(dejavu_test_OBJECTS) $(dejavu_test_LDADD) $(LIBS)

bitarray.test: $(bitarray_test_OBJECTS) $(bitarray_test_DEPENDENCIES)
	@rm -f bitarray.test
	$(LINK) $(bitarray_test_LDFLAGS) $(bitarray_test_OBJECTS) $(bitarray_test_LDADD) $(LIBS)

splaytree.test: $(splaytree_test_OBJECTS) $(splaytree_test_DEPENDENCIES)
	@rm -f splaytree.test
	$(LINK) $(splaytree_test_LDFLAGS) $(splaytree_test_OBJECTS) $(splaytree_test_LDADD) $(LIBS)

noitree.test: $(noitree_test_OBJECTS) $(noitree_test_DEPENDENCIES)
	@rm -f noitree.test
	$(LINK) $(noitree_test_LDFLAGS) $(noitree_test_OBJECTS) $(noitree_test_LDADD) $(LIBS)

tags: TAGS

ID: $(HEADERS) $(SOURCES) $(LISP)
	list='$(SOURCES) $(HEADERS)'; \
	unique=`for i in $$list; do echo $$i; done | \
	  awk '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	here=`pwd` && cd $(srcdir) \
	  && mkid -f$$here/ID $$unique $(LISP)

TAGS:  $(HEADERS) $(SOURCES)  $(TAGS_DEPENDENCIES) $(LISP)
	tags=; \
	here=`pwd`; \
	list='$(SOURCES) $(HEADERS)'; \
	unique=`for i in $$list; do echo $$i; done | \
	  awk '    { files[$$0] = 1; } \
	       END { for (i in files) print i; }'`; \
	test -z "$(ETAGS_ARGS)$$unique$(LISP)$$tags" \
	  || (cd $(srcdir) && etags -o $$here/TAGS $(ETAGS_ARGS) $$tags  $$unique $(LISP))

mostlyclean-tags:

clean-tags:

distclean-tags:
	-rm -f TAGS ID

maintainer-clean-tags:

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)

subdir = src/struct

distdir: $(DISTFILES)
	@for file in $(DISTFILES); do \
	  d=$(srcdir); \
	  if test -d $$d/$$file; then \
	    cp -pr $$d/$$file $(distdir)/$$file; \
	  else \
	    test -f $(distdir)/$$file \
	    || ln $$d/$$file $(distdir)/$$file 2> /dev/null \
	    || cp -p $$d/$$file $(distdir)/$$file || :; \
	  fi; \
	done
check-TESTS: $(TESTS)
	@failed=0; all=0; \
	srcdir=$(srcdir); export srcdir; \
	for tst in $(TESTS); do \
	  if test -f $$tst; then dir=.; \
	  else dir="$(srcdir)"; fi; \
	  if $(TESTS_ENVIRONMENT) $$dir/$$tst; then \
	    all=`expr $$all + 1`; \
	    echo "PASS: $$tst"; \
	  elif test $$? -ne 77; then \
	    all=`expr $$all + 1`; \
	    failed=`expr $$failed + 1`; \
	    echo "FAIL: $$tst"; \
	  fi; \
	done; \
	if test "$$failed" -eq 0; then \
	  banner="All $$all tests passed"; \
	else \
	  banner="$$failed of $$all tests failed"; \
	fi; \
	dashes=`echo "$$banner" | sed s/./=/g`; \
	echo "$$dashes"; \
	echo "$$banner"; \
	echo "$$dashes"; \
	test "$$failed" -eq 0
info-am:
info: info-am
dvi-am:
dvi: dvi-am
check-am: all-am
	$(MAKE) $(AM_MAKEFLAGS) check-TESTS
check: check-am
installcheck-am:
installcheck: installcheck-am
install-exec-am:
install-exec: install-exec-am

install-data-am:
install-data: install-data-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am
install: install-am
uninstall-am:
uninstall: uninstall-am
all-am: Makefile $(PROGRAMS) $(HEADERS)
all-redirect: all-am
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) AM_INSTALL_PROGRAM_FLAGS=-s install
installdirs:


mostlyclean-generic:

clean-generic:

distclean-generic:
	-rm -f Makefile $(CONFIG_CLEAN_FILES)
	-rm -f config.cache config.log stamp-h stamp-h[0-9]*

maintainer-clean-generic:
	-test -z "$(MAINTAINERCLEANFILES)" || rm -f $(MAINTAINERCLEANFILES)
mostlyclean-am:  mostlyclean-noinstPROGRAMS mostlyclean-compile \
		mostlyclean-tags mostlyclean-generic

mostlyclean: mostlyclean-am

clean-am:  clean-noinstPROGRAMS clean-compile clean-tags clean-generic \
		mostlyclean-am

clean: clean-am

distclean-am:  distclean-noinstPROGRAMS distclean-compile distclean-tags \
		distclean-generic clean-am

distclean: distclean-am

maintainer-clean-am:  maintainer-clean-noinstPROGRAMS \
		maintainer-clean-compile maintainer-clean-tags \
		maintainer-clean-generic distclean-am
	@echo "This command is intended for maintainers to use;"
	@echo "it deletes files that may require special tools to rebuild."

maintainer-clean: maintainer-clean-am

.PHONY: mostlyclean-noinstPROGRAMS distclean-noinstPROGRAMS \
clean-noinstPROGRAMS maintainer-clean-noinstPROGRAMS \
mostlyclean-compile distclean-compile clean-compile \
maintainer-clean-compile tags mostlyclean-tags distclean-tags \
clean-tags maintainer-clean-tags distdir check-TESTS info-am info \
dvi-am dvi check check-am installcheck-am installcheck install-exec-am \
install-exec install-data-am install-data install-am install \
uninstall-am uninstall all-redirect all-am all installdirs \
mostlyclean-generic distclean-generic clean-generic \
maintainer-clean-generic clean mostlyclean distclean maintainer-clean


# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
