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

TESTS = region.test c4.test alignment.test codegen.test optimal.test          layout.test viterbi.test opair.test subopt.test cgutil.test


noinst_PROGRAMS = $(TESTS)

INCLUDES = -I$(top_srcdir)/src/struct                           -I$(top_srcdir)/src/comparison                       -I$(top_srcdir)/src/sequence                         -I$(top_srcdir)/src/general                          -DSOURCE_ROOT_DIR="\"/loginhome/jzucker/src/jamg-code/3rd_party/exonerate\""            -DGLIB_CFLAGS="\"-I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include  \""                    -DHOSTTYPE="\"x86_64-unknown-linux-gnu\""


noinst_HEADERS = region.h c4.h alignment.h codegen.h optimal.h                   layout.h viterbi.h opair.h subopt.h                  cgutil.h


region_test_SOURCES = region.test.c region.c

c4_test_SOURCES = c4.test.c c4.c region.c
c4_test_LDADD = $(top_srcdir)/src/struct/matrix.o

SEQUENCE_OBJ = $(top_srcdir)/src/struct/sparsecache.o                 $(top_srcdir)/src/struct/matrix.o                      $(top_srcdir)/src/sequence/sequence.o                  $(top_srcdir)/src/sequence/alphabet.o                  $(top_srcdir)/src/sequence/splice.o                    $(top_srcdir)/src/sequence/translate.o                 $(top_srcdir)/src/general/argument.o                   $(top_srcdir)/src/general/lineparse.o                  -lm


ALIGNMENT_OBJ = $(top_srcdir)/src/comparison/match.o                     $(top_srcdir)/src/sequence/codonsubmat.o                 $(top_srcdir)/src/sequence/submat.o                      $(SEQUENCE_OBJ)

alignment_test_SOURCES = alignment.test.c alignment.c c4.c region.c
alignment_test_LDADD = $(ALIGNMENT_OBJ)

codegen_test_SOURCES = codegen.test.c codegen.c
codegen_test_LDADD = $(top_srcdir)/src/general/argument.o

cgutil_test_SOURCES = cgutil.test.c cgutil.c codegen.c c4.c region.c
cgutil_test_LDADD = $(top_srcdir)/src/general/argument.o                     $(top_srcdir)/src/struct/matrix.o


viterbi_test_SOURCES = viterbi.test.c viterbi.c c4.c region.c                        alignment.c codegen.c layout.c subopt.c                        cgutil.c

viterbi_test_LDADD = $(top_srcdir)/src/struct/slist.o                           $(top_srcdir)/src/struct/recyclebin.o                      $(top_srcdir)/src/struct/rangetree.o                       $(ALIGNMENT_OBJ)


optimal_test_SOURCES = optimal.test.c optimal.c c4.c alignment.c                        codegen.c region.c layout.c viterbi.c                            subopt.c cgutil.c

optimal_test_LDADD = $(top_srcdir)/src/struct/slist.o                           $(top_srcdir)/src/struct/recyclebin.o                      $(top_srcdir)/src/struct/rangetree.o                       $(ALIGNMENT_OBJ)


layout_test_SOURCES = layout.test.c layout.c
layout_test_LDADD = $(top_srcdir)/src/struct/matrix.o

opair_test_SOURCES = opair.test.c opair.c optimal.c c4.c alignment.c                      codegen.c region.c layout.c viterbi.c subopt.c                       cgutil.c

opair_test_LDADD = $(top_srcdir)/src/struct/slist.o                         $(top_srcdir)/src/struct/recyclebin.o                    $(top_srcdir)/src/struct/rangetree.o                     $(ALIGNMENT_OBJ)


subopt_test_SOURCES = subopt.test.c subopt.c region.c
subopt_test_LDADD = $(top_srcdir)/src/struct/rangetree.o                     $(top_srcdir)/src/struct/recyclebin.o


# Files to clear away

MAINTAINERCLEANFILES = Makefile.in
mkinstalldirs = $(SHELL) $(top_srcdir)/mkinstalldirs
CONFIG_CLEAN_FILES = 
PROGRAMS =  $(noinst_PROGRAMS)


DEFS = -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE=\"exonerate\" -DVERSION=\"2.2.0\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIMITS_H=1 -DHAVE_ERRNO_H=1 -DHAVE_MATH_H=1 -DHAVE_STRINGS_H=1 -DHAVE_CTYPE_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_TIME_H=1 -DHAVE_UNISTD_H=1 -I. -I$(srcdir) 
CPPFLAGS = 
LDFLAGS =  -lpthread
LIBS =  -lglib-2.0  
region_test_OBJECTS =  region.test.o region.o
region_test_LDADD = $(LDADD)
region_test_DEPENDENCIES = 
region_test_LDFLAGS = 
c4_test_OBJECTS =  c4.test.o c4.o region.o
c4_test_DEPENDENCIES =  $(top_srcdir)/src/struct/matrix.o
c4_test_LDFLAGS = 
alignment_test_OBJECTS =  alignment.test.o alignment.o c4.o region.o
alignment_test_DEPENDENCIES =  $(top_srcdir)/src/comparison/match.o \
$(top_srcdir)/src/sequence/codonsubmat.o \
$(top_srcdir)/src/sequence/submat.o \
$(top_srcdir)/src/struct/sparsecache.o \
$(top_srcdir)/src/struct/matrix.o $(top_srcdir)/src/sequence/sequence.o \
$(top_srcdir)/src/sequence/alphabet.o \
$(top_srcdir)/src/sequence/splice.o \
$(top_srcdir)/src/sequence/translate.o \
$(top_srcdir)/src/general/argument.o \
$(top_srcdir)/src/general/lineparse.o
alignment_test_LDFLAGS = 
codegen_test_OBJECTS =  codegen.test.o codegen.o
codegen_test_DEPENDENCIES =  $(top_srcdir)/src/general/argument.o
codegen_test_LDFLAGS = 
optimal_test_OBJECTS =  optimal.test.o optimal.o c4.o alignment.o \
codegen.o region.o layout.o viterbi.o subopt.o cgutil.o
optimal_test_DEPENDENCIES =  $(top_srcdir)/src/struct/slist.o \
$(top_srcdir)/src/struct/recyclebin.o \
$(top_srcdir)/src/struct/rangetree.o \
$(top_srcdir)/src/comparison/match.o \
$(top_srcdir)/src/sequence/codonsubmat.o \
$(top_srcdir)/src/sequence/submat.o \
$(top_srcdir)/src/struct/sparsecache.o \
$(top_srcdir)/src/struct/matrix.o $(top_srcdir)/src/sequence/sequence.o \
$(top_srcdir)/src/sequence/alphabet.o \
$(top_srcdir)/src/sequence/splice.o \
$(top_srcdir)/src/sequence/translate.o \
$(top_srcdir)/src/general/argument.o \
$(top_srcdir)/src/general/lineparse.o
optimal_test_LDFLAGS = 
layout_test_OBJECTS =  layout.test.o layout.o
layout_test_DEPENDENCIES =  $(top_srcdir)/src/struct/matrix.o
layout_test_LDFLAGS = 
viterbi_test_OBJECTS =  viterbi.test.o viterbi.o c4.o region.o \
alignment.o codegen.o layout.o subopt.o cgutil.o
viterbi_test_DEPENDENCIES =  $(top_srcdir)/src/struct/slist.o \
$(top_srcdir)/src/struct/recyclebin.o \
$(top_srcdir)/src/struct/rangetree.o \
$(top_srcdir)/src/comparison/match.o \
$(top_srcdir)/src/sequence/codonsubmat.o \
$(top_srcdir)/src/sequence/submat.o \
$(top_srcdir)/src/struct/sparsecache.o \
$(top_srcdir)/src/struct/matrix.o $(top_srcdir)/src/sequence/sequence.o \
$(top_srcdir)/src/sequence/alphabet.o \
$(top_srcdir)/src/sequence/splice.o \
$(top_srcdir)/src/sequence/translate.o \
$(top_srcdir)/src/general/argument.o \
$(top_srcdir)/src/general/lineparse.o
viterbi_test_LDFLAGS = 
opair_test_OBJECTS =  opair.test.o opair.o optimal.o c4.o alignment.o \
codegen.o region.o layout.o viterbi.o subopt.o cgutil.o
opair_test_DEPENDENCIES =  $(top_srcdir)/src/struct/slist.o \
$(top_srcdir)/src/struct/recyclebin.o \
$(top_srcdir)/src/struct/rangetree.o \
$(top_srcdir)/src/comparison/match.o \
$(top_srcdir)/src/sequence/codonsubmat.o \
$(top_srcdir)/src/sequence/submat.o \
$(top_srcdir)/src/struct/sparsecache.o \
$(top_srcdir)/src/struct/matrix.o $(top_srcdir)/src/sequence/sequence.o \
$(top_srcdir)/src/sequence/alphabet.o \
$(top_srcdir)/src/sequence/splice.o \
$(top_srcdir)/src/sequence/translate.o \
$(top_srcdir)/src/general/argument.o \
$(top_srcdir)/src/general/lineparse.o
opair_test_LDFLAGS = 
subopt_test_OBJECTS =  subopt.test.o subopt.o region.o
subopt_test_DEPENDENCIES =  $(top_srcdir)/src/struct/rangetree.o \
$(top_srcdir)/src/struct/recyclebin.o
subopt_test_LDFLAGS = 
cgutil_test_OBJECTS =  cgutil.test.o cgutil.o codegen.o c4.o region.o
cgutil_test_DEPENDENCIES =  $(top_srcdir)/src/general/argument.o \
$(top_srcdir)/src/struct/matrix.o
cgutil_test_LDFLAGS = 
CFLAGS = -g -O2 -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include   -DG_DISABLE_ASSERT -O3 -finline-functions -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DUSE_PTHREADS
COMPILE = $(CC) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)
CCLD = $(CC)
LINK = $(CCLD) $(AM_CFLAGS) $(CFLAGS) $(LDFLAGS) -o $@
HEADERS =  $(noinst_HEADERS)

DIST_COMMON =  Makefile.am Makefile.in


DISTFILES = $(DIST_COMMON) $(SOURCES) $(HEADERS) $(TEXINFOS) $(EXTRA_DIST)

TAR = tar
GZIP_ENV = --best
SOURCES = $(region_test_SOURCES) $(c4_test_SOURCES) $(alignment_test_SOURCES) $(codegen_test_SOURCES) $(optimal_test_SOURCES) $(layout_test_SOURCES) $(viterbi_test_SOURCES) $(opair_test_SOURCES) $(subopt_test_SOURCES) $(cgutil_test_SOURCES)
OBJECTS = $(region_test_OBJECTS) $(c4_test_OBJECTS) $(alignment_test_OBJECTS) $(codegen_test_OBJECTS) $(optimal_test_OBJECTS) $(layout_test_OBJECTS) $(viterbi_test_OBJECTS) $(opair_test_OBJECTS) $(subopt_test_OBJECTS) $(cgutil_test_OBJECTS)

all: all-redirect
.SUFFIXES:
.SUFFIXES: .S .c .o .s
$(srcdir)/Makefile.in: Makefile.am $(top_srcdir)/configure.in $(ACLOCAL_M4) 
	cd $(top_srcdir) && $(AUTOMAKE) --gnu --include-deps src/c4/Makefile

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

region.test: $(region_test_OBJECTS) $(region_test_DEPENDENCIES)
	@rm -f region.test
	$(LINK) $(region_test_LDFLAGS) $(region_test_OBJECTS) $(region_test_LDADD) $(LIBS)

c4.test: $(c4_test_OBJECTS) $(c4_test_DEPENDENCIES)
	@rm -f c4.test
	$(LINK) $(c4_test_LDFLAGS) $(c4_test_OBJECTS) $(c4_test_LDADD) $(LIBS)

alignment.test: $(alignment_test_OBJECTS) $(alignment_test_DEPENDENCIES)
	@rm -f alignment.test
	$(LINK) $(alignment_test_LDFLAGS) $(alignment_test_OBJECTS) $(alignment_test_LDADD) $(LIBS)

codegen.test: $(codegen_test_OBJECTS) $(codegen_test_DEPENDENCIES)
	@rm -f codegen.test
	$(LINK) $(codegen_test_LDFLAGS) $(codegen_test_OBJECTS) $(codegen_test_LDADD) $(LIBS)

optimal.test: $(optimal_test_OBJECTS) $(optimal_test_DEPENDENCIES)
	@rm -f optimal.test
	$(LINK) $(optimal_test_LDFLAGS) $(optimal_test_OBJECTS) $(optimal_test_LDADD) $(LIBS)

layout.test: $(layout_test_OBJECTS) $(layout_test_DEPENDENCIES)
	@rm -f layout.test
	$(LINK) $(layout_test_LDFLAGS) $(layout_test_OBJECTS) $(layout_test_LDADD) $(LIBS)

viterbi.test: $(viterbi_test_OBJECTS) $(viterbi_test_DEPENDENCIES)
	@rm -f viterbi.test
	$(LINK) $(viterbi_test_LDFLAGS) $(viterbi_test_OBJECTS) $(viterbi_test_LDADD) $(LIBS)

opair.test: $(opair_test_OBJECTS) $(opair_test_DEPENDENCIES)
	@rm -f opair.test
	$(LINK) $(opair_test_LDFLAGS) $(opair_test_OBJECTS) $(opair_test_LDADD) $(LIBS)

subopt.test: $(subopt_test_OBJECTS) $(subopt_test_DEPENDENCIES)
	@rm -f subopt.test
	$(LINK) $(subopt_test_LDFLAGS) $(subopt_test_OBJECTS) $(subopt_test_LDADD) $(LIBS)

cgutil.test: $(cgutil_test_OBJECTS) $(cgutil_test_DEPENDENCIES)
	@rm -f cgutil.test
	$(LINK) $(cgutil_test_LDFLAGS) $(cgutil_test_OBJECTS) $(cgutil_test_LDADD) $(LIBS)

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

subdir = src/c4

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
