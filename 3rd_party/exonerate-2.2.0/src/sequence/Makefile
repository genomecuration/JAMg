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

TESTS = sequence.test submat.test translate.test splice.test         alphabet.test codonsubmat.test

noinst_PROGRAMS = $(TESTS)

INCLUDES = -I$(top_srcdir)/src/struct            -I$(top_srcdir)/src/general


noinst_HEADERS = sequence.h submat.h translate.h splice.h alphabet.h                  codonsubmat.h


splice_test_SOURCES = splice.test.c splice.c
SPLICE_OBJ = $(top_srcdir)/src/struct/matrix.o                   $(top_srcdir)/src/struct/sparsecache.o              $(top_srcdir)/src/general/argument.o                $(top_srcdir)/src/general/lineparse.o               -lm

splice_test_LDADD = $(SPLICE_OBJ)

sequence_test_SOURCES = sequence.test.c sequence.c alphabet.c translate.c                         splice.c

sequence_test_LDADD = $(SPLICE_OBJ)
submat_test_SOURCES = submat.test.c submat.c

translate_test_SOURCES = translate.test.c translate.c
translate_test_LDADD = $(top_srcdir)/src/general/argument.o

alphabet_test_SOURCES = alphabet.test.c alphabet.c
alphabet_test_LDADD = $(top_srcdir)/src/general/argument.o

codonsubmat_test_SOURCES = codonsubmat.test.c codonsubmat.c
codonsubmat_test_LDADD = $(top_srcdir)/src/general/argument.o -lm

# Files to clear away

MAINTAINERCLEANFILES = Makefile.in
mkinstalldirs = $(SHELL) $(top_srcdir)/mkinstalldirs
CONFIG_CLEAN_FILES = 
PROGRAMS =  $(noinst_PROGRAMS)


DEFS = -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE=\"exonerate\" -DVERSION=\"2.2.0\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIMITS_H=1 -DHAVE_ERRNO_H=1 -DHAVE_MATH_H=1 -DHAVE_STRINGS_H=1 -DHAVE_CTYPE_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_TIME_H=1 -DHAVE_UNISTD_H=1 -I. -I$(srcdir) 
CPPFLAGS = 
LDFLAGS =  -lpthread
LIBS =  -lglib-2.0  
sequence_test_OBJECTS =  sequence.test.o sequence.o alphabet.o \
translate.o splice.o
sequence_test_DEPENDENCIES =  $(top_srcdir)/src/struct/matrix.o \
$(top_srcdir)/src/struct/sparsecache.o \
$(top_srcdir)/src/general/argument.o \
$(top_srcdir)/src/general/lineparse.o
sequence_test_LDFLAGS = 
submat_test_OBJECTS =  submat.test.o submat.o
submat_test_LDADD = $(LDADD)
submat_test_DEPENDENCIES = 
submat_test_LDFLAGS = 
translate_test_OBJECTS =  translate.test.o translate.o
translate_test_DEPENDENCIES =  $(top_srcdir)/src/general/argument.o
translate_test_LDFLAGS = 
splice_test_OBJECTS =  splice.test.o splice.o
splice_test_DEPENDENCIES =  $(top_srcdir)/src/struct/matrix.o \
$(top_srcdir)/src/struct/sparsecache.o \
$(top_srcdir)/src/general/argument.o \
$(top_srcdir)/src/general/lineparse.o
splice_test_LDFLAGS = 
alphabet_test_OBJECTS =  alphabet.test.o alphabet.o
alphabet_test_DEPENDENCIES =  $(top_srcdir)/src/general/argument.o
alphabet_test_LDFLAGS = 
codonsubmat_test_OBJECTS =  codonsubmat.test.o codonsubmat.o
codonsubmat_test_DEPENDENCIES =  $(top_srcdir)/src/general/argument.o
codonsubmat_test_LDFLAGS = 
CFLAGS = -g -O2 -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include   -DG_DISABLE_ASSERT -O3 -finline-functions -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DUSE_PTHREADS
COMPILE = $(CC) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)
CCLD = $(CC)
LINK = $(CCLD) $(AM_CFLAGS) $(CFLAGS) $(LDFLAGS) -o $@
HEADERS =  $(noinst_HEADERS)

DIST_COMMON =  Makefile.am Makefile.in


DISTFILES = $(DIST_COMMON) $(SOURCES) $(HEADERS) $(TEXINFOS) $(EXTRA_DIST)

TAR = tar
GZIP_ENV = --best
SOURCES = $(sequence_test_SOURCES) $(submat_test_SOURCES) $(translate_test_SOURCES) $(splice_test_SOURCES) $(alphabet_test_SOURCES) $(codonsubmat_test_SOURCES)
OBJECTS = $(sequence_test_OBJECTS) $(submat_test_OBJECTS) $(translate_test_OBJECTS) $(splice_test_OBJECTS) $(alphabet_test_OBJECTS) $(codonsubmat_test_OBJECTS)

all: all-redirect
.SUFFIXES:
.SUFFIXES: .S .c .o .s
$(srcdir)/Makefile.in: Makefile.am $(top_srcdir)/configure.in $(ACLOCAL_M4) 
	cd $(top_srcdir) && $(AUTOMAKE) --gnu --include-deps src/sequence/Makefile

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

sequence.test: $(sequence_test_OBJECTS) $(sequence_test_DEPENDENCIES)
	@rm -f sequence.test
	$(LINK) $(sequence_test_LDFLAGS) $(sequence_test_OBJECTS) $(sequence_test_LDADD) $(LIBS)

submat.test: $(submat_test_OBJECTS) $(submat_test_DEPENDENCIES)
	@rm -f submat.test
	$(LINK) $(submat_test_LDFLAGS) $(submat_test_OBJECTS) $(submat_test_LDADD) $(LIBS)

translate.test: $(translate_test_OBJECTS) $(translate_test_DEPENDENCIES)
	@rm -f translate.test
	$(LINK) $(translate_test_LDFLAGS) $(translate_test_OBJECTS) $(translate_test_LDADD) $(LIBS)

splice.test: $(splice_test_OBJECTS) $(splice_test_DEPENDENCIES)
	@rm -f splice.test
	$(LINK) $(splice_test_LDFLAGS) $(splice_test_OBJECTS) $(splice_test_LDADD) $(LIBS)

alphabet.test: $(alphabet_test_OBJECTS) $(alphabet_test_DEPENDENCIES)
	@rm -f alphabet.test
	$(LINK) $(alphabet_test_LDFLAGS) $(alphabet_test_OBJECTS) $(alphabet_test_LDADD) $(LIBS)

codonsubmat.test: $(codonsubmat_test_OBJECTS) $(codonsubmat_test_DEPENDENCIES)
	@rm -f codonsubmat.test
	$(LINK) $(codonsubmat_test_LDFLAGS) $(codonsubmat_test_OBJECTS) $(codonsubmat_test_LDADD) $(LIBS)

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

subdir = src/sequence

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
