# Makefile.in generated automatically by automake 1.4-p6 from Makefile.am

# Copyright (C) 1994, 1995-8, 1999, 2001 Free Software Foundation, Inc.
# This Makefile.in is free software; the Free Software Foundation
# gives unlimited permission to copy and/or distribute it,
# with or without modifications, as long as this notice is preserved.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY, to the extent permitted by law; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A
# PARTICULAR PURPOSE.


SHELL = /bin/sh

srcdir = .
top_srcdir = ../..

prefix = /home/pap056/workspace/transcripts4community/jamg/3rd_party/exonerate-2.2.0/../
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

ACLOCAL = /home/pap056/workspace/transcripts4community/jamg/3rd_party/exonerate-2.2.0/missing aclocal-1.4
AUTOCONF = autoconf
AUTOMAKE = /home/pap056/workspace/transcripts4community/jamg/3rd_party/exonerate-2.2.0/missing automake-1.4
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
MAKEINFO = /home/pap056/workspace/transcripts4community/jamg/3rd_party/exonerate-2.2.0/missing makeinfo
PACKAGE = exonerate
PKG_CONFIG = /usr/bin/pkg-config
VERSION = 2.2.0
codegen_extra_ldadd = viterbi.o scheduler.o c4_model_archive.a
codegen_extra_sources = c4_model_archive.h
datarootdir = ${prefix}/share
glib_cflags = -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include  
glib_libs = -lglib-2.0  
host = x86_64-unknown-linux-gnu
installed_util_list =                                                fastaclean fastaclip fastachecksum fastacomposition            fastadiff fastaexplode fastafetch fastahardmask                fastaindex fastalength fastanrdb fastaoverlap                  fastareformat fastaremove fastarevcomp fastasoftmask           fastasort fastasplit fastasubseq fastatranslate                fastavalidcds fasta2esd fastaannotatecdna esd2esi
source_root_dir = /home/pap056/workspace/transcripts4community/jamg/3rd_party/exonerate-2.2.0
use_pthreads = 

TESTS = gam.test analysis.test bsam.test
noinst_PROGRAMS = $(TESTS)

INCLUDES = -I$(top_srcdir)/src/sequence              -I$(top_srcdir)/src/database              -I$(top_srcdir)/src/comparison            -I$(top_srcdir)/src/struct                -I$(top_srcdir)/src/c4                    -I$(top_srcdir)/src/bsdp                  -I$(top_srcdir)/src/sdp                   -I$(top_srcdir)/src/model                 -I$(top_srcdir)/src/general               -DHOSTTYPE="\"x86_64-unknown-linux-gnu\""                   -DSOURCE_ROOT_DIR="\"/home/pap056/workspace/transcripts4community/jamg/3rd_party/exonerate-2.2.0\""


noinst_HEADERS = gam.h analysis.h bsam.h

C4_OBJECTS = $(top_srcdir)/src/c4/c4.o                             $(top_srcdir)/src/c4/codegen.o                        $(top_srcdir)/src/c4/cgutil.o                         $(top_srcdir)/src/c4/opair.o                          $(top_srcdir)/src/c4/alignment.o                      $(top_srcdir)/src/c4/optimal.o                        $(top_srcdir)/src/c4/viterbi.o                        $(top_srcdir)/src/c4/layout.o                         $(top_srcdir)/src/c4/region.o                         $(top_srcdir)/src/c4/subopt.o                         $(top_srcdir)/src/bsdp/bsdp.o                         $(top_srcdir)/src/bsdp/sar.o                          $(top_srcdir)/src/bsdp/heuristic.o                    $(top_srcdir)/src/bsdp/hpair.o                        $(top_srcdir)/src/sdp/boundary.o                      $(top_srcdir)/src/sdp/sdp.o                           $(top_srcdir)/src/sdp/scheduler.o                     $(top_srcdir)/src/sdp/straceback.o                    $(top_srcdir)/src/sdp/lookahead.o                     $(top_srcdir)/src/struct/matrix.o                     $(top_srcdir)/src/struct/slist.o                      $(top_srcdir)/src/struct/pqueue.o                     $(top_srcdir)/src/struct/recyclebin.o                 $(top_srcdir)/src/struct/rangetree.o                  $(top_srcdir)/src/sequence/sequence.o                 $(top_srcdir)/src/sequence/alphabet.o                 $(top_srcdir)/src/sequence/translate.o                $(top_srcdir)/src/sequence/splice.o                   $(top_srcdir)/src/model/affine.o                      $(top_srcdir)/src/model/est2genome.o                  $(top_srcdir)/src/model/ner.o                         $(top_srcdir)/src/model/protein2dna.o                 $(top_srcdir)/src/model/protein2genome.o              $(top_srcdir)/src/model/coding2coding.o               $(top_srcdir)/src/model/coding2genome.o               $(top_srcdir)/src/model/cdna2genome.o                 $(top_srcdir)/src/model/genome2genome.o               $(top_srcdir)/src/model/ungapped.o                    $(top_srcdir)/src/model/intron.o                      $(top_srcdir)/src/model/frameshift.o                  $(top_srcdir)/src/model/phase.o                       $(top_srcdir)/src/model/modeltype.o                   $(top_srcdir)/src/general/argument.o                  $(top_srcdir)/src/general/lineparse.o                 $(top_srcdir)/src/comparison/match.o                  $(top_srcdir)/src/sequence/submat.o                   $(top_srcdir)/src/sequence/codonsubmat.o              -lm


gam_test_SOURCES = gam.test.c gam.c
gam_test_LDADD = $(top_srcdir)/src/comparison/hspset.o                      $(top_srcdir)/src/comparison/comparison.o                  $(top_srcdir)/src/comparison/wordhood.o                    $(top_srcdir)/src/struct/sparsecache.o                     $(top_srcdir)/src/general/compoundfile.o                   $(C4_OBJECTS)


bsam_test_SOURCES = bsam.test.c bsam.c
bsam_test_LDADD = $(top_srcdir)/src/general/argument.o                         $(top_srcdir)/src/general/lineparse.o                        $(top_srcdir)/src/comparison/hspset.o                        $(top_srcdir)/src/comparison/comparison.o                    $(top_srcdir)/src/comparison/match.o                         $(top_srcdir)/src/comparison/wordhood.o                      $(top_srcdir)/src/sequence/sequence.o                        $(top_srcdir)/src/sequence/alphabet.o                        $(top_srcdir)/src/sequence/translate.o                       $(top_srcdir)/src/sequence/submat.o                          $(top_srcdir)/src/sequence/codonsubmat.o                     $(top_srcdir)/src/sequence/splice.o                          $(top_srcdir)/src/struct/dejavu.o                            $(top_srcdir)/src/struct/pqueue.o                            $(top_srcdir)/src/struct/sparsecache.o                       $(top_srcdir)/src/struct/recyclebin.o                        $(top_srcdir)/src/struct/matrix.o                            -lm

analysis_test_SOURCES = analysis.test.c analysis.c gam.c bsam.c
analysis_test_LDADD = $(top_srcdir)/src/comparison/hspset.o                           $(top_srcdir)/src/comparison/wordhood.o                         $(top_srcdir)/src/comparison/seeder.o                           $(top_srcdir)/src/comparison/comparison.o                       $(top_srcdir)/src/database/fastapipe.o                          $(top_srcdir)/src/database/fastadb.o                            $(top_srcdir)/src/struct/dejavu.o                               $(top_srcdir)/src/struct/fsm.o                                  $(top_srcdir)/src/struct/vfsm.o                                 $(top_srcdir)/src/struct/sparsecache.o                          $(top_srcdir)/src/general/compoundfile.o                        $(top_srcdir)/src/general/socket.o                              $(C4_OBJECTS)


# Files to clear away

MAINTAINERCLEANFILES = Makefile.in
mkinstalldirs = $(SHELL) $(top_srcdir)/mkinstalldirs
CONFIG_CLEAN_FILES = 
PROGRAMS =  $(noinst_PROGRAMS)


DEFS = -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE=\"exonerate\" -DVERSION=\"2.2.0\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIMITS_H=1 -DHAVE_ERRNO_H=1 -DHAVE_MATH_H=1 -DHAVE_STRINGS_H=1 -DHAVE_CTYPE_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_TIME_H=1 -DHAVE_UNISTD_H=1 -I. -I$(srcdir) 
CPPFLAGS =  -fPIC 
LDFLAGS =  -lpthread
LIBS =  -lglib-2.0  
gam_test_OBJECTS =  gam.test.o gam.o
gam_test_DEPENDENCIES =  $(top_srcdir)/src/comparison/hspset.o \
$(top_srcdir)/src/comparison/comparison.o \
$(top_srcdir)/src/comparison/wordhood.o \
$(top_srcdir)/src/struct/sparsecache.o \
$(top_srcdir)/src/general/compoundfile.o $(top_srcdir)/src/c4/c4.o \
$(top_srcdir)/src/c4/codegen.o $(top_srcdir)/src/c4/cgutil.o \
$(top_srcdir)/src/c4/opair.o $(top_srcdir)/src/c4/alignment.o \
$(top_srcdir)/src/c4/optimal.o $(top_srcdir)/src/c4/viterbi.o \
$(top_srcdir)/src/c4/layout.o $(top_srcdir)/src/c4/region.o \
$(top_srcdir)/src/c4/subopt.o $(top_srcdir)/src/bsdp/bsdp.o \
$(top_srcdir)/src/bsdp/sar.o $(top_srcdir)/src/bsdp/heuristic.o \
$(top_srcdir)/src/bsdp/hpair.o $(top_srcdir)/src/sdp/boundary.o \
$(top_srcdir)/src/sdp/sdp.o $(top_srcdir)/src/sdp/scheduler.o \
$(top_srcdir)/src/sdp/straceback.o $(top_srcdir)/src/sdp/lookahead.o \
$(top_srcdir)/src/struct/matrix.o $(top_srcdir)/src/struct/slist.o \
$(top_srcdir)/src/struct/pqueue.o $(top_srcdir)/src/struct/recyclebin.o \
$(top_srcdir)/src/struct/rangetree.o \
$(top_srcdir)/src/sequence/sequence.o \
$(top_srcdir)/src/sequence/alphabet.o \
$(top_srcdir)/src/sequence/translate.o \
$(top_srcdir)/src/sequence/splice.o $(top_srcdir)/src/model/affine.o \
$(top_srcdir)/src/model/est2genome.o $(top_srcdir)/src/model/ner.o \
$(top_srcdir)/src/model/protein2dna.o \
$(top_srcdir)/src/model/protein2genome.o \
$(top_srcdir)/src/model/coding2coding.o \
$(top_srcdir)/src/model/coding2genome.o \
$(top_srcdir)/src/model/cdna2genome.o \
$(top_srcdir)/src/model/genome2genome.o \
$(top_srcdir)/src/model/ungapped.o $(top_srcdir)/src/model/intron.o \
$(top_srcdir)/src/model/frameshift.o $(top_srcdir)/src/model/phase.o \
$(top_srcdir)/src/model/modeltype.o \
$(top_srcdir)/src/general/argument.o \
$(top_srcdir)/src/general/lineparse.o \
$(top_srcdir)/src/comparison/match.o \
$(top_srcdir)/src/sequence/submat.o \
$(top_srcdir)/src/sequence/codonsubmat.o
gam_test_LDFLAGS = 
analysis_test_OBJECTS =  analysis.test.o analysis.o gam.o bsam.o
analysis_test_DEPENDENCIES =  $(top_srcdir)/src/comparison/hspset.o \
$(top_srcdir)/src/comparison/wordhood.o \
$(top_srcdir)/src/comparison/seeder.o \
$(top_srcdir)/src/comparison/comparison.o \
$(top_srcdir)/src/database/fastapipe.o \
$(top_srcdir)/src/database/fastadb.o $(top_srcdir)/src/struct/dejavu.o \
$(top_srcdir)/src/struct/fsm.o $(top_srcdir)/src/struct/vfsm.o \
$(top_srcdir)/src/struct/sparsecache.o \
$(top_srcdir)/src/general/compoundfile.o \
$(top_srcdir)/src/general/socket.o $(top_srcdir)/src/c4/c4.o \
$(top_srcdir)/src/c4/codegen.o $(top_srcdir)/src/c4/cgutil.o \
$(top_srcdir)/src/c4/opair.o $(top_srcdir)/src/c4/alignment.o \
$(top_srcdir)/src/c4/optimal.o $(top_srcdir)/src/c4/viterbi.o \
$(top_srcdir)/src/c4/layout.o $(top_srcdir)/src/c4/region.o \
$(top_srcdir)/src/c4/subopt.o $(top_srcdir)/src/bsdp/bsdp.o \
$(top_srcdir)/src/bsdp/sar.o $(top_srcdir)/src/bsdp/heuristic.o \
$(top_srcdir)/src/bsdp/hpair.o $(top_srcdir)/src/sdp/boundary.o \
$(top_srcdir)/src/sdp/sdp.o $(top_srcdir)/src/sdp/scheduler.o \
$(top_srcdir)/src/sdp/straceback.o $(top_srcdir)/src/sdp/lookahead.o \
$(top_srcdir)/src/struct/matrix.o $(top_srcdir)/src/struct/slist.o \
$(top_srcdir)/src/struct/pqueue.o $(top_srcdir)/src/struct/recyclebin.o \
$(top_srcdir)/src/struct/rangetree.o \
$(top_srcdir)/src/sequence/sequence.o \
$(top_srcdir)/src/sequence/alphabet.o \
$(top_srcdir)/src/sequence/translate.o \
$(top_srcdir)/src/sequence/splice.o $(top_srcdir)/src/model/affine.o \
$(top_srcdir)/src/model/est2genome.o $(top_srcdir)/src/model/ner.o \
$(top_srcdir)/src/model/protein2dna.o \
$(top_srcdir)/src/model/protein2genome.o \
$(top_srcdir)/src/model/coding2coding.o \
$(top_srcdir)/src/model/coding2genome.o \
$(top_srcdir)/src/model/cdna2genome.o \
$(top_srcdir)/src/model/genome2genome.o \
$(top_srcdir)/src/model/ungapped.o $(top_srcdir)/src/model/intron.o \
$(top_srcdir)/src/model/frameshift.o $(top_srcdir)/src/model/phase.o \
$(top_srcdir)/src/model/modeltype.o \
$(top_srcdir)/src/general/argument.o \
$(top_srcdir)/src/general/lineparse.o \
$(top_srcdir)/src/comparison/match.o \
$(top_srcdir)/src/sequence/submat.o \
$(top_srcdir)/src/sequence/codonsubmat.o
analysis_test_LDFLAGS = 
bsam_test_OBJECTS =  bsam.test.o bsam.o
bsam_test_DEPENDENCIES =  $(top_srcdir)/src/general/argument.o \
$(top_srcdir)/src/general/lineparse.o \
$(top_srcdir)/src/comparison/hspset.o \
$(top_srcdir)/src/comparison/comparison.o \
$(top_srcdir)/src/comparison/match.o \
$(top_srcdir)/src/comparison/wordhood.o \
$(top_srcdir)/src/sequence/sequence.o \
$(top_srcdir)/src/sequence/alphabet.o \
$(top_srcdir)/src/sequence/translate.o \
$(top_srcdir)/src/sequence/submat.o \
$(top_srcdir)/src/sequence/codonsubmat.o \
$(top_srcdir)/src/sequence/splice.o $(top_srcdir)/src/struct/dejavu.o \
$(top_srcdir)/src/struct/pqueue.o \
$(top_srcdir)/src/struct/sparsecache.o \
$(top_srcdir)/src/struct/recyclebin.o $(top_srcdir)/src/struct/matrix.o
bsam_test_LDFLAGS = 
CFLAGS = -g -O2 -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include   -DG_DISABLE_ASSERT -O3 -finline-functions -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DUSE_PTHREADS
COMPILE = $(CC) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)
CCLD = $(CC)
LINK = $(CCLD) $(AM_CFLAGS) $(CFLAGS) $(LDFLAGS) -o $@
HEADERS =  $(noinst_HEADERS)

DIST_COMMON =  Makefile.am Makefile.in


DISTFILES = $(DIST_COMMON) $(SOURCES) $(HEADERS) $(TEXINFOS) $(EXTRA_DIST)

TAR = tar
GZIP_ENV = --best
SOURCES = $(gam_test_SOURCES) $(analysis_test_SOURCES) $(bsam_test_SOURCES)
OBJECTS = $(gam_test_OBJECTS) $(analysis_test_OBJECTS) $(bsam_test_OBJECTS)

all: all-redirect
.SUFFIXES:
.SUFFIXES: .S .c .o .s
$(srcdir)/Makefile.in: Makefile.am $(top_srcdir)/configure.in $(ACLOCAL_M4) 
	cd $(top_srcdir) && $(AUTOMAKE) --gnu --include-deps src/hub/Makefile

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

gam.test: $(gam_test_OBJECTS) $(gam_test_DEPENDENCIES)
	@rm -f gam.test
	$(LINK) $(gam_test_LDFLAGS) $(gam_test_OBJECTS) $(gam_test_LDADD) $(LIBS)

analysis.test: $(analysis_test_OBJECTS) $(analysis_test_DEPENDENCIES)
	@rm -f analysis.test
	$(LINK) $(analysis_test_LDFLAGS) $(analysis_test_OBJECTS) $(analysis_test_LDADD) $(LIBS)

bsam.test: $(bsam_test_OBJECTS) $(bsam_test_DEPENDENCIES)
	@rm -f bsam.test
	$(LINK) $(bsam_test_LDFLAGS) $(bsam_test_OBJECTS) $(bsam_test_LDADD) $(LIBS)

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

subdir = src/hub

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
