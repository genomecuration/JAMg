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

bin_PROGRAMS = exonerate ipcress exonerate-server

INCLUDES = -I$(top_srcdir)/src/sequence              -I$(top_srcdir)/src/database              -I$(top_srcdir)/src/comparison            -I$(top_srcdir)/src/struct                -I$(top_srcdir)/src/c4                    -I$(top_srcdir)/src/bsdp                  -I$(top_srcdir)/src/sdp                   -I$(top_srcdir)/src/model                 -I$(top_srcdir)/src/hub                   -I$(top_srcdir)/src/general               -DHOSTTYPE="\"x86_64-unknown-linux-gnu\""                   -DSOURCE_ROOT_DIR="\"/home/pap056/workspace/transcripts4community/jamg/3rd_party/exonerate-2.2.0\""


# Missing viterb.o and scheduler.o so codegen versions may be substituted
C4_OBJECTS = $(top_srcdir)/src/c4/c4.o                             $(top_srcdir)/src/c4/codegen.o                        $(top_srcdir)/src/c4/cgutil.o                         $(top_srcdir)/src/c4/opair.o                          $(top_srcdir)/src/c4/alignment.o                      $(top_srcdir)/src/c4/optimal.o                        $(top_srcdir)/src/c4/layout.o                         $(top_srcdir)/src/c4/region.o                         $(top_srcdir)/src/c4/subopt.o                         $(top_srcdir)/src/bsdp/bsdp.o                         $(top_srcdir)/src/bsdp/sar.o                          $(top_srcdir)/src/bsdp/heuristic.o                    $(top_srcdir)/src/bsdp/hpair.o                        $(top_srcdir)/src/sdp/boundary.o                      $(top_srcdir)/src/sdp/sdp.o                           $(top_srcdir)/src/sdp/lookahead.o                     $(top_srcdir)/src/sdp/straceback.o                    $(top_srcdir)/src/struct/matrix.o                     $(top_srcdir)/src/struct/slist.o                      $(top_srcdir)/src/struct/pqueue.o                     $(top_srcdir)/src/struct/recyclebin.o                 $(top_srcdir)/src/struct/rangetree.o                  $(top_srcdir)/src/sequence/sequence.o                 $(top_srcdir)/src/sequence/alphabet.o                 $(top_srcdir)/src/sequence/translate.o                $(top_srcdir)/src/sequence/splice.o                   $(top_srcdir)/src/model/affine.o                      $(top_srcdir)/src/model/est2genome.o                  $(top_srcdir)/src/model/ner.o                         $(top_srcdir)/src/model/protein2dna.o                 $(top_srcdir)/src/model/protein2genome.o              $(top_srcdir)/src/model/coding2coding.o               $(top_srcdir)/src/model/coding2genome.o               $(top_srcdir)/src/model/cdna2genome.o                 $(top_srcdir)/src/model/genome2genome.o               $(top_srcdir)/src/model/ungapped.o                    $(top_srcdir)/src/model/intron.o                      $(top_srcdir)/src/model/frameshift.o                  $(top_srcdir)/src/model/phase.o                       $(top_srcdir)/src/model/modeltype.o                   $(top_srcdir)/src/general/argument.o                  $(top_srcdir)/src/general/lineparse.o                 $(top_srcdir)/src/general/socket.o                    $(top_srcdir)/src/comparison/match.o                  $(top_srcdir)/src/sequence/submat.o                   $(top_srcdir)/src/sequence/codonsubmat.o


EXTRA_exonerate_SOURCES = c4_model_archive.a c4_model_archive.h

BUILT_SOURCES = c4_model_archive.a
BOOTSTRAPPER = $(top_srcdir)/src/model/bootstrapper

exonerate_COMPONENTS = viterbi.o scheduler.o c4_model_archive.a                                            $(top_srcdir)/src/comparison/wordhood.o                          $(top_srcdir)/src/comparison/hspset.o                            $(top_srcdir)/src/comparison/seeder.o                            $(top_srcdir)/src/comparison/comparison.o                        $(top_srcdir)/src/database/fastapipe.o                           $(top_srcdir)/src/database/fastadb.o                             $(top_srcdir)/src/struct/dejavu.o                                $(top_srcdir)/src/struct/fsm.o                                   $(top_srcdir)/src/struct/vfsm.o                                  $(top_srcdir)/src/struct/sparsecache.o                           $(top_srcdir)/src/hub/gam.o                                      $(top_srcdir)/src/hub/bsam.o                                     $(top_srcdir)/src/hub/analysis.o                                 $(top_srcdir)/src/general/compoundfile.o                         $(C4_OBJECTS)


exonerate_SOURCES = exonerate.c
exonerate_DEPENDENCIES = c4_model_archive.h Makefile                          $(exonerate_COMPONENTS)

exonerate_LDADD = $(exonerate_COMPONENTS) -lm

ipcress_SOURCES = ipcress.c
ipcress_LDADD = $(top_srcdir)/src/comparison/pcr.o                       $(top_srcdir)/src/comparison/wordhood.o                  $(top_srcdir)/src/database/fastadb.o                     $(top_srcdir)/src/sequence/alphabet.o                    $(top_srcdir)/src/sequence/sequence.o                    $(top_srcdir)/src/sequence/splice.o                      $(top_srcdir)/src/sequence/translate.o                   $(top_srcdir)/src/sequence/submat.o                      $(top_srcdir)/src/sequence/codonsubmat.o                 $(top_srcdir)/src/struct/fsm.o                           $(top_srcdir)/src/struct/slist.o                         $(top_srcdir)/src/struct/matrix.o                        $(top_srcdir)/src/struct/recyclebin.o                    $(top_srcdir)/src/struct/sparsecache.o                   $(top_srcdir)/src/general/lineparse.o                    $(top_srcdir)/src/general/argument.o                     $(top_srcdir)/src/general/compoundfile.o                 -lm

# Files to clear away

exonerate_server_SOURCES = exonerate-server.c
exonerate_server_LDADD = $(top_srcdir)/src/general/socket.o                                  $(top_srcdir)/src/general/argument.o                                $(top_srcdir)/src/comparison/hspset.o                               $(top_srcdir)/src/comparison/match.o                                $(top_srcdir)/src/comparison/wordhood.o                             $(top_srcdir)/src/database/dataset.o                                $(top_srcdir)/src/database/index.o                                  $(top_srcdir)/src/database/fastadb.o                                $(top_srcdir)/src/struct/bitarray.o                                 $(top_srcdir)/src/struct/sparsecache.o                              $(top_srcdir)/src/struct/matrix.o                                   $(top_srcdir)/src/struct/splaytree.o                                $(top_srcdir)/src/struct/noitree.o                                  $(top_srcdir)/src/struct/vfsm.o                                     $(top_srcdir)/src/struct/pqueue.o                                   $(top_srcdir)/src/struct/rangetree.o                                $(top_srcdir)/src/struct/recyclebin.o                               $(top_srcdir)/src/sequence/sequence.o                               $(top_srcdir)/src/sequence/alphabet.o                               $(top_srcdir)/src/sequence/translate.o                              $(top_srcdir)/src/sequence/splice.o                                 $(top_srcdir)/src/sequence/submat.o                                 $(top_srcdir)/src/sequence/codonsubmat.o                            $(top_srcdir)/src/general/compoundfile.o                            $(top_srcdir)/src/general/lineparse.o                               -lm


CLEANFILES = $(EXTRA_exonerate_SOURCES) c4_model_archive.h              @codegen_extra_ld_add@

MAINTAINERCLEANFILES = Makefile.in
mkinstalldirs = $(SHELL) $(top_srcdir)/mkinstalldirs
CONFIG_CLEAN_FILES = 
PROGRAMS =  $(bin_PROGRAMS)


DEFS = -DPACKAGE_NAME=\"\" -DPACKAGE_TARNAME=\"\" -DPACKAGE_VERSION=\"\" -DPACKAGE_STRING=\"\" -DPACKAGE_BUGREPORT=\"\" -DPACKAGE=\"exonerate\" -DVERSION=\"2.2.0\" -DSTDC_HEADERS=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_STDLIB_H=1 -DHAVE_STRING_H=1 -DHAVE_MEMORY_H=1 -DHAVE_STRINGS_H=1 -DHAVE_INTTYPES_H=1 -DHAVE_STDINT_H=1 -DHAVE_UNISTD_H=1 -DHAVE_LIMITS_H=1 -DHAVE_ERRNO_H=1 -DHAVE_MATH_H=1 -DHAVE_STRINGS_H=1 -DHAVE_CTYPE_H=1 -DHAVE_SYS_STAT_H=1 -DHAVE_SYS_TYPES_H=1 -DHAVE_TIME_H=1 -DHAVE_UNISTD_H=1 -I. -I$(srcdir) 
CPPFLAGS =  -fPIC 
LDFLAGS =  -lpthread
LIBS =  -lglib-2.0  
exonerate_OBJECTS =  exonerate.o
exonerate_LDFLAGS = 
ipcress_OBJECTS =  ipcress.o
ipcress_DEPENDENCIES =  $(top_srcdir)/src/comparison/pcr.o \
$(top_srcdir)/src/comparison/wordhood.o \
$(top_srcdir)/src/database/fastadb.o \
$(top_srcdir)/src/sequence/alphabet.o \
$(top_srcdir)/src/sequence/sequence.o \
$(top_srcdir)/src/sequence/splice.o \
$(top_srcdir)/src/sequence/translate.o \
$(top_srcdir)/src/sequence/submat.o \
$(top_srcdir)/src/sequence/codonsubmat.o $(top_srcdir)/src/struct/fsm.o \
$(top_srcdir)/src/struct/slist.o $(top_srcdir)/src/struct/matrix.o \
$(top_srcdir)/src/struct/recyclebin.o \
$(top_srcdir)/src/struct/sparsecache.o \
$(top_srcdir)/src/general/lineparse.o \
$(top_srcdir)/src/general/argument.o \
$(top_srcdir)/src/general/compoundfile.o
ipcress_LDFLAGS = 
exonerate_server_OBJECTS =  exonerate-server.o
exonerate_server_DEPENDENCIES =  $(top_srcdir)/src/general/socket.o \
$(top_srcdir)/src/general/argument.o \
$(top_srcdir)/src/comparison/hspset.o \
$(top_srcdir)/src/comparison/match.o \
$(top_srcdir)/src/comparison/wordhood.o \
$(top_srcdir)/src/database/dataset.o $(top_srcdir)/src/database/index.o \
$(top_srcdir)/src/database/fastadb.o \
$(top_srcdir)/src/struct/bitarray.o \
$(top_srcdir)/src/struct/sparsecache.o \
$(top_srcdir)/src/struct/matrix.o $(top_srcdir)/src/struct/splaytree.o \
$(top_srcdir)/src/struct/noitree.o $(top_srcdir)/src/struct/vfsm.o \
$(top_srcdir)/src/struct/pqueue.o $(top_srcdir)/src/struct/rangetree.o \
$(top_srcdir)/src/struct/recyclebin.o \
$(top_srcdir)/src/sequence/sequence.o \
$(top_srcdir)/src/sequence/alphabet.o \
$(top_srcdir)/src/sequence/translate.o \
$(top_srcdir)/src/sequence/splice.o $(top_srcdir)/src/sequence/submat.o \
$(top_srcdir)/src/sequence/codonsubmat.o \
$(top_srcdir)/src/general/compoundfile.o \
$(top_srcdir)/src/general/lineparse.o
exonerate_server_LDFLAGS = 
CFLAGS = -g -O2 -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include   -DG_DISABLE_ASSERT -O3 -finline-functions -D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -DUSE_PTHREADS
COMPILE = $(CC) $(DEFS) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)
CCLD = $(CC)
LINK = $(CCLD) $(AM_CFLAGS) $(CFLAGS) $(LDFLAGS) -o $@
DIST_COMMON =  Makefile.am Makefile.in


DISTFILES = $(DIST_COMMON) $(SOURCES) $(HEADERS) $(TEXINFOS) $(EXTRA_DIST)

TAR = tar
GZIP_ENV = --best
SOURCES = $(exonerate_SOURCES) $(EXTRA_exonerate_SOURCES) $(ipcress_SOURCES) $(exonerate_server_SOURCES)
OBJECTS = $(exonerate_OBJECTS) $(ipcress_OBJECTS) $(exonerate_server_OBJECTS)

all: all-redirect
.SUFFIXES:
.SUFFIXES: .S .c .o .s
$(srcdir)/Makefile.in: Makefile.am $(top_srcdir)/configure.in $(ACLOCAL_M4) 
	cd $(top_srcdir) && $(AUTOMAKE) --gnu --include-deps src/program/Makefile

Makefile: $(srcdir)/Makefile.in  $(top_builddir)/config.status
	cd $(top_builddir) \
	  && CONFIG_FILES=$(subdir)/$@ CONFIG_HEADERS= $(SHELL) ./config.status


mostlyclean-binPROGRAMS:

clean-binPROGRAMS:
	-test -z "$(bin_PROGRAMS)" || rm -f $(bin_PROGRAMS)

distclean-binPROGRAMS:

maintainer-clean-binPROGRAMS:

install-binPROGRAMS: $(bin_PROGRAMS)
	@$(NORMAL_INSTALL)
	$(mkinstalldirs) $(DESTDIR)$(bindir)
	@list='$(bin_PROGRAMS)'; for p in $$list; do \
	  if test -f $$p; then \
	    echo "  $(INSTALL_PROGRAM) $$p $(DESTDIR)$(bindir)/`echo $$p|sed 's/$(EXEEXT)$$//'|sed '$(transform)'|sed 's/$$/$(EXEEXT)/'`"; \
	     $(INSTALL_PROGRAM) $$p $(DESTDIR)$(bindir)/`echo $$p|sed 's/$(EXEEXT)$$//'|sed '$(transform)'|sed 's/$$/$(EXEEXT)/'`; \
	  else :; fi; \
	done

uninstall-binPROGRAMS:
	@$(NORMAL_UNINSTALL)
	list='$(bin_PROGRAMS)'; for p in $$list; do \
	  rm -f $(DESTDIR)$(bindir)/`echo $$p|sed 's/$(EXEEXT)$$//'|sed '$(transform)'|sed 's/$$/$(EXEEXT)/'`; \
	done

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

exonerate: $(exonerate_OBJECTS) $(exonerate_DEPENDENCIES)
	@rm -f exonerate
	$(LINK) $(exonerate_LDFLAGS) $(exonerate_OBJECTS) $(exonerate_LDADD) $(LIBS)

ipcress: $(ipcress_OBJECTS) $(ipcress_DEPENDENCIES)
	@rm -f ipcress
	$(LINK) $(ipcress_LDFLAGS) $(ipcress_OBJECTS) $(ipcress_LDADD) $(LIBS)

exonerate-server: $(exonerate_server_OBJECTS) $(exonerate_server_DEPENDENCIES)
	@rm -f exonerate-server
	$(LINK) $(exonerate_server_LDFLAGS) $(exonerate_server_OBJECTS) $(exonerate_server_LDADD) $(LIBS)

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

subdir = src/program

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
	$(MAKE) $(AM_MAKEFLAGS) top_distdir="$(top_distdir)" distdir="$(distdir)" dist-hook
info-am:
info: info-am
dvi-am:
dvi: dvi-am
check-am: all-am
check: check-am
installcheck-am:
installcheck: installcheck-am
install-exec-am: install-binPROGRAMS
install-exec: install-exec-am

install-data-am:
install-data: install-data-am

install-am: all-am
	@$(MAKE) $(AM_MAKEFLAGS) install-exec-am install-data-am
install: install-am
uninstall-am: uninstall-binPROGRAMS
uninstall: uninstall-am
all-am: Makefile $(PROGRAMS)
all-redirect: all-am
install-strip:
	$(MAKE) $(AM_MAKEFLAGS) AM_INSTALL_PROGRAM_FLAGS=-s install
installdirs:
	$(mkinstalldirs)  $(DESTDIR)$(bindir)


mostlyclean-generic:

clean-generic:
	-test -z "$(CLEANFILES)" || rm -f $(CLEANFILES)

distclean-generic:
	-rm -f Makefile $(CONFIG_CLEAN_FILES)
	-rm -f config.cache config.log stamp-h stamp-h[0-9]*

maintainer-clean-generic:
	-test -z "$(BUILT_SOURCES)$(MAINTAINERCLEANFILES)" || rm -f $(BUILT_SOURCES) $(MAINTAINERCLEANFILES)
mostlyclean-am:  mostlyclean-binPROGRAMS mostlyclean-compile \
		mostlyclean-tags mostlyclean-generic

mostlyclean: mostlyclean-am

clean-am:  clean-binPROGRAMS clean-compile clean-tags clean-generic \
		mostlyclean-am

clean: clean-am

distclean-am:  distclean-binPROGRAMS distclean-compile distclean-tags \
		distclean-generic clean-am

distclean: distclean-am

maintainer-clean-am:  maintainer-clean-binPROGRAMS \
		maintainer-clean-compile maintainer-clean-tags \
		maintainer-clean-generic distclean-am
	@echo "This command is intended for maintainers to use;"
	@echo "it deletes files that may require special tools to rebuild."

maintainer-clean: maintainer-clean-am

.PHONY: mostlyclean-binPROGRAMS distclean-binPROGRAMS clean-binPROGRAMS \
maintainer-clean-binPROGRAMS uninstall-binPROGRAMS install-binPROGRAMS \
mostlyclean-compile distclean-compile clean-compile \
maintainer-clean-compile tags mostlyclean-tags distclean-tags \
clean-tags maintainer-clean-tags distdir info-am info dvi-am dvi check \
check-am installcheck-am installcheck install-exec-am install-exec \
install-data-am install-data install-am install uninstall-am uninstall \
all-redirect all-am all installdirs mostlyclean-generic \
distclean-generic clean-generic maintainer-clean-generic clean \
mostlyclean distclean maintainer-clean

c4_model_archive.h: $(BOOTSTRAPPER)
	$(BOOTSTRAPPER) --compiled no
c4_model_archive.a: $(BOOTSTRAPPER)
	$(BOOTSTRAPPER) --compiled no

# Special viterbi.o and scheduler.o for compiled model linking
viterbi.o: $(top_srcdir)/src/c4/viterbi.c
	$(CC) $(CFLAGS) $(INCLUDES) -I.   \
        -DUSE_COMPILED_MODELS             \
        -c $(top_srcdir)/src/c4/viterbi.c
scheduler.o: $(top_srcdir)/src/sdp/scheduler.c
	$(CC) $(CFLAGS) $(INCLUDES) -I.   \
        -DUSE_COMPILED_MODELS             \
        -c $(top_srcdir)/src/sdp/scheduler.c

# Remove built sources from distribution
# (will be able to use nodist_exonerate_SOURCES in automake 1.5)
# nodist_exonerate_SOURCES = c4_model_archive.a
dist-hook:
	for file in $(EXTRA_exonerate_SOURCES) ; \
            do rm -f $(distdir)/$$file ; done ;

# Tell versions [3.59,3.63) of GNU make to not export all variables.
# Otherwise a system limit (for SysV at least) may be exceeded.
.NOEXPORT:
