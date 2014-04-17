
TESTS = argument.test lineparse.test compoundfile.test socket.test
noinst_PROGRAMS = $(TESTS)

noinst_HEADERS = argument.h lineparse.h compoundfile.h socket.h

INCLUDES = -DSOURCE_ROOT_DIR="\"@source_root_dir@\""

argument_test_SOURCES = argument.test.c argument.c
lineparse_test_SOURCES = lineparse.test.c lineparse.c
compoundfile_test_SOURCES = compoundfile.test.c compoundfile.c
socket_test_SOURCES = socket.test.c socket.c

# Files to clear away

MAINTAINERCLEANFILES = Makefile.in

