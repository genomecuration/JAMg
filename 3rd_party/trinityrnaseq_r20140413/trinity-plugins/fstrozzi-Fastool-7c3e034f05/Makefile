CC     = gcc
CFLAGS += -O2 -std=c99
LDFLAGS +=

all:kseq.h fastool.c
	$(CC) $(CFLAGS) $(LDFLAGS) fastool.c -o fastool -lz

clean:
	rm -f *.o fastool
