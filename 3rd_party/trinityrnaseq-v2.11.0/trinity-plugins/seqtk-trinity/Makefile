CC=gcc
CFLAGS=-g -Wall -O2 -Wno-unused-function 

all:seqtk-trinity

seqtk-trinity:seqtk.c khash.h kseq.h
		$(CC) $(CFLAGS) seqtk.c -o $@ -lz -lm

clean:
		rm -fr gmon.out *.o ext/*.o a.out seqtk trimadap *~ *.a *.dSYM session*


test:
	./seqtk-trinity seq -A -R 1 testing/oldformat_1.fq > test.fa
	./seqtk-trinity seq -A -R 2 testing/oldformat_2.fq > test.fa	

	./seqtk-trinity seq -A -R 1 testing/newformat_1.fq > test.fa
	./seqtk-trinity seq -A -R 2 testing/newformat_2.fq > test.fa	

	./seqtk-trinity seq -A -R 1 testing/stripped.fq > test.fa
	./seqtk-trinity seq -A -R 2 testing/stripped.fq > test.fa

	./seqtk-trinity seq -A -R 1 testing/forward_reverse_notation_1.fq > test.fa
	./seqtk-trinity seq -A -R 2 testing/forward_reverse_notation_2.fq > test.fa


