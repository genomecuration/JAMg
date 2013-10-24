#Makefile

CC= gcc
CFLAGS= -O
IFLAGS= -I.
LFLAGS= -lm

BIN= .

.KEEP_STATE:

PROGRAMS= 	pictogram

.c.o:
	$(CC) $(CFLAGS) $(IFLAGS) -c $<

all: $(PROGRAMS)

pictogram :	pictogram.c
		$(CC) $(CFLAGS) $(IFLAGS) pictogram.c $(LFLAGS) -o $(BIN)/pictogram



