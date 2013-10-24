# Makefile for SSgff

INCLUDE= ./include
CDIR= ./src
OBJ = ./objects
HEADERS = $(INCLUDE)/SSgff.h 
PROGRAM= SSgff
PRODUCT= bin/$(PROGRAM)
CC= gcc
OPTS=-I$(INCLUDE) -Wall 

#######

OBJECTS = $(OBJ)/input.o $(OBJ)/output.o $(OBJ)/scan.o $(OBJ)/Dictionary.o $(OBJ)/readargv.o 

#######

$(PRODUCT): $(OBJ)/$(PROGRAM).o $(OBJECTS) $(HEADERS) 
	$(CC) $(OPTS) -o $(PRODUCT) $(OBJ)/$(PROGRAM).o $(OBJECTS) -lm

$(OBJ)/$(PROGRAM).o : $(CDIR)/$(PROGRAM).c  $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/$(PROGRAM).c -o $(OBJ)/$(PROGRAM).o 

$(OBJ)/input.o :  $(CDIR)/input.c $(HEADERS) 
	$(CC) -c $(OPTS) $(CDIR)/input.c -o $(OBJ)/input.o

$(OBJ)/output.o :  $(CDIR)/output.c $(HEADERS) 
	$(CC) -c $(OPTS) $(CDIR)/output.c -o $(OBJ)/output.o

$(OBJ)/scan.o : $(CDIR)/scan.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/scan.c -o $(OBJ)/scan.o

$(OBJ)/Dictionary.o : $(CDIR)/Dictionary.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/Dictionary.c -o $(OBJ)/Dictionary.o

$(OBJ)/readargv.o : $(CDIR)/readargv.c $(HEADERS)
	$(CC) -c $(OPTS) $(CDIR)/readargv.c -o $(OBJ)/readargv.o


clean:
	rm -f $(OBJ)/*.o $(PRODUCT) *~ $(INCLUDE)/*~ $(CDIR)/*~ core

