######################
# Makefile for SNAP  #
######################

LIB = -lm
INC = -IZoe

OBJECTS = \
	Zoe/zoeAlignment.o\
	Zoe/zoeCDS.o\
	Zoe/zoeCounter.o\
	Zoe/zoeDNA.o\
	Zoe/zoeDistribution.o\
	Zoe/zoeDuration.o\
	Zoe/zoeFastaFile.o\
	Zoe/zoeFeature.o\
	Zoe/zoeFeatureFactory.o\
	Zoe/zoeFeatureTable.o\
	Zoe/zoeHMM.o\
	Zoe/zoeIsochore.o\
	Zoe/zoeMath.o\
	Zoe/zoeModel.o\
	Zoe/zoePhasePref.o\
	Zoe/zoeProtein.o\
	Zoe/zoeScanner.o\
	Zoe/zoeState.o\
	Zoe/zoeTools.o\
	Zoe/zoeTransition.o\
	Zoe/zoeTrellis.o\

APP = snap
SRC = snap.c
OBJ = snap.o

APP2 = fathom
SRC2 = fathom.c
OBJ2 = fathom.o

APP3 = forge
SRC3 = forge.c
OBJ3 = forge.o

APP4 = hmm-info
SRC4 = hmm-info.c
OBJ4 = hmm-info.o

APP5 = exonpairs
SRC5 = exonpairs.c
OBJ5 = exonpairs.o

DATE = $(shell date +\%Y-\%m-\%d)

###########
# Targets #
###########

default:
	make gcc

$(APP): $(OBJ) $(OBJECTS)
	$(CC) -o $(APP) $(CFLAGS) $(OBJ) $(OBJECTS) $(LIB)

$(APP2): $(OBJ2) $(OBJECTS)
	$(CC) -o $(APP2) $(CFLAGS) $(OBJ2) $(OBJECTS) $(LIB)

$(APP3): $(OBJ3) $(OBJECTS)
	$(CC) -o $(APP3) $(CFLAGS) $(OBJ3) $(OBJECTS) $(LIB)

$(APP4): $(OBJ4) $(OBJECTS)
	$(CC) -o $(APP4) $(CFLAGS) $(OBJ4) $(OBJECTS) $(LIB)

$(APP5): $(OBJ5) $(OBJECTS)
	$(CC) -o $(APP5) $(CFLAGS) $(OBJ5) $(OBJECTS) $(LIB)

clean:
	rm -f *.o $(APP) $(APP2) $(APP3) $(APP4) $(APP5)
	cd Zoe; make clean

depend: $(OBJECTS:.o=.c)
	gcc $(INC) -MM $^ > $@

tar:
	rm -rf /tmp/$(APP)
	mkdir /tmp/$(APP)
	cp -r * /tmp/$(APP)
	cd /tmp/$(APP); make clean;  rm -rf CVS */CVS
	cd /tmp; tar -zcvf $(APP)-$(DATE).tar.gz $(APP)
	rm -rf /tmp/$(APP)


#################
# Architectures #
#################

gcc:
	cd Zoe; make;
	make $(APP)  CC="gcc" CFLAGS="-O2 -Wall -Werror"
	make $(APP2) CC="gcc" CFLAGS="-O2 -Wall -Werror"
	make $(APP3) CC="gcc" CFLAGS="-O2 -Wall -Werror"
	make $(APP4) CC="gcc" CFLAGS="-O2 -Wall -Werror"
	make $(APP5) CC="gcc" CFLAGS="-O2 -Wall -Werror"


###################
# Inference Rules #
###################

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

################
# Dependancies #
################

include depend

