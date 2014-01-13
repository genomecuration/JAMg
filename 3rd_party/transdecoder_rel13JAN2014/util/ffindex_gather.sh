#!/bin/bash

# if you don't specify an argument file then it will look for *.all.db outputs and strip the 0 delimiting character
# to make a proper text file output

if [ $1 ]; then
        if [ -e $1.db ]; then
	   ffindex_build -as $1.all.db $1.all.idx -d $1.db -i $1.db.idx
	   rm -f $1.db $1.db.idx
	fi
	for i in {0..100}; \
         do 
          if [ -e $1.db.$i ]; then 
	   ffindex_build -as $1.all.db $1.all.idx -d $1.db.$i -i $1.db.idx.$i
           rm -f $1.db.$i $1.db.idx.$i
          fi
	 done
        ffindex_resume.pl *orig *all.idx
else
    cat  *.all.db |tr -d '\000' >  results.all.db.hhr
fi

