#!/bin/bash

LOWMEM_CPUS=8
HIGHMEM_CPUS=4

rm -f commands.large commands.small

trim_path=$(which trim_fasta_all.pl)
if [ ! $trim_path ] ; then
   echo "Cannot find trim_fasta_all.pl in the path. Have you installed JAMg?"
   exit
fi

ffindex_path=$(which ffindex_from_fasta)
if [ ! $ffindex_path ] ; then
   echo "Cannot find ffindex_from_fasta in the path. Have you installed JAMg?"
   exit
fi

# THIS WILL FAIL IF IDs are TOO long >32
PROTEIN_FILE=$1

if [ ! $PROTEIN_FILE ]; then
 echo Please provide a protein file
 exit
fi

DBDIR=/dev/shm/$USER/
if [ ! $DBDIR ]; then
 DBDIR=$PWD
fi

if [ ! $DBDIR/swissprot20_2016_02 ]; then
 echo Cannot find $DBDIR/swissprot20_2016_02*
 exit
fi

sed -i '~s/\*$//' $PROTEIN_FILE

$trim_path -le 5000 -over -in $PROTEIN_FILE -out $PROTEIN_FILE.le5000

#SMALL
rm -f $PROTEIN_FILE.le5000.ff{data,index}
ffindex_from_fasta $PROTEIN_FILE.le5000.ff{data,index} $PROTEIN_FILE.le5000 > /dev/null
shuf $PROTEIN_FILE.le5000.ffindex > $PROTEIN_FILE.le5000.ffindex.
mv -f $PROTEIN_FILE.le5000.ffindex. $PROTEIN_FILE.le5000.ffindex

split -d -a 3 -l 100 $PROTEIN_FILE.le5000.ffindex $PROTEIN_FILE.le5000.ffindex.
let NUMBERSPROCESSES=`ls -l $PROTEIN_FILE.le5000.ffindex.???|wc -l`-1

for (( i=0; i<=$NUMBERSPROCESSES; i++ )); do
	INDEX_PADDED=`printf "%03d" $i`
	IDX=$PROTEIN_FILE.le5000.$INDEX_PADDED.ffindex
	mv -f $PROTEIN_FILE.le5000.ffindex.$INDEX_PADDED $IDX
	DTX=$PROTEIN_FILE.le5000.$INDEX_PADDED.ffdata
	ln -s $PROTEIN_FILE.le5000.ffdata $DTX
	echo hhblits_omp -maxmem 3 -d $DBDIR/swissprot20_2016_02 -mact 0.3 -cpu $LOWMEM_CPUS -i $PROTEIN_FILE.le5000.$INDEX_PADDED -o $PROTEIN_FILE.le5000.$INDEX_PADDED.out \
         -id 100 -noprefilt -maxres 15000 -e 0.0001 -p 60 -E 1E-03 -z 0 -b 0 -v 0 -n 2 >> commands.small
done

	echo See commands.small

if [ -f $PROTEIN_FILE.le5000.discard ]; then
	#LARGE
	rm -f $PROTEIN_FILE.le5000.discard.ff{data,index}
	ffindex_from_fasta $PROTEIN_FILE.le5000.discard.ff{data,index} $PROTEIN_FILE.le5000.discard > /dev/null
	shuf $PROTEIN_FILE.le5000.discard.ffindex > $PROTEIN_FILE.le5000.discard.ffindex.
	mv -f $PROTEIN_FILE.le5000.discard.ffindex. $PROTEIN_FILE.le5000.discard.ffindex

	split -d -a 3 -l 10 $PROTEIN_FILE.le5000.discard.ffindex $PROTEIN_FILE.le5000.discard.ffindex.
	let NUMBERSPROCESSES=`ls -l $PROTEIN_FILE.le5000.discard.ffindex.???|wc -l`-1

	for (( i=0; i<=$NUMBERSPROCESSES; i++ )); do
		INDEX_PADDED=`printf "%03d" $i`
		IDX=$PROTEIN_FILE.le5000.discard.$INDEX_PADDED.ffindex
		mv -f $PROTEIN_FILE.le5000.discard.ffindex.$INDEX_PADDED $IDX
		DTX=$PROTEIN_FILE.le5000.discard.$INDEX_PADDED.ffdata
		ln -s $PROTEIN_FILE.le5000.discard.ffdata $DTX
		echo hhblits_omp -maxmem 10 -d $DBDIR/swissprot20_2016_02 -mact 0.3 -cpu $HIGHMEM_CPUS -i $PROTEIN_FILE.le5000.discard.$INDEX_PADDED -o $PROTEIN_FILE.le5000.discard.$INDEX_PADDED.out \
        	 -id 100 -noprefilt -maxres 431000 -e 0.0001 -p 60 -E 1E-03 -z 0 -b 0 -v 0 -n 2 >> commands.large
	done

	echo See commands.large
fi

