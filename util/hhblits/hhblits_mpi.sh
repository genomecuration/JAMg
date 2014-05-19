#!/bin/bash

# THIS is an example script to run PFAM on the grid with MPI. PBS code follows but ought to be straightforward to convert to LSF
# You will need ffindex 0.9.9.3+ installed in your path.

#number of nodes to use (1 CPU per node)
NUMBERSPROCESSES=40

# protein file to process (no *s as stop codons at the end please, e.g sed ~s/\*$// PROTEIN_FILE
PROTEIN_FILE=

if [ ! -e $PROTEIN_FILE.db ]; then
 sed -i '~s/\*$//' $PROTEIN_FILE
 ffindex_from_fasta $PROTEIN_FILE.db $PROTEIN_FILE.db.idx $PROTEIN_FILE
 mv $PROTEIN_FILE.db.idx $PROTEIN_FILE.db.idx.orig ; cp $PROTEIN_FILE.db.idx.orig $PROTEIN_FILE.db.idx.orig.notdone; ln -s $PROTEIN_FILE.db.idx.orig.notdone $PROTEIN_FILE.db.idx
fi

qsub -l select=$NUMBERSPROCESSES:ncpus=1:mpiprocs=1:mem=7gb:NodeType=any -l walltime=48:00:00 -V -A sf-CSIRO -r n -N hb -- $PWD/hhblits_mpi.pbs $PROTEIN_FILE.db $NUMBERSPROCESSES


