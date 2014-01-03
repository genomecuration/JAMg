#!/bin/bash

# THIS is an example script to run PFAM on the grid with MPI. PBS code follows but ought to be straightforward to convert to LSF
# You will need ffindex 0.9.9.3+ installed in your path.

# Copy this and pfam_mpi.pbs in your WORKING directory that and type ./pfam_mpi.sh my_protein.fasta (or see PROTEIN_FILE variable)
# You do NOT need to use partioned FASTA (ffindex will take care of all that). So just provide the longest_orfs.pep file
PROTEIN_FILE=longest_orfs.pep

#number of nodes to use (1 CPU per node)
NUMBERSPROCESSES=60

if [ ! -e $PROTEIN_FILE.db ]; then
 ffindex_from_fasta $PROTEIN_FILE.db $PROTEIN_FILE.db.idx $PROTEIN_FILE
 mv $PROTEIN_FILE.db.idx $PROTEIN_FILE.db.idx.orig ; cp $PROTEIN_FILE.db.idx.orig $PROTEIN_FILE.db.idx.orig.notdone; ln -s $PROTEIN_FILE.db.idx.orig.notdone $PROTEIN_FILE.db.idx
fi

#Check the QSUB string
qsub -l select=$NUMBERSPROCESSES:ncpus=1:mpiprocs=1:mem=7gb:NodeType=any -l walltime=20:00:00 -V -A sf-CSIRO -r n -N pfam -- $PWD/pfam_mpi.pbs $FILE $NUMBERSPROCESSES


# Once the run is complete, you can use ffindex_gather.sh to see if there are unfinished parts and remove the \0 delimiting character
# If there are, ffindex_gather will simply make a new index and you can just relaunch this script without changing anything
