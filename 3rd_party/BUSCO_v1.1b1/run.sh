#!/bin/bash


source /home/aap599/software/BUSCO_v1.1b1/vars.export
BUSCO_v1.1b1.py -o BUSCO -in *fa -l $HOME/software/BUSCO_v1.1b1/arthropoda -m genome --cpu $NCPUS -sp fly

