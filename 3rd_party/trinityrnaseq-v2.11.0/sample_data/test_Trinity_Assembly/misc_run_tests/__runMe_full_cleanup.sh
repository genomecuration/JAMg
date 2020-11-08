#!/bin/bash -ve

#######################################################
##  Run Trinity to Generate Transcriptome Assemblies ##
#######################################################

## use jellyfish
../../Trinity --seqType fq --max_memory 2G --left reads.left.fq.gz --right reads.right.fq.gz --SS_lib_type RF --CPU 2 --full_cleanup

##### Done Running Trinity #####

if [ ! $* ]; then
    exit 0
fi

