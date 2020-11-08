#!/bin/bash
/usr/lib/trinityrnaseq/util/bowtie_PE_separate_then_join.pl --seqType fq --left ../reads.left.fq --right ../reads.right.fq --target refSeqs.fa --aligner bowtie -- -p 4 --all --best --strata -m 300
