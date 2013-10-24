#!/bin/bash

# the data where prepared using transdecoder
echo "Predicting proteins"
../3rd_party/transdecoder/transcripts_to_best_scoring_ORFs.pl -t test.fsa

# prepare genome. Softmasked using RepeatMasker's GFF and bedtools
bunzip2 -k dmel-X-r5.53.repeat.hard.bz2 dmel-X-r5.53.repeat.soft.bz2

# predict genes. Don't produce genome annotation (not installed by default)
echo "Predicting genes"
../bin/prepare_golden_genes_for_predictors.pl -genome dmel-X-r5.53.repeat.hard -softmasked dmel-X-r5.53.repeat.soft -same_species \
 -transdecoder_gff test.fsa.transdecoder.gff3 -transdecoder_peptides test.fsa.transdecoder.pep -transdecoder_assembly test.fsa \
 -threads 2 -augustus ../3rd_party/augustus -norefine # options for speed/this demo

find . -type f -empty -delete

echo "Completed."
