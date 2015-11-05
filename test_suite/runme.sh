#!/bin/bash

# the data where prepared using transdecoder
echo "Normally we use PASA but for this demo we will just use TransDecoder for quickly predicting proteins"
../3rd_party/transdecoder/TransDecoder.LongOrfs -t test.fsa
../3rd_party/transdecoder/TransDecoder.Predict -t test.fsa

# prepare genome. Softmasked using RepeatMasker's GFF and bedtools
bunzip2 -k dmel-X-r5.53.repeat.hard.bz2 dmel-X-r5.53.repeat.soft.bz2

# predict genes. Don't produce genome annotation (not installed by default)
echo "Predicting genes"
../bin/prepare_golden_genes_for_predictors.pl -genome dmel-X-r5.53.repeat.hard -softmasked dmel-X-r5.53.repeat.soft -same_species \
 -peptide test.fsa.transdecoder.pep \
 -threads 2 -augustus ../3rd_party/augustus -norefine

find . -type f -empty -delete

echo "Completed."
