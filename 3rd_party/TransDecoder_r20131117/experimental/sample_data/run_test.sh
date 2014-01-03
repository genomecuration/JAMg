#!/bin/bash
tar xzf Dmel_refseq.tar.gz
tar xzf mito_aa_inv_80.fsa.trim.tar.gz
tar xzf rDNA_nt_inv.fsa_nr.tar.gz
time ../extract_cds_on_frameshifts.pl -in test_seq.fsa -mt mito_aa_inv_80.fsa.trim -rrna rDNA_nt_inv.fsa_nr -nuc Dmel_refseq -gencon_mt 5 -gencon_nuc 1 -debug 2>&1| tee test.log
time ../../transcripts_to_best_scoring_ORFs.pl -t test_seq.fsa.untranslated.fsa --train test_seq.fsa_protpred_nuclear_nuc.fsa
