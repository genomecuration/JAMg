# Gene Ontology Enrichment using Trinotate and GOseq 

The following describes how to use Trinotate and GOseq to explore functional enrichment among gene sets.

>This combined Trinotate/GOseq system hasn't been rigorously benchmarked, so use for exploratory purposes. Alternatives such as [Blast2GO](http://www.blast2go.com/b2ghome) should be considered as potentially more useful and/or accurate alternatives, yet the system described here can be quite useful and is readily accessible.

## Run Trinotate to Generate an Annotation Report

Use [Trinotate](https://github.com/Trinotate/Trinotate.github.io/wiki) to generate an annotation report 'trinotate.xls', which includes top-matching blast matches from SwissProt (and optionally UniRef) and any corresponding Gene Ontology (GO) assignments from the TrEMBL/SwissProt databases.

## Extract GO assignments per gene

Extract all GO assignments for each gene feature, including all parent terms within the GO DAG, using a script included in **Trinotate** (not Trinity) like so:

      ${TRINOTATE_HOME}/util/extract_GO_assignments_from_Trinotate_xls.pl \
                             --Trinotate_xls trinotate.xls \
                             -G --include_ancestral_terms \
                             > go_annotations.txt


## Run GOseq

The Bioconductor package [GOseq](http://www.bioconductor.org/packages/release/bioc/html/goseq.html) can then be used to perform functional enrichment tests, like so:

     ${TRINITY_HOME}/Analysis/DifferentialExpression/run_GOseq.pl \
                           --factor_labeling  factor_labeling.txt \
                           --GO_assignments go_annotations.txt \
                           --lengths gene.lengths.txt

Note, the above script is here in the **Trinity** software package (not Trinotate). Yes, it can be confusing going between the two.

The 'factor_labeling.txt' file should be of format:

     factor (tab) gene_id 

where factor is a string describing that subset of genes.

For example:

     diff_expressed_cond_X_Y (tab) my_gene_A 
     diff_expressed_cond_X_Y (tab) my_gene_B 
     ...
     diff_cond_W_Z (tab) my_gene_M
     diff_cond_W_Z (tab) my_gene_N
     ...

You can come up with whatever gene subset you want and call it whatever you want.  The enrichment tests will be performed separately for each factor defined.

The gene.lengths.txt file has the format

     ${gene_name} (tab) ${gene_length}

To create a gene lengths file, first create a file containing the transcript lengths:

    %  ${TRINITY_HOME}/util/misc/fasta_seq_length.pl  Trinity.fasta > Trinity.fasta.seq_lens

and then use that file when creating the gene lengths file:

    %  ${TRINITY_HOME}/util/misc/TPM_weighted_gene_length.py  \
             --gene_trans_map trinity_out_dir/Trinity.fasta.gene_trans_map \
             --trans_lengths Trinity.fasta.seq_lens \
             --TPM_matrix isoforms.TMM.EXPR.matrix > Trinity.gene_lengths.txt


## GO-Seq on Differentially Expressed Features

See [DE analysis](Trinity-Differential-Expression) for instructions on automatically running GO-Seq on significantly DE features.


## Output format

Two outputs will be generated for each set of genes tested for functional enrichment, one containing the enriched categories, and another containing the depleted categories (see .enriched and .depleted files).    The .enriched file will contain those categories that are found to have enriched representation among that set of genes.  Similarly, the .depleted file will contain those functional categories that are depleted (under-represented) among that very same set of target genes.

Each file is formatted like so: 

     category    over_represented_pvalue  under_represented_pvalue  numDEInCat  numInCat  term                 ontology       over_represented_FDR  go_term             gene_ids
     GO:0055114  2.40039957410791e-05     0.999993887954075         29          47        oxidation-reduction  process        BP                    0.0938076153561372  BP                 oxidation-reduction  process            comp1017_c0_seq1,  comp1172_c0_seq1,  comp1251_c0_seq1,  comp1322_c0_seq1,  comp1738_c0_seq1,  comp2681_c0_seq1,  comp2806_c0_seq1,  comp2852_c0_seq1,  comp2866_c0_seq1,  comp2870_c0_seq1,  comp2955_c0_seq1,  comp3495_c0_seq1,  comp4962_c0_seq1,  comp4991_c0_seq1,  comp4997_c0_seq1,  comp5085_c0_seq1,  comp5094_c0_seq1,  comp5115_c0_seq1,  comp5186_c0_seq1,  comp5201_c0_seq1,  comp5225_c0_seq1,  comp5237_c0_seq1,  comp5270_c0_seq1,  comp5317_c0_seq1,  comp5351_c0_seq1,  comp5402_c0_seq1,  comp5770_c0_seq1,  comp5963_c0_seq1,  comp935_c0_seq1
     GO:0005975  0.000281178108282779     0.999895382246112         34          72        carbohydrate metabolic process               BP                  0.474138636506016  BP                   carbohydrate metabolic process            comp1012_c0_seq1,  comp1017_c0_seq1,  comp1172_c0_seq1,  comp1249_c0_seq1,  comp1588_c0_seq1,  comp1738_c0_seq1,  comp2805_c0_seq1,  comp2806_c0_seq1,  comp2840_c0_seq1,  comp2852_c0_seq1,  comp2870_c0_seq1,  comp3495_c0_seq1,  comp4027_c0_seq2,  comp4900_c0_seq1,  comp4985_c0_seq6,  comp4991_c0_seq1,  comp5073_c0_seq2,  comp5073_c0_seq4,  comp5115_c0_seq1,  comp5150_c0_seq1,  comp5186_c0_seq1,  comp5187_c0_seq1,  comp5191_c0_seq1,  comp5211_c0_seq1,  comp5222_c0_seq1,  comp5290_c0_seq1,  comp5351_c0_seq1,  comp5455_c0_seq1,  comp5477_c0_seq1,  comp5495_c0_seq1,  comp5575_c0_seq1,  comp5963_c0_seq1,  comp6074_c0_seq1,  comp7332_c0_seq1
     GO:0044723  0.000435275555708771     0.999856357770813         27          51        single-organism      carbohydrate metabolic process             BP                 0.474138636506016    BP                 single-organism    carbohydrate       metabolic          process            comp1012_c0_seq1,  comp1017_c0_seq1,  comp1172_c0_seq1,  comp1249_c0_seq1,  comp1588_c0_seq1,  comp1738_c0_seq1,  comp2805_c0_seq1,  comp2806_c0_seq1,  comp2840_c0_seq1,  comp2852_c0_seq1,  comp2870_c0_seq1,  comp3495_c0_seq1,  comp4027_c0_seq2,  comp4900_c0_seq1,  comp4991_c0_seq1,  comp5073_c0_seq2,  comp5073_c0_seq4,  comp5186_c0_seq1,  comp5211_c0_seq1,  comp5222_c0_seq1,  comp5290_c0_seq1,  comp5351_c0_seq1,  comp5455_c0_seq1,  comp5477_c0_seq1,  comp5575_c0_seq1,  comp5963_c0_seq1,  comp6074_c0_seq1
     GO:0005355  0.00109192623555633      1                         7           7         glucose transmembrane transporter activity            MF                 0.474138636506016    MF                 glucose transmembrane transporter activity           comp5073_c0_seq2,  comp5073_c0_seq3,  comp5073_c0_seq4,  comp5075_c0_seq1,  comp5075_c0_seq3,  comp5075_c1_seq1,  comp5075_c2_seq1
     GO:0008645  0.00109192623555633      1                         7           7         hexose transport      BP                    0.474138636506016   BP                 hexose transport          comp5073_c0_seq2,  comp5073_c0_seq3,  comp5073_c0_seq4,  comp5075_c0_seq1,  comp5075_c0_seq3,  comp5075_c1_seq1,  comp5075_c2_seq1


## Example Data

See 'TRINITY_HOME/sample_data/test_GOSeq_trinotate_pipe/Spombe_runGOseqOnly/'
which includes sample input data and a 'runMe.sh' script for execution.

