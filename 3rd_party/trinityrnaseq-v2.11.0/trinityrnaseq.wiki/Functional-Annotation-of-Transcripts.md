# Functional Annotation of Trinity Transcriptome Assemblies

We developed a functional annotation protocol and supporting software for functionally annotating Trinity de novo transcriptome assemblies called [Trinotate](http://trinotate.github.io), available on GitHub at <http://trinotate.github.io>.  Visit the [Trinotate website](http://trinotate.github.io) for documentation and software.

## Add Annotations to Expression Matrices

It's often useful to include a brief annotation string along with the transcript (or gene) identifier so that it's carried through all downstream expression analyses, such as included in heatmaps, etc.   One way to do this is to take your count matrix and update the feature ID (left-most column of the matrix) to include annotation identifiers. Using the Trinotate (see above) report file, you can do the following to integrate functional annotations:

>Scripts for doing this are in either the Trinotate or Trinity package, so be sure to take notice. Script names are prefixed with either ${TRINOTATE_HOME} or ${TRINITY_HOME} accordingly.

First, using the Trinotate report file, generate a map of feature identifier to an annotated feature identifier like so:

    ${TRINOTATE_HOME}/util/Trinotate_get_feature_name_encoding_attributes.pl \
                      Trinotate_report.xls  > annot_feature_map.txt

Looking at the first few lines of the output file, you'll see a formatting that resembles the following:

    comp0_c0	comp0_c0^YFY8_SCHPO^Steroid_dh^Tm3
    comp0_c0_seq1	comp0_c0_seq1^YFY8_SCHPO^Steroid_dh^Tm3
    comp1000_c0	comp1000_c0^ARD1_SCHPO^Acetyltransf_1
    comp1000_c0_seq1	comp1000_c0_seq1^ARD1_SCHPO^Acetyltransf_1
    comp10023_c0	comp10023_c0^TF29_SCHPO^Peptidase_A2E
    comp10023_c0_seq1	comp10023_c0_seq1^TF29_SCHPO^Peptidase_A2E
    comp1002_c0	comp1002_c0^CCS1_SCHPO
    comp1002_c0_seq1	comp1002_c0_seq1^CCS1_SCHPO
    comp1002_c0_seq2	comp1002_c0_seq2^CCS1_SCHPO
    comp10033_c0	comp10033_c0^RSM1_SCHPO^zf-C3HC

The above maps the original feature identifier (Trinity transcript or gene identifier) to a version that includes the accessions of top blast hits and Pfam identifiers, transmembrane domains, and and signal peptides.

Given this annotation mapping file, you can then update your expression matrix. For example, given a counts matrix 'Trinity_trans_counts.matrix', which looks like so:

              ds_rep1 hs_rep1 log_rep1        plat_rep1
    comp4806_c1_seq1        0.00    1.00    0.00    2.00
    comp3938_c1_seq1        6.00    4.00    9.00    4.00
    comp31964_c0_seq1       0.00    0.00    0.00    0.00
    comp6295_c0_seq1        95.00   77.00   95.00   86.00
    comp3109_c1_seq1        22.00   8.00    28.00   47.00
    comp3439_c1_seq1        4.00    1.00    3.00    6.00

we can integrate functional annotations like so:

    ${TRINITY_HOME}/Analysis/DifferentialExpression/rename_matrix_feature_identifiers.pl \
        Trinity_trans.counts.matrix annot_feature_map.txt > Trinity_trans.counts.wAnnot.matrix

and looking at the top entries of our updated matrix:

                ds_rep1	hs_rep1	log_rep1	plat_rep1
    comp4806_c1_seq1	0.00	1.00	0.00	2.00
    comp3938_c1_seq1^POPI_SCHPO^POP1	6.00	4.00	9.00	4.00
    comp31964_c0_seq1^NAA25_SCHPO	0.00	0.00	0.00	0.00
    comp6295_c0_seq1^NAGS_SCHPO^DUF619	95.00	77.00	95.00	86.00
    comp3109_c1_seq1^DEP1_SCHPO^Sds3	22.00	8.00	28.00	47.00
    comp3439_c1_seq1^CUT12_SCHPO	4.00	1.00	3.00	6.00
    comp4603_c2_seq1^MAP4_SCHPO	1.00	0.00	0.00	0.00
    comp596_c0_seq1^ELL1_SCHPO^ELL	174.00	81.00	58.00	203.01
    comp4091_c0_seq1^BYR1_SCHPO^Pkinase_Tyr	139.00	51.00	66.00	203.01


This updated matrix can be used for downstream analysis steps, such as differential expression analysis.

