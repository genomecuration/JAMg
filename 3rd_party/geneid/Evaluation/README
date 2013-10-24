README for evaluation'01: Enrique Blanco.

1. Evaluation:

evaluation is a program to assess the accuracy of gene predictions at nucleotide, exon and gene level. Two input 
gff-files are required: one for predictions and other for annotations (real genes). There are several measures
to print about these comparisons:

* Nucleotide level:
NuR    NuP      TP      TN      FP      FN    SN      SP      AC      CC 

* Exon level:
ExR    ExP     TPe      ME      WE  raME    raWE     SNe     SPe    SNSP
 
* Gene level: 
GeR    GeP     TPg      MG      WG  raMG    raWG     SNg     SPg    SNSPg   raJG    raSG
  
New gene level statistics are about split and join genes while remain classical and already known measures
such as specifity or sensibility. Ratios (ra) are percents in Missing/Wrong exons and genes.

--------------------------

2. Options:

Options are about printing:

        -v: Verbose. Print all messages
        -a: Average. Print average stats (more than 1 sequence)
        -t: Total. Print total stats (more than 1 sequence)
        -s: Short. Print a short output                  

* Short: Short version from the original standard output only with percent values at all the levels.

* Average: Arithmetic mean using the already computed values in the input sequences.

* Total: means as whether we built one artificial sequence with the previous read sequences and measured again
the percents, preserving the absolute values such as TN, TP, FN, FP, ...


NAME
        evaluation - a program to measure gene prediction accuracy
SYNOPSIS
        evaluation [-vats] <Predict_genes> <Real_genes>             

--------------------------

3. Format for real and prediction files:

In both files, there is a token to notice end/begin of data for more than one sequence ($)

-> Annotations:
At the beginning it is necessary one gff-line to put the length of sequence.

For instance:

Annotation filename:
AB009080        ANNOTA  Sequence        1       4875    .       .       .       1
AB009080        ANNOTATED       First   354     1352    1.000000        +       0       1
AB009080        ANNOTATED       Internal        1489    1581    1.000000        +       0       1
AB009080        ANNOTATED       Internal        1669    2485    1.000000        +       0       1
AB009080        ANNOTATED       Terminal        2596    4859    1.000000        +       2       1
$
AB016609        ANNOTA  Sequence        1       2836    .       .       .       1
AB016609        ANNOTATED       First   200     552     1.000000        +       0       1
AB016609        ANNOTATED       Internal        703     1233    1.000000        +       1       1
AB016609        ANNOTATED       Internal        1328    1547    1.000000        +       1       1
AB016609        ANNOTATED       Terminal        1630    2766    1.000000        +       0       1   

Prediction filename:
AB009080        geneid_dd       First   354     1352    184.47  +       0       gene_1  #       AA      1:333
AB009080        geneid_dd       Internal        1489    1581    14.70   +       0       gene_1  #       AA      334:364
AB009080        geneid_dd       Internal        1669    2485    139.27  +       0       gene_1  #       AA      365:637
AB009080        geneid_dd       Terminal        2596    4859    314.64  +       2       gene_1  #       AA      637:1391
$
AB016609        geneid_dd       First   200     552     30.68   +       0       gene_1  #       AA      1:118
AB016609        geneid_dd       Internal        703     1227    47.60   +       1       gene_1  #       AA      118:293
AB016609        geneid_dd       Internal        1328    1547    22.06   +       1       gene_1  #       AA      293:366
AB016609        geneid_dd       Terminal        1630    2766    136.21  +       0       gene_1  #       AA      367:745    

Here one script to insert $ into a gff-file is kindly provided (by G.Parra)

Script to insert "$" between sequences in a gff file with more than one sequence:
sort +0 -1 +3n  input.gff | gawk '{if (NR==1) ant=$1; if ($1!=ant) {print "$";ant=$1}; print }' > output.gff

--------------------------

4. Directories:

src: source code
samples: artifTests/ realTests/ -> sets of filenames to test the program
ancient: old code for evaluation (CGI version)
include: headers
bin: place for binary
objects: place for obj. files


