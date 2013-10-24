#!/bin/bash 

./cleanMe.pl

echo
echo
echo "*************************************************"
echo "running AAT nucleotide spliced alignment pipeline"
echo "*************************************************"
echo 
echo

cmd="../bin/AAT.pl -X -N -q arab.genomicSeq -s arab.cdna --dds '-f 100 -i 20 -o 75 -p 70 -a 2000' --filter '-c 10' --gap2 '-x 1' "

echo $cmd
eval $cmd


echo "**********************************************"
echo "running AAT protein spliced alignment pipeline"
echo "**********************************************"
echo
echo

 
cmd="../bin/AAT.pl -X -P -q arab.genomicSeq -s arab.pep --dps '-f 100 -i 30 -a 200' --filter '-c 10' --nap '-x 10' "
echo $cmd
eval $cmd

echo "************************************************"
echo "combining results into a multiple alignment file"
echo "************************************************"
echo
echo

cmd="../bin/show *gap2 *nap > multalignment.show.txt"
echo $cmd
eval $cmd

echo
echo
echo "done.  See 'multalignment.show.txt' to examine all spliced alignments."
echo


## run exonerate based on the filter output for comparison  (Thanks, Alexie!)

exonerate_util/run_exonerate.pl -in arab.pep -ref arab.genomicSeq -thread 4 -same -protein -aat *filter -separ
cat arab.pep_queries/*results > exonerate.result
exonerate_util/get_golden_exonerate.pl -exon exonerate.result -genome arab.genomicSeq

