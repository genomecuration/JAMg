#!/bin/bash
# GNU bash, version 2.03.6(1)-release (i386-redhat-linux-gnu)
# $Id: tetraodon.nw,v 0.2 2001/07/13 14:01:56 gparra Exp gparra $
#
SECONDS=0 # Reset Timing
# This program adds a score to a specified type of exons, runs geneid and
# evaluates the results with the annotations.

IeWF='-8.5'                              # Initial Weigth First     value -7
deWF='0.5'                               # Delta Weigth First       value 0.5
FeWF='-1.0'                              # Final Weigth First       value -1.5
IeWFini=$IeWF                          # Initial value

IoWF='0.05'                              # Initial Weigth First     value 0.30
doWF='0.05'                             # Delta Weigth First       value 0.05
FoWF='0.70'                              # Final Weigth First       value 0.70

export PATH=$JAMG_PATH/bin:$JAMG_PATH/3rd_party/geneid/bin/:$JAMG_PATH/3rd_party/geneid/Evaluation/bin/:$PATH

PARAM=$1     # Param file for geneid
SEQFILE=$2
CDSFILE=$3
PREFIX=$4

if [[ ! $PARAM || ! $SEQFILE || ! $CDSFILE || ! -f $PARAM || ! -f $SEQFILE || ! -f $CDSFILE || ! $PREFIX ]];then
	echo No input! Give: parameter_file genome_fasta gff_annotations prefix
	exit 255
fi
mkdir /tmp/$USER 2>/dev/null

echo "
SN  = sensitivity nucleotide level
SP  = specificity nucleotide level
CC  = correlation SN/SP
SNe  = sensitivity exon level
SPe  = specificity exon level
SNSP  = correlation SNe/Spe
raME  = ratio missing exons
raWE  = ratio wrong exons
SNSPg = correlation SNg/SPg
raMG  = ratio missing genes
raWG  = ratio wrong genes

oWF	eWF	SN	SP	CC	SNe	SPe	SNSP	raME	raWE	raMG	raWG" > header

ef=`gawk "BEGIN{print ($IeWF <= $FeWF) ? 1 : 0}"`

of=`gawk "BEGIN{print ($IoWF <= $FoWF) ? 1 : 0}"`

while [ $of -ne 0 ]
do
        IeWF=$IeWFini
        ef=`gawk "BEGIN{print ($IeWF <= $FeWF) ? 1 : 0}"`
        while [ $ef -ne 0 ]
        do

echo "gawk '{if (NR==27)  {print 1-$IoWF+0.2, 1-$IoWF-0.1, 1-$IoWF-0.1, 1-$IoWF+0.2;} else if (NR==30)  {print $IoWF, $IoWF, $IoWF, $IoWF;} else if (NR==36) {print $IeWF, $IeWF, $IeWF, $IeWF;} else  print}' $PARAM > /tmp/$USER/tmp.$$.$IoWF.$IeWF   ;   geneid -GP /tmp/$USER/tmp.$$.$IoWF.$IeWF $SEQFILE |egrep -v 'exon|^#' | gawk '{if (NR==1) ant=\$1; if (\$1!=ant) {print \"#\$\";ant=\$1}; print }' > /tmp/$USER/Predictions.$$.$IoWF.$IeWF.gff; purge.geneid.real.gff.pl $CDSFILE /tmp/$USER/Predictions.$$.$IoWF.$IeWF.gff  ;  evaluation -sta /tmp/$USER/Predictions.$$.$IoWF.$IeWF.gff /tmp/$USER/Predictions.$$.$IoWF.$IeWF.gff.valid | tail -2 | head -1 |  gawk '{printf \"%6.2f	%6.2f	%6.2f	%6.2f	%6.2f	%6.2f	%6.2f	%6.2f	%6.2f	%6.2f	%6.2f	%6.2f\\n\", $IoWF, $IeWF, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$12, \$13}' > "$PREFIX"_result_$IoWF.$IeWF.txt ;   rm -f  /tmp/$USER/Predictions.$$.$IoWF.$IeWF.gff* "

          IeWF=$(echo $IeWF + $deWF|bc)
          ef=`gawk "BEGIN{print ($IeWF <= $FeWF) ? 1 : 0}"`
        done

    IoWF=$(echo $IoWF + $doWF|bc)
    of=`gawk "BEGIN{print ($IoWF <= $FoWF) ? 1 : 0}"`
  done
exit
rm /tmp/$USER/tmp.$$
{ echo "###"; echo "### Execution time for [$0] : $SECONDS secs";
  echo "$L$L$L$L";
  echo ""; } 1>&2;
#
exit 0
