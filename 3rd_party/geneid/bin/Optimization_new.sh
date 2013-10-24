#!/bin/bash
# GNU bash, version 2.03.6(1)-release (i386-redhat-linux-gnu)
# $Id: tetraodon.nw,v 0.2 2001/07/13 14:01:56 gparra Exp gparra $
#
SECONDS=0 # Reset Timing
# Which script are we running...
L="####################"
{ echo "$L$L$L$L";
  echo "### RUNNING [$0]";
  echo "### Current date:`date`";
  echo "###"; } 1>&2;
# This program adds a score to a specified type of exons, runs geneid and
# evaluates the results with the annotations.

IeWF=-4.5                              # Initial Weigth First     value -7
deWF=0.5                               # Delta Weigth First       value 0.5
FeWF=-1.5                              # Final Weigth First       value -1.5
IeWFini=$IeWF                          # Initial value

IoWF=0.20                              # Initial Weigth First     value 0.30
doWF=0.05                              # Delta Weigth First       value 0.05
FoWF=0.45                              # Final Weigth First       value 0.70


GENEID=$HOME/software/geneid/bin/
EVAL=$HOME/software/geneid/Evaluation/bin/

PARAM=$1     # Param file for geneid
SEQFILE=$2
CDSFILE=$3

if [[ ! $PARAM || ! $SEQFILE || ! $CDSFILE || ! -f $PARAM || ! -f $SEQFILE || ! -f $CDSFILE ]];then
	echo No input! Give: parameter_file genome_fasta gff_annotations
	exit 255
fi

echo "   oWF    eWF    SN     SP     CC    SNe    SPe   SNSP    raME   raWE \
  raMG    raWG"

ef=`gawk "BEGIN{print ($IeWF <= $FeWF) ? 1 : 0}"`

of=`gawk "BEGIN{print ($IoWF <= $FoWF) ? 1 : 0}"`


while [ $of -ne 0 ]
do
        IeWF=$IeWFini
        ef=`gawk "BEGIN{print ($IeWF <= $FeWF) ? 1 : 0}"`

        while [ $ef -ne 0 ]
        do

              # update exon score in geneid.param.multi

             gawk "{if (NR==27)  {print 1-$IoWF+0.2, 1-$IoWF-0.1, 1-$IoWF-0.1, 1-$IoWF+0.2;} else if (NR==30)  {print $IoWF, $IoWF, $IoWF, $IoWF;} else if (NR==36) {print $IeWF, $IeWF, $IeWF, $IeWF;} else  print}" $PARAM > /tmp/tmp.$$

              # make gene predictions $BIN/geneid/bin/geneid

            $GENEID/geneid -GP /tmp/tmp.$$ $SEQFILE |egrep -v 'exon|^#' |sort +0 -1 +3n | gawk '{if (NR==1) ant=$1; if ($1!=ant) {print "#$";ant=$1}; print }' > /tmp/Predictions.$$.gff
            
              # evaluate gene predictions

             $EVAL/evaluation -sta /tmp/Predictions.$$.gff $CDSFILE | tail -2 | head -1 |  gawk "{printf \"%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\\n\", $IoWF, $IeWF, \$1, \$2, \$3, \$4, \$5, \$6, \$7, \$8, \$12, \$13}"

              ## remove files

              rm /tmp/Predictions.$$.gff

          IeWF=`echo $IeWF $deWF | gawk '{print $1+$2}'`
          ef=`gawk "BEGIN{print ($IeWF <= $FeWF) ? 1 : 0}"`
        done

    IoWF=`echo $IoWF $doWF | gawk '{print $1+$2}'`
    of=`gawk "BEGIN{print ($IoWF <= $FoWF) ? 1 : 0}"`
  done

rm /tmp/tmp.$$
{ echo "###"; echo "### Execution time for [$0] : $SECONDS secs";
  echo "$L$L$L$L";
  echo ""; } 1>&2;
#
exit 0
