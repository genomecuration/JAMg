#!/bin/bash

###########################
# CHECK THESE
MAX_CONTIGS=3 #NB if there are more than 9999 scaffold, then change the number
PREFIX_GENOME="EG6_Chr" # after prefix a number is experted, 0 to $MAX_CONTIGS
###########################

# pass genome and hints files as arguments
FASTA_GENOME=$1  # genome fasta
HINTS=$2   # merged hints file (full path or relative to current dir)
AUGUSTUS_EXTRINSIC_CFG=$3
RESULTS=$4

MODE=complete  # partial or complete
UTR=off   # on or off
ALT_FROM_EVIDENCE=false
#set -beEu -o pipefail

if [ ! $SPECIES ]; then
	echo "It seems that the environmental variable SPECIES is not set. You must set this environmental variables as per JAMG documentation."
	exit
fi

if [[ ! $FASTA_GENOME || ! -s $FASTA_GENOME ]]; then 
	echo "Please provide the genome fasta, hints file and extrinsic cfg file (in that order)"
	exit
fi

if [[ ! $HINTS || ! -s $HINTS ]]; then 
	echo "Please provide the genome fasta, hints file and extrinsic cfg file (in that order)"
	exit
fi

if [[ ! $AUGUSTUS_EXTRINSIC_CFG || ! -s $AUGUSTUS_EXTRINSIC_CFG ]]; then 
	echo "Please provide the genome fasta, hints file and extrinsic cfg file (in that order)"
	exit
fi

if [ ! $RESULTS ]; then
    RESULTS=results
fi

SOURCE="$(dirname "$(test -L "$0" && readlink "$0" || echo "$0")")"
export PATH=$SOURCE:$JAMG_PATH:$PATH

if [ ! -d $FASTA_GENOME.split ]; then
 splitfasta.pl -i $FASTA_GENOME -depth 1 -dir $FASTA_GENOME.split
fi



echo preparing genome commands
mkdir -p $RESULTS
rm -f $RESULTS.cmds

for (( i=1; i<=$MAX_CONTIGS; i++ ));
  do
  if [ -e $FASTA_GENOME.split/$PREFIX_GENOME$i ]; then
    grep "^$PREFIX_GENOME$i\b" $HINTS > $RESULTS/$PREFIX_GENOME$i.hints
    echo "augustus --UTR=$UTR --gff3=on --species=$SPECIES --uniqueGeneId=true --genemodel=$MODE  --alternatives-from-evidence=$ALT_FROM_EVIDENCE --maxtracks=10 --extrinsicCfgFile=$AUGUSTUS_EXTRINSIC_CFG  --hintsfile=$RESULTS/$PREFIX_GENOME$i.hints $FASTA_GENOME.split/$PREFIX_GENOME$i \
   2> $RESULTS/$PREFIX_GENOME$i.augustus.log > $RESULTS/$PREFIX_GENOME$i.augustus.result & " >> $RESULTS.cmds
  fi
done;

chmod +x $RESULTS.cmds
# if shuf exists
#shuf augustus.commands > augustus.commands.
#if [ -s augustus.commands. ]; then mv -f augustus.commands. augustus.commands; fi
#split -a 3 -d -l 400 augustus.commands augustus.commands.

echo "
 Done, see $RESULTS.cmds
 What to do next? 
 #0 move to compute node and:
 ## import any environmental variables
 ##    export AUGUSTUS_CONFIG_PATH=\$JAMG_PATH/3rd_party/augustus/config
 #1 run with Parafly (either augustus.commands or the split files augustus.commands.* for each compute node)
 #2 grep -hv '^#' $RESULTS/*result | grep -h '\bAUGUSTUS\b' > augustus_results.gtf
 #3 augustus_gtf_2_gtf_proper.pl augustus_results.gtf
 #4 gtf_to_gff3_format.pl augustus_results.gtf.gtf > augustus_results.gff3
"
