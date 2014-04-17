#!/bin/bash

###########################
# CHECK THESE
MAX_CONTIGS=9999 #NB there are more than 9999 scaffold, then change the number
PREFIX_GENOME="scaffold_"
###########################

# pass genome and hints files as arguments
FASTA_GENOME=$1  # genome fasta
HINTS=$2   # merged hints file (full path or relative to current dir)
AUGUSTUS_EXTRINSIC_CFG=$3


MODE=complete  # partial or complete
UTR=on   # on or off


if [[ ! $FASTA_GENOME || ! -s $FASTA_GENOME ]]; then 
	echo Please provide the genome fasta, hints file and extrinsic cfg file (in that order)
	exit
fi

if [[ ! $HINTS || ! -s $HINTS ]]; then 
	echo Please provide the genome fasta, hints file and extrinsic cfg file (in that order)
	exit
fi

if [[ ! $AUGUSTUS_EXTRINSIC_CFG || ! -s $AUGUSTUS_EXTRINSIC_CFG ]]; then 
	echo Please provide the genome fasta, hints file and extrinsic cfg file (in that order)
	exit
fi

SOURCE="$(dirname "$(test -L "$0" && readlink "$0" || echo "$0")")"
export PATH=$SOURCE:$JAMG_PATH:$PATH

if [ ! -d $FASTA_GENOME.split ]; then
 splitfasta.pl -i $FASTA_GENOME -depth 1 -dir $FASTA_GENOME.split
fi

echo preparing HINTS
split_augustus_hints.pl $HINTS

echo preparing genome commands
mkdir results
rm -f augustus.commands

for (( i=0; i<=$MAX_CONTIGS; i++ ));
  do
  if [ -e *_dir1/$PREFIX_GENOME$i ]; then
    echo "augustus --UTR=$UTR --species=H.armigera --uniqueGeneId=true --genemodel=$MODE  --alternatives-from-evidence=true --maxtracks=10 --extrinsicCfgFile=$AUGUSTUS_EXTRINSIC_CFG  --hintsfile=hints/$PREFIX_GENOME$i.hints $DIR/$PREFIX_GENOME$i \
   2> results/$PREFIX_GENOME$i.augustus.log > results/$PREFIX_GENOME$i.augustus.result; grep -v '^#' results/$PREFIX_GENOME$i.augustus.result|grep '\bAUGUSTUS\b' > results/$PREFIX_GENOME$i.result.gtf " >> augustus.commands
  fi
done;

# if shuf exists
shuf augustus.commands > augustus.commands.
if [ -s augustus.commands. ]; then mv -f augustus.commands. augustus.commands; fi
split -a 3 -d -l 400 augustus.commands augustus.commands.

# What to do next?
# move to compute node and:
#1 run with Parafly
#2 then cat results/*gtf > augustus_results.gtf
#3 augustus_gtf_2_gtf_proper.pl augustus_results.gtf
#4 gtf_to_gff3_format.pl augustus_results.gtf.gtf > augustus_results.gff3
