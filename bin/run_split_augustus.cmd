#!/bin/bash

###########################
# CHECK THESE
MAX_CONTIGS=9999 #NB if there are more than 9999 scaffold, then change the number.
 # Doesn't matter if it is too high (script will be slightly slower)
PREFIX_GENOME="scaffold_" # after prefix a number is experted, 0 to $MAX_CONTIGS
###########################

# pass genome and hints files as arguments
FASTA_GENOME=$1  # genome fasta
HINTS=$2   # merged hints file (full path or relative to current dir)
AUGUSTUS_EXTRINSIC_CFG=$3


MODE=complete  # partial or complete
UTR=on   # on or off


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

SOURCE="$(dirname "$(test -L "$0" && readlink "$0" || echo "$0")")"
export PATH=$SOURCE:$JAMG_PATH/bin/:$PATH

if [ ! -d $FASTA_GENOME.split ]; then
 splitfasta.pl -i $FASTA_GENOME -depth 1 -dir $FASTA_GENOME.split
fi

echo preparing HINTS
split_augustus_hints.pl $HINTS

echo preparing genome commands
mkdir results
rm -f augustus.commands augustus.commands.???

for (( i=0; i<=$MAX_CONTIGS; i++ ));
  do
  if [ -e $FASTA_GENOME.split/$PREFIX_GENOME$i ]; then
    touch hints/$PREFIX_GENOME$i.hints
    echo "augustus --UTR=$UTR --species=$SPECIES --uniqueGeneId=true --genemodel=$MODE  --alternatives-from-evidence=true --maxtracks=10 --extrinsicCfgFile=$AUGUSTUS_EXTRINSIC_CFG  --hintsfile=hints/$PREFIX_GENOME$i.hints $FASTA_GENOME.split/$PREFIX_GENOME$i \
   2> results/$PREFIX_GENOME$i.augustus.log > results/$PREFIX_GENOME$i.augustus.result " >> augustus.commands
  fi
done;

# if shuf exists
shuf augustus.commands > augustus.commands.
if [ -s augustus.commands. ]; then mv -f augustus.commands. augustus.commands; fi
split -a 3 -d -l 400 augustus.commands augustus.commands.

echo "
 Done, see augustus.commands or augustus.commands.*
 What to do next? 
 #0 move to compute node (optionally) and:
 ## import any environmental variables
 ##    export AUGUSTUS_CONFIG_PATH=\$JAMG_PATH/3rd_party/augustus/config
 #1 run with Parafly (either augustus.commands or the split files augustus.commands.* for each compute node)
 ##    \$JAMG_PATH/3rd_party/bin/ParaFly -c augustus.commands -CPU \$LOCAL_CPUS -v
 #2 grep -hv '^#' results/*result | grep -h '\bAUGUSTUS\b' > augustus_results.gtf
 #3 augustus_gtf_2_gtf_proper.pl augustus_results.gtf
 # ignore the errors regarding CDS length (an Augustus issue) here:
 #4 gtf_to_gff3_format.pl augustus_results.gtf.gtf $FASTA_GENOME > augustus_results.gff3 2> err
 ## Pls see instructions.txt for this instructions again
"  > instructions.txt

cat instructions.txt
