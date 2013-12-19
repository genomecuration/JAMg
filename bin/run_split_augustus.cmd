#!/bin/bash

# pass genome and hints files as arguments
GENOME=$1  # genome fasta
HINTS=$2   # merged hints file (full path or relative to current dir)
MODE=complete  # partial or complete
UTR=on   # on or off


if [[ ! $GENOME || ! -s $GENOME ]]; then 
	echo provide the genome fasta and hints file
	exit
fi

if [[ ! $HINTS || ! -s $HINTS ]]; then 
	echo provide the genome fasta and hints file
	exit
fi

SOURCE="$(dirname "$(test -L "$0" && readlink "$0" || echo "$0")")"
export PATH=$SOURCE:$PATH

if [ ! -d $GENOME.split ]; then
 splitfasta.pl -i $GENOME -depth 1 -dir $GENOME.split
fi

echo preparing HINTS
split_augustus.pl $HINTS

echo preparing genome commands
mkdir results
rm -f augustus.commands

#NB there are more than 9999 scaffold, then change the number below (but why would you have 9999 scaffolds?!)
for i in {0..9999};
  do
  if [ -e *_dir1/scaffold_$i ]; then
    echo "augustus --UTR=$UTR --species=H.armigera --uniqueGeneId=true --genemodel=$MODE  --alternatives-from-evidence=true --maxtracks=10 --extrinsicCfgFile=extrinsic.all.cfg  --hintsfile=hints/scaffold_$i.hints $DIR/scaffold_$i \
   2> results/scaffold_$i.augustus.log > results/scaffold_$i.augustus.result; grep -v '^#' results/scaffold_$i.augustus.result|grep '\bAUGUSTUS\b' > results/scaffold_$i.result.gtf " >> augustus.commands
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
