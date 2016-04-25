#!/bin/bash
./cleanme.sh

printf "Test 1: Doing straight trimming\n" |tee -a errors
../preprocess_illumina.pl -paired -noadaptor test?.fastq.bz2 2>> errors

rm *trim* *zip -f
printf "\n\nTest 2: Using also Trimmomatic and adaptor searching\n" |tee -a errors
../preprocess_illumina.pl -paired test?.fastq.bz2 2>>errors

rm *trim* *zip -f
printf "\n\nTest 3: Doing straight trimming with Allpaths deduplication\n" |tee -a errors
../preprocess_illumina.pl -paired -dedup test -noadaptor test?.fastq.bz2 2>> errors

rm *trim* *zip -f
printf "\n\nTest 4: Doing straight trimming with Native deduplication\n" |tee -a errors
../preprocess_illumina.pl -paired -dedup 1 -noadaptor test?.fastq.bz2 2>> errors

printf "\n\nDONE! See errors file for any errors.\n"
