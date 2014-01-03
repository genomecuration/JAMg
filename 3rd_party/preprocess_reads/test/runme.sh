#!/bin/bash
rm *trim* *zip errors -f

echo "Test 1: Doing straight trimming" |tee -a errors
../preprocess_illumina.pl -paired -gdna -no_screen test?.fastq.bz2 2>> errors


rm *trim* *zip -f
printf "\n\nTest 2: Using also Trimmomatic and adaptor searching\n" |tee -a errors
../preprocess_illumina.pl -paired -gdna -no_screen -trimmo ~/software/Trimmomatic/*jar  test?.fastq.bz2 2>>errors

printf "\n\nDONE! See errors file for any errors.\n"
