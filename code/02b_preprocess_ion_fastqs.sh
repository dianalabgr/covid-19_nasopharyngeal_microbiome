#!/usr/bin/env bash
# example run: ./preprocess_ion_fastqs.sh ./ FASTQ_FILTERED

RUN_PATH=$1
OUTPUT_PATH=$PWD/$2
mkdir $OUTPUT_PATH
echo $RUN_PATH
cd $RUN_PATH

for file in $(ls |grep '.fastq')
do
  id=`basename --suffix .fastq $file`
  echo $id
  cutadapt --cores=24 \
  -q 25 \
  --length 300 \
  --minimum-length 100 \
  --output=$OUTPUT_PATH/$id"_filt.fastq" \
  $file	
done

