#!/usr/bin/env bash
# example run: ./run_bamtofastq_on_many_files.sh RAW_DATA/16S_metagenomics_experiment_05112021/ubam_05112021/ ./FASTQ

START_PATH=$PWD
RUN_PATH=$1
OUTPUT_PATH=$2
mkdir $OUTPUT_PATH
echo $RUN_PATH
cd $RUN_PATH

for file in $(ls |grep '.bam')
do
  id=`basename --suffix .bam $file`
  echo $id
  bamToFastq -i $file \
  -fq $START_PATH/$OUTPUT_PATH/$id.fastq
done