conda activate qiime2_2021.02
source tab-qiime
#conda info
#conda install -c conda-forge deicode


PROJECT_FOLDER="/mnt/raid1/philipos/corona_project/16S_analysis"
QIIME_OUTDIR=$PROJECT_FOLDER/data/intermediate/02_qiime_denoise
IMPORT_LIST=$PROJECT_FOLDER/data/intermediate/01_quality_filtered_fastq_files/filtered_fastq_list.tsv
METADATA_FILE=$PROJECT_FOLDER/data/16S_3batches_metadata.tsv
GREENGENES_FOLDER=$PROJECT_FOLDER/data/reference/gg_13_8_otus
cd $PROJECT_FOLDER

#STEP 1 ubam to fastq
./code/02a_run_bamtofastq_on_many_files.sh "data/raw/first_batch/bam_files_22022021" "data/raw/first_batch/fastq"
./code/02a_run_bamtofastq_on_many_files.sh "data/raw/second_batch/bam_files_10052021" "data/raw/second_batch/fastq"
./code/02a_run_bamtofastq_on_many_files.sh "data/raw/third_batch/bam_files_05112021" "data/raw/third_batch/fastq"


#STEP 2 preprocess fastq files
./code/02b_preprocess_ion_fastqs.sh data/raw/first_batch/fastq data/intermediate/01_quality_filtered_fastq_files/first_batch
./code/02b_preprocess_ion_fastqs.sh data/raw/second_batch/fastq data/intermediate/01_quality_filtered_fastq_files/second_batch
./code/02b_preprocess_ion_fastqs.sh data/raw/third_batch/fastq data/intermediate/01_quality_filtered_fastq_files/third_batch



#STEP 3 QIIME2  DENOISING & BASIC ANALYSIS
#import to qiime2
cd $QIIME_OUTDIR
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path $IMPORT_LIST \
  --output-path demux.qza \
  --input-format SingleEndFastqManifestPhred33V2
qiime demux summarize   --i-data demux.qza   --o-visualization demux.qzv

#denoise and find features (group reads) but use pyro and trim 15 left
qiime dada2 denoise-pyro \
  --i-demultiplexed-seqs demux.qza \
  --p-trim-left 15 \
  --p-trunc-len 0 \
  --o-representative-sequences rep-seqs-dada2.qza \
  --o-table table-dada2.qza \
  --o-denoising-stats stats-dada2.qza \
  --p-n-threads 40

qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o-visualization stats-dada2.qzv

mv rep-seqs-dada2.qza rep-seqs.qza
mv table-dada2.qza table.qza

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime feature-table summarize \
  --i-table table.qza \
  --m-sample-metadata-file $METADATA_FILE \
  --o-visualization table.qzv

