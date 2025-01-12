conda activate qiime2_2021.02
source tab-qiime
#conda info
#conda install -c conda-forge deicode

PROJECT_FOLDER="/mnt/raid1/philipos/corona_project/16S_analysis"
OUTDIR=$PROJECT_FOLDER/data/intermediate/03_qiime_otu_clustering
GREENGENES_FOLDER=$PROJECT_FOLDER/data/reference/gg_13_8_otus
METADATA_FILE=$PROJECT_FOLDER/data/16S_3batches_metadata.tsv

#Closed-reference clustering (Classification)
qiime vsearch cluster-features-closed-reference \
     --i-table $PROJECT_FOLDER/data/intermediate/02_qiime_denoise/table.qza \
     --i-sequences $PROJECT_FOLDER/data/intermediate/02_qiime_denoise/rep-seqs.qza \
     --i-reference-sequences $GREENGENES_FOLDER/99_otus.qza \
     --p-perc-identity 0.99 \
	 --p-threads 20 \
	 --p-strand 'both' \
     --o-clustered-table $OUTDIR/table-cr-99.qza \
     --o-clustered-sequences $OUTDIR/rep-seqs-cr-99.qza \
     --o-unmatched-sequences $OUTDIR/unmatched-cr-99.qza
	 
qiime feature-table tabulate-seqs \
  --i-data $OUTDIR/rep-seqs-cr-99.qza \
  --o-visualization $OUTDIR/rep-seqs-cr-99.qzv
#qiime tools view $OUTDIR/rep-seqs-cr-99.qzv

qiime feature-table summarize \
  --i-table $OUTDIR/table-cr-99.qza \
  --m-sample-metadata-file $METADATA_FILE \
  --o-visualization $OUTDIR/table-cr-99.qzv
#qiime tools view $OUTDIR/table-cr-99.qzv

#visualize taxonomy
qiime metadata tabulate \
  --m-input-file $GREENGENES_FOLDER/99_otus_taxonomy.qza \
  --m-input-file $OUTDIR/rep-seqs-cr-99.qza \
  --o-visualization $GREENGENES_FOLDER/99_otus_taxonomy.qzv
