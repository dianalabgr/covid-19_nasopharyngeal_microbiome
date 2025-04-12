conda activate qiime2_2021.02
source tab-qiime
#conda info

PROJECT_FOLDER="/mnt/raid1/philipos/corona_project/16S_analysis"
OUTDIR=$PROJECT_FOLDER/data/intermediate/06_frequency_filtering
GREENGENES_FOLDER=$PROJECT_FOLDER/data/reference/gg_13_8_otus
METADATA_FILE=$PROJECT_FOLDER/data/covid19_study_metadata.tsv



#This filter can be applied to the feature axis to remove low abundance features from a table. 
#here we remove all features with a total abundance (summed across all samples) of less than 2 (remove singletons).
qiime feature-table filter-features \
 --i-table $PROJECT_FOLDER/data/intermediate/03_qiime_otu_clustering/table-cr-99.qza \
 --p-min-frequency 2 \
 --o-filtered-table $OUTDIR/features_filtered_table.qza
#Visualize
qiime feature-table summarize \
  --i-table $OUTDIR/features_filtered_table.qza \
  --m-sample-metadata-file $METADATA_FILE \
  --o-visualization $OUTDIR/features_filtered_table.qzv
  
#qiime tools view $OUTDIR/features_filtered_table.qzv


#how to remove low count samples
#First, see alpha-rarefaction plots to decide rarefaction/sampling depth
qiime diversity alpha-rarefaction \
	--i-table  $OUTDIR/features_filtered_table.qza \
	--i-phylogeny $GREENGENES_FOLDER/99_otus_tree.qza \
	--p-max-depth 95000 \
	--p-steps 100 \
	--m-metadata-file $METADATA_FILE \
	--o-visualization $OUTDIR/alpha-rarefaction.qzv
#qiime tools view alpha-rarefaction.qzv

#Then filter samples by total frequency (total-frequency-based filtering)
#Total-frequency-based filtering is used to filter samples or features based on how frequently they are represented in the feature table.

#filter samples whose total frequency is an outlier in the distribution of sample frequencies. 
#remove samples based on their minimum total frequency (i.e., total number of sequences obtained for the sample). This can be achieved as follows (samples with a total frequency less than 21432 will be filtered).
qiime feature-table filter-samples \
	--i-table $OUTDIR/features_filtered_table.qza  \
	--p-min-frequency 21432 \
	--o-filtered-table $OUTDIR/filtered_table.qza
#Visualize
qiime feature-table summarize \
  --i-table $OUTDIR/filtered_table.qza \
  --m-sample-metadata-file $METADATA_FILE \
  --o-visualization $OUTDIR/filtered_table.qzv

#create collapsed  tables !
#levels: 1 Kingdom, 2 Phylum, 3 Class, 4 Order, 5 Family, 6 Genus, 7 Species
qiime taxa collapse \
  --i-table $OUTDIR/samples_features_filtered_table.qza \
  --i-taxonomy $GREENGENES_FOLDER/99_otus_taxonomy.qza \
  --p-level 2 \
  --o-collapsed-table $OUTDIR/filtered_table_l2.qza

qiime taxa collapse \
  --i-table $OUTDIR/samples_features_filtered_table.qza \
  --i-taxonomy $GREENGENES_FOLDER/99_otus_taxonomy.qza \
  --p-level 5 \
  --o-collapsed-table $OUTDIR/filtered_table_l5.qza
  
qiime taxa collapse \
  --i-table $OUTDIR/samples_features_filtered_table.qza \
  --i-taxonomy $GREENGENES_FOLDER/99_otus_taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table $OUTDIR/filtered_table_l6.qza

qiime tools export \
  --input-path $OUTDIR/filtered_table_l6.qza \
  --output-path $OUTDIR/filtered_table_l6.tsv
  