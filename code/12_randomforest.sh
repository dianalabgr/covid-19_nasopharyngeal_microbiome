conda activate qiime2_2021.02
source tab-qiime
#conda info
#conda install -c conda-forge deicode

PROJECT_FOLDER="/mnt/raid1/philipos/corona_project/16S_analysis"
OUTDIR=$PROJECT_FOLDER/data/intermediate/10_randomforest
OUTDIR2=$PROJECT_FOLDER/data/intermediate/11_hierarchical_clustering
mkdir $OUTDIR2
INDIR=$PROJECT_FOLDER/data/intermediate/06_frequency_filtering
METADATA_FILE=$PROJECT_FOLDER/data/16S_3batches_metadata.tsv
mkdir $OUTDIR
#use genus level only
qiime sample-classifier classify-samples \
      --i-table $INDIR/filtered_table_l6.qza \
      --m-metadata-file $METADATA_FILE \
      --m-metadata-column group_1 \
      --p-optimize-feature-selection \
      --p-parameter-tuning \
      --p-estimator RandomForestClassifier \
      --p-n-estimators 100 \
      --p-n-jobs 40 \
      --p-random-state 1 \
      --output-dir $OUTDIR/classifier_1_l6
qiime tools view $OUTDIR/classifier_1_l6/accuracy_results.qzv  
  
qiime metadata tabulate \
  --m-input-file $OUTDIR/classifier_1_l6/feature_importance.qza \
  --o-visualization $OUTDIR/classifier_1_l6/feature_importance.qzv
#qiime tools view $OUTDIR/classifier_1_l6/feature_importance.qzv 
#qiime tools view $OUTDIR/classifier_1_l6/accuracy_results.qzv
#qiime tools view $OUTDIR/classifier_1_l6/roc_plot.qzv


qiime feature-table summarize \
  --i-table $INDIR/filtered_table_l6.qza \
  --m-sample-metadata-file $METADATA_FILE \
  --o-visualization $INDIR/filtered_table_l6.qzv
#qiime tools view $INDIR/filtered_table_l6.qzv
  
qiime tools export \
  --input-path feature-table.qza \
  --output-path exported-feature-table
  
#########################
#Make a heatmap
qiime feature-table heatmap \
	--i-table $INDIR/samples_features_filtered_table_l6.qza \
	--o-visualization $OUTDIR2/heatmap_l6.qzv
