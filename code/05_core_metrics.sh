conda activate qiime2_2021.02
source tab-qiime
#conda info

PROJECT_FOLDER="/mnt/raid1/philipos/corona_project/16S_analysis"
OUTDIR=$PROJECT_FOLDER/data/intermediate/08_core_metrics
INDIR=$PROJECT_FOLDER/data/intermediate/06_frequency_filtering
GREENGENES_FOLDER=$PROJECT_FOLDER/data/reference/gg_13_8_otus
METADATA_FILE=$PROJECT_FOLDER/data/covid19_study_metadata.tsv
mkdir $OUTDIR

#ALPHA AND BETA DIVERSITY
#decide the sampling depth after exploring this plot:  qiime tools view table.qzv
#we will remove the 4 samples with the lowest features/sequences
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny $GREENGENES_FOLDER/99_otus_tree.qza \
  --i-table $INDIR/filtered_table.qza \
  --p-sampling-depth  21432 \
  --m-metadata-file $METADATA_FILE \
  --output-dir $OUTDIR/core_metrics_phylogenetic \
  --p-n-jobs-or-threads 40

#alpha diversity metrics
##Shannon’s diversity index (a quantitative measure of community richness)
##Observed OTUs (a qualitative measure of community richness)
##Faith’s Phylogenetic Diversity (a qualitiative measure of community richness that incorporates phylogenetic relationships between the features)
##Evenness (or Pielou’s Evenness; a measure of community evenness)
qiime diversity alpha-group-significance \
	--i-alpha-diversity $OUTDIR/core_metrics_phylogenetic/shannon_vector.qza \
	--m-metadata-file $METADATA_FILE \
	--o-visualization $OUTDIR/core_metrics_phylogenetic/shannon-significance.qzv
 #qiime tools view $OUTDIR/core_metrics_phylogenetic/shannon-significance.qzv
 
#also calculates Simpson’s index
#Measures the relative abundance of the different species making up the sample richness
qiime diversity alpha \
	--i-table $INDIR/filtered_table.qza \
	--p-metric simpson \
	--o-alpha-diversity $OUTDIR/core_metrics_phylogenetic/simpson_index_vector.qza
  
#create alpha diversity plots
for a in shannon observed_features faith_pd evenness simpson_index
do
qiime diversity alpha-group-significance \
  --i-alpha-diversity $OUTDIR/core_metrics_phylogenetic/$a"_vector.qza" \
  --m-metadata-file $METADATA_FILE \
  --o-visualization $OUTDIR/core_metrics_phylogenetic/$a"_significance.qzv"
done
#qiime tools view $OUTDIR/core_metrics_phylogenetic/shannon_significance.qzv
#qiime tools view $OUTDIR/core_metrics_phylogenetic/observed_otus-significance.qzv
#qiime tools view $OUTDIR/core_metrics_phylogenetic/faith_pd_significance.qzv
#qiime tools view $OUTDIR/core_metrics_phylogenetic/evenness-significance.qzv
#qiime tools view $OUTDIR/core_metrics_phylogenetic/simpson_index_significance.qzv

#see shannon correlation with numeric metadata:
qiime diversity alpha-correlation \
      --i-alpha-diversity $OUTDIR/core_metrics_phylogenetic/shannon_vector.qza \
      --m-metadata-file $METADATA_FILE \
      --o-visualization $OUTDIR/core_metrics_phylogenetic/shannon_correlation.qzv 
#qiime tools view $OUTDIR/core_metrics_phylogenetic/shannon_correlation.qzv 


#beta diversity metrics and plots , grouped by "type" metadata column
for b in unweighted_unifrac weighted_unifrac jaccard bray_curtis
do
qiime diversity beta-group-significance \
  --i-distance-matrix $OUTDIR/core_metrics_phylogenetic/$b"_distance_matrix.qza" \
  --m-metadata-file $METADATA_FILE \
  --m-metadata-column type \
  --o-visualization $OUTDIR/core_metrics_phylogenetic/$b"_type_significance_beta.qzv" \
  --p-pairwise
done

#qiime tools view core-metrics-results/jaccard_emperor.qzv
#qiime tools view core-metrics-results/bray_curtis_emperor.qzv
#qiime tools view core-metrics-results/unweighted_unifrac_emperor.qzv
#qiime tools view core-metrics-results/weighted_unifrac_emperor.qzv

#qiime tools view core-metrics-results/unweighted-unifrac-type-significance_beta.qzv
#qiime tools view core-metrics-results/weighted-unifrac-type-significance_beta.qzv
#qiime tools view core-metrics-results/jaccard-type-significance_beta.qzv
#qiime tools view core-metrics-results/bray-curtis-type-significance_beta.qzv

