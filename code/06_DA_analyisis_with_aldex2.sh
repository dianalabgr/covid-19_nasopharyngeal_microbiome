conda activate qiime2_2021.02
source tab-qiime
#conda info

PROJECT_FOLDER="/mnt/raid1/philipos/corona_project/16S_analysis"
OUTDIR=$PROJECT_FOLDER/data/intermediate/09_DA_analysis
SE_NE=$OUTDIR"/DA_severe_vs_negative"
SE_AS=$OUTDIR"/DA_severe_vs_asymptomatic"
SE_MI=$OUTDIR"/DA_severe_vs_mild"
MI_NE=$OUTDIR"/DA_mild_vs_negative"
MI_AS=$OUTDIR"/DA_mild_vs_asymptomatic"
AS_NE=$OUTDIR"/DA_asymptomatic_vs_negative"

INDIR=$PROJECT_FOLDER/data/intermediate/06_frequency_filtering
GREENGENES_FOLDER=$PROJECT_FOLDER/data/reference/gg_13_8_otus
METADATA_FILE=$PROJECT_FOLDER/data/16S_3batches_metadata.tsv
mkdir $SE_NE
mkdir $SE_AS
mkdir $SE_MI
mkdir $MI_NE
mkdir $MI_AS
mkdir $AS_NE

######################
#COMPARE severity groups
#l2 = phylum level
#l5 = family level
#l6 = genus level

X=mild
Y=asymptomatic
METADATA_COLUMN="type2"
for X in severe mild asymptomatic
do
  for Y in mild asymptomatic negative
  do
	if [[ $X == $Y ]]
	then
	  continue
	elif [[ $X == "asymptomatic" ]] && [[ $Y == "mild" ]]
	then
	  continue
	else
	  #l2
    qiime feature-table filter-samples \
      --i-table $INDIR/filtered_table_l2.qza \
      --m-metadata-file $METADATA_FILE \
      --p-where "[type]='"$X"' OR [type]='"$Y"'" \
      --o-filtered-table $OUTDIR"/DA_"$X"_vs_"$Y"/filtered_table_l2_"$X"_"$Y".qza"
    qiime aldex2 aldex2 \
      --i-table $OUTDIR"/DA_"$X"_vs_"$Y"/filtered_table_l2_"$X"_"$Y".qza" \
      --m-metadata-file $METADATA_FILE \
      --m-metadata-column $METADATA_COLUMN \
      --o-differentials $OUTDIR"/DA_"$X"_vs_"$Y"/differentials_l2.qza"
    qiime aldex2 extract-differences \
      --i-table $OUTDIR"/DA_"$X"_vs_"$Y"/differentials_l2.qza" \
      --o-differentials $OUTDIR"/DA_"$X"_vs_"$Y"/extracted_sig_differences_l2" \
      --p-sig-threshold 1 \
      --p-effect-threshold 0 \
      --p-difference-threshold 0
    #The tab separated file of differentially called features can be exported.
    mkdir $OUTDIR"/DA_"$X"_vs_"$Y"/differentials"
    qiime tools export \
      --input-path $OUTDIR"/DA_"$X"_vs_"$Y"/extracted_sig_differences_l2.qza" \
      --output-path $OUTDIR"/DA_"$X"_vs_"$Y"/differentials"
    mv $OUTDIR"/DA_"$X"_vs_"$Y"/differentials/differentials.tsv" $OUTDIR"/DA_"$X"_vs_"$Y"/differentials/l2_differentials.tsv"
  
    #l5
    qiime feature-table filter-samples \
      --i-table $INDIR/filtered_table_l5.qza \
      --m-metadata-file $METADATA_FILE \
      --p-where "[type]='"$X"' OR [type]='"$Y"'" \
      --o-filtered-table $OUTDIR"/DA_"$X"_vs_"$Y"/filtered_table_l5_"$X"_"$Y".qza"
    qiime aldex2 aldex2 \
      --i-table $OUTDIR"/DA_"$X"_vs_"$Y"/filtered_table_l5_"$X"_"$Y".qza" \
      --m-metadata-file $METADATA_FILE \
      --m-metadata-column $METADATA_COLUMN \
      --o-differentials $OUTDIR"/DA_"$X"_vs_"$Y"/differentials_l5.qza"
    qiime aldex2 extract-differences \
      --i-table $OUTDIR"/DA_"$X"_vs_"$Y"/differentials_l5.qza" \
      --o-differentials $OUTDIR"/DA_"$X"_vs_"$Y"/extracted_sig_differences_l5" \
      --p-sig-threshold 0.05 \
      --p-effect-threshold 0 \
      --p-difference-threshold 0
    #The tab separated file of differentially called features can be exported.
    mkdir $OUTDIR"/DA_"$X"_vs_"$Y"/differentials"
    qiime tools export \
      --input-path $OUTDIR"/DA_"$X"_vs_"$Y"/extracted_sig_differences_l5.qza" \
      --output-path $OUTDIR"/DA_"$X"_vs_"$Y"/differentials"
    mv $OUTDIR"/DA_"$X"_vs_"$Y"/differentials/differentials.tsv" $OUTDIR"/DA_"$X"_vs_"$Y"/differentials/l5_differentials.tsv"
  
    #l6
    qiime feature-table filter-samples \
      --i-table $INDIR/filtered_table_l6.qza \
      --m-metadata-file $METADATA_FILE \
      --p-where "[type]='"$X"' OR [type]='"$Y"'" \
      --o-filtered-table $OUTDIR"/DA_"$X"_vs_"$Y"/filtered_table_l6_"$X"_"$Y".qza"
    qiime aldex2 aldex2 \
      --i-table $OUTDIR"/DA_"$X"_vs_"$Y"/filtered_table_l6_"$X"_"$Y".qza" \
      --m-metadata-file $METADATA_FILE \
      --m-metadata-column $METADATA_COLUMN \
      --o-differentials $OUTDIR"/DA_"$X"_vs_"$Y"/differentials_l6.qza"
    qiime aldex2 extract-differences \
      --i-table $OUTDIR"/DA_"$X"_vs_"$Y"/differentials_l6.qza" \
      --o-differentials $OUTDIR"/DA_"$X"_vs_"$Y"/extracted_sig_differences_l6" \
      --p-sig-threshold 0.05 \
      --p-effect-threshold 0 \
      --p-difference-threshold 0
    #The tab separated file of differentially called features can be exported.
    mkdir $OUTDIR"/DA_"$X"_vs_"$Y"/differentials"
    qiime tools export \
      --input-path $OUTDIR"/DA_"$X"_vs_"$Y"/extracted_sig_differences_l6.qza" \
      --output-path $OUTDIR"/DA_"$X"_vs_"$Y"/differentials"
    mv $OUTDIR"/DA_"$X"_vs_"$Y"/differentials/differentials.tsv" $OUTDIR"/DA_"$X"_vs_"$Y"/differentials/l6_differentials.tsv"
	  
	fi
  done
done
