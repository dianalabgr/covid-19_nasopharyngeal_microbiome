conda activate qiime2_2021.02
source tab-qiime
#conda info

PROJECT_FOLDER="/mnt/raid1/philipos/corona_project/16S_analysis"
OUTDIR=$PROJECT_FOLDER/results/DA2_analysis
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
level=l2
METADATA_COLUMN="type2"
for X in severe mild asymptomatic
do
  for Y in mild asymptomatic negative
  do
	for level in l2 l5 l6
	do
		
		if [[ $X == $Y ]]
		then
		  continue
		elif [[ $X == "asymptomatic" ]] && [[ $Y == "mild" ]]
		then
		  continue
		else
		qiime feature-table filter-samples \
		  --i-table $INDIR"/filtered_table_"$level".qza" \
		  --m-metadata-file $METADATA_FILE \
		  --p-where "[type]='"$X"' OR [type]='"$Y"'" \
		  --o-filtered-table $OUTDIR"/DA_"$X"_vs_"$Y"/filtered_table_"$level"_"$X"_"$Y".qza"
		qiime aldex2 aldex2 \
		  --i-table $OUTDIR"/DA_"$X"_vs_"$Y"/filtered_table_"$level"_"$X"_"$Y".qza" \
		  --m-metadata-file $METADATA_FILE \
		  --m-metadata-column $METADATA_COLUMN \
		  --o-differentials $OUTDIR"/DA_"$X"_vs_"$Y"/differentials_"$level".qza"
		qiime aldex2 extract-differences \
		  --i-table $OUTDIR"/DA_"$X"_vs_"$Y"/differentials_"$level".qza" \
		  --o-differentials $OUTDIR"/DA_"$X"_vs_"$Y"/extracted_differences_"$level \
		  --p-sig-threshold 1 \
		  --p-effect-threshold 0 \
		  --p-difference-threshold 0
		#The tab separated file of differentially called features can be exported.
		qiime tools export \
		  --input-path $OUTDIR"/DA_"$X"_vs_"$Y"/extracted_differences_"$level".qza" \
		  --output-path $OUTDIR"/DA_"$X"_vs_"$Y
		mv $OUTDIR"/DA_"$X"_vs_"$Y"/differentials.tsv" $OUTDIR"/DA_"$X"_vs_"$Y"/"$level"_ALDEx2_table.tsv"
		rm $OUTDIR"/DA_"$X"_vs_"$Y"/filtered_table_"$level"_"$X"_"$Y".qza"
		fi
	done
  done
done
