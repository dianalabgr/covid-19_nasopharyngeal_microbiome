conda activate qiime2_2021.02
source tab-qiime
#conda info
GREENGENES_FOLDER="/mnt/raid1/philipos/corona_project/16S_analysis/data/reference/gg_13_8_otus"

#Import greengenes reference sequences
qiime tools import \
--type 'FeatureData[Sequence]' \
--input-path $GREENGENES_FOLDER/rep_set/99_otus.fasta \
--output-path $GREENGENES_FOLDER/99_otus.qza

qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path $GREENGENES_FOLDER/taxonomy/99_otu_taxonomy.txt \
--output-path $GREENGENES_FOLDER/99_otus_taxonomy.qza

#convert tree to .qza
cd $GREENGENES_FOLDER/trees
python 00a_convert_tree_to_qza.py 99_otus.tree
mv 99_otus.qza ../ 99_otus_tree.qza

#convert it to .tree
cd $GREENGENES_FOLDER
qiime tools export \
     --input-path 99_otus_tree.qza \
     --output-path ./

sed 's\; \|\g' tree.nwk > 99_otus_tree.nwk
	 
	 