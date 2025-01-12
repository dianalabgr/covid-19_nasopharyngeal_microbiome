rm(list=ls())
library(qiime2R)
library(this.path)
library(dplyr)
library(tibble)

project_dir=this.path::here(..=2)
setwd(project_dir)
out_dir=file.path(project_dir, "results/abundances/")
dir.create(out_dir, showWarnings = FALSE)


feature_table=paste0(project_dir,"/data/intermediate/03_qiime_otu_clustering/table-cr-99.qza")
feature_table=paste0(project_dir,"/data/intermediate/06_frequency_filtering/filtered_table_l6.qza")


                     #06_frequency_filtering")

OTUs<-read_qza(feature_table)$data
OTUs<-apply(OTUs, 2, function(x) x/sum(x)*100) 
OTUs<-round(OTUs,2)

#convert to percent
OTUs= OTUs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID")


write.table(OTUs, file = paste0(out_dir,"Metagenomics_relative_abundances_l6.tsv"), sep = "\t", row.names = FALSE, quote=FALSE)
