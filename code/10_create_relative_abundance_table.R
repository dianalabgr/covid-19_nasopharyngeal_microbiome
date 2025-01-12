rm(list=ls())
library(qiime2R)
library(this.path)
library(dplyr)
library(tibble)
library(tidyverse)  


project_dir=this.path::here(..=1)
setwd(project_dir)
out_dir=file.path(project_dir, "results/abundance_tables/")
dir.create(out_dir, showWarnings = FALSE)


#SPECIES LEVEL
feature_table=paste0(project_dir,"/data/intermediate/06_frequency_filtering/filtered_table_l6.qza")
OTUs<-read_qza(feature_table)$data
OTUs<-apply(OTUs, 2, function(x) x/sum(x)*100) 
OTUs<-round(OTUs,2)
#convert to percent
OTUs= OTUs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID")
# #check if they sum to 100%
# otusmatrix = data.matrix(OTUs, rownames.force = NA)
# colSums(otusmatrix)
write.table(OTUs, file = paste0(out_dir,"Relative_percent_abundances_species.tsv"), sep = "\t", row.names = FALSE, quote=FALSE)

#FAMILY LEVEL
feature_table=paste0(project_dir,"/data/intermediate/06_frequency_filtering/filtered_table_l5.qza")
OTUs<-read_qza(feature_table)$data
OTUs<-apply(OTUs, 2, function(x) x/sum(x)*100) 
OTUs<-round(OTUs,2)
#convert to percent
OTUs= OTUs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID")
# #check if they sum to 100%
# otusmatrix = data.matrix(OTUs, rownames.force = NA)
# colSums(otusmatrix)
write.table(OTUs, file = paste0(out_dir,"Relative_percent_abundances_family.tsv"), sep = "\t", row.names = FALSE, quote=FALSE)

#PHYLUM LEVEL
feature_table=paste0(project_dir,"/data/intermediate/06_frequency_filtering/filtered_table_l2.qza")
OTUs<-read_qza(feature_table)$data
OTUs<-apply(OTUs, 2, function(x) x/sum(x)*100) 
OTUs<-round(OTUs,2)
#convert to percent
OTUs= OTUs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID")
# #check if they sum to 100%
# otusmatrix = data.matrix(OTUs, rownames.force = NA)
# colSums(otusmatrix)
write.table(OTUs, file = paste0(out_dir,"Relative_percent_abundances_phylum.tsv"), sep = "\t", row.names = FALSE, quote=FALSE)


#ALL TAXONOMY LEVELS
feature_table=paste0(project_dir,"/data/intermediate/06_frequency_filtering/filtered_table.qza")
OTUs<-read_qza(feature_table)$data
OTUs<-apply(OTUs, 2, function(x) x/sum(x)*100) 
OTUs<-round(OTUs,2)
#convert to percent
OTUs= OTUs %>%
  as.data.frame() %>%
  rownames_to_column("Feature.ID")
# #check if they sum to 100%
# otusmatrix = data.matrix(OTUs, rownames.force = NA)
# colSums(otusmatrix)
write.table(OTUs, file = paste0(out_dir,"Relative_percent_abundances_all_levels.tsv"), sep = "\t", row.names = FALSE, quote=FALSE)
