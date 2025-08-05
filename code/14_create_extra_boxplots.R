rm(list=ls())
library(qiime2R)
library(tidyverse)
library(xlsx)
library(tidyr)
library(ggplot2)
library(data.table)
library(wesanderson)

project_dir=this.path::here(..=1)
# setwd(project_dir)
figures_dir=file.path(project_dir, "results/figures")
supp_tables_dir=file.path(project_dir, "results/supplementary_material")
other_figures_dir = file.path(project_dir, "results/supplementary_figures")
dir.create(figures_dir, showWarnings = FALSE)
dir.create(supp_tables_dir, showWarnings = FALSE)
dir.create(other_figures_dir, showWarnings = FALSE)

metadata_file = paste0(project_dir, "/data/covid19_study_metadata.tsv")
metadata<-read_q2metadata(metadata_file) %>% filter(group_2!="asymptomatic recovery phase") %>% filter(!SampleID %in% c("SE-04", "SE-15","NE-27","SE-08"))

taxonomy_file="/mnt/raid1/philipos/corona_project/16S_analysis/data/reference/gg_13_8_otus/99_otus_taxonomy.qza"
in_dir="/mnt/raid1/philipos/corona_project/16S_analysis/data/intermediate/06_frequency_filtering"
rf_results_folder="/mnt/raid1/philipos/corona_project/16S_analysis/data/intermediate/10_randomforest/classifier_1_l6"
important_taxa = read_table(paste0(rf_results_folder,"/feature_importance_table.tsv"), col_names = TRUE)[2:20,]

taxonomy_level="genus"
table_name="filtered_table_l6.qza"


process_taxonomy <- function(taxonomy) {
  # Split taxonomy into levels
  levels <- unlist(strsplit(taxonomy, ";"))
  
  # Find last assigned index
  last_assigned_index <- -1
  for (i in seq_along(levels)) {
    # Split prefix and name
    parts <- strsplit(levels[i], "__")[[1]]
    if (length(parts) == 2) {
      name <- trimws(parts[2])
      # Check if assigned (non-empty name)
      if (nchar(name) > 0) {
        last_assigned_index <- i
      }
    }
  }
  if (last_assigned_index == -1) {
    # No assigned level found; return original
    return(taxonomy)
  }
  # Extract the last assigned level without prefix
  last_level <- levels[last_assigned_index]
  last_name <- sub("^[a-z]__","", last_level) # Remove prefix
    # Extract remaining unassigned levels including their prefixes
  if (last_assigned_index < length(levels)) {
    unassigned_levels <- levels[(last_assigned_index + 1):length(levels)]
    # Join unassigned levels keeping prefixes, separated by ;
    unassigned_str <- paste(unassigned_levels, collapse = ";")
    # Combine last assigned name + ";" + unassigned parts (if any)
    combined <- if (nchar(unassigned_str) > 0) {
      paste0(last_name, ";", unassigned_str)
    } else {
      last_name
    }
    return(combined)
  } else {
    # No unassigned suffix levels
    return(last_name)
  }
}

OTUs<-read_qza(paste0(in_dir,"/",table_name))$data
OTUs<-apply(OTUs, 2, function(x) x/sum(x)*100) #convert to percent
OTUs= OTUs %>% as.data.frame() %>%
  rownames_to_column("Feature.ID")  %>% 
  filter(Feature.ID %in% important_taxa$id)

OTUs_composition_table=
  OTUs %>%
  as.data.frame() %>%
  gather(-Feature.ID, key="SampleID", value="Abundance") %>%
  filter(Feature.ID %in% important_taxa$id) %>%
  group_by(SampleID, Feature.ID) %>%
  inner_join(metadata[,1:25]) %>%
  mutate(Taxon=Feature.ID)
OTUs_composition_table$group_1=factor(OTUs_composition_table$group_1, c("negative","asymptomatic","mild","severe"))
OTUs_composition_table$Taxon <- sapply(OTUs_composition_table$Taxon, process_taxonomy)

 
plot_rf_taxa = ggplot(data=OTUs_composition_table,
                   aes(x=Taxon, y= Abundance)) +
  geom_boxplot(width = 0.8,aes(fill=group_1), outlier.shape=NA) +
  geom_point(position = position_jitterdodge(jitter.width=0.4, dodge.width = 0.8, seed = 1234), 
             aes(group=group_1, fill=group_1), pch=21, alpha= 0.7) +
  scale_fill_manual(values=c("#56B4E9","#999999","orange","brown")) +
  labs(fill="COVID-19 group") +
  theme(axis.text.x = element_text(angle = 55, hjust = 1),
        legend.text = element_text(size=12),
        legend.position = "top")  +
  guides(fill=guide_legend(ncol=4))+
  scale_y_sqrt(  breaks = seq(0, 100, 10)) +
  ylab("Relative Abundance (%)")


ggsave(paste0(other_figures_dir, "/Supplementary_Figure_S3.pdf"), width=26.5, height=12, units="cm", device="pdf") 
ggsave(paste0(other_figures_dir, "/Supplementary_Figure_S3.jpg"), width=26.5, height=12, units="cm", device="jpg") 




