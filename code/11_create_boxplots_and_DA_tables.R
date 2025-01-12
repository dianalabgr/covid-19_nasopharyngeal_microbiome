rm(list=ls())
library(qiime2R)
library(tidyverse)
library(xlsx)
library(tidyr)
library(ggplot2)
library(data.table)
library(wesanderson)

project_dir=this.path::here(..=1)
setwd(project_dir)
plots_dir=file.path(project_dir, "results/plots")
supp2_plots_dir=file.path(project_dir,"results/supplementary_plots_2")
supp_tables_dir=file.path(project_dir, "results/supplementary_tables")
dir.create(plots_dir, showWarnings = FALSE)
dir.create(supp2_plots_dir, showWarnings = FALSE)
dir.create(supp_tables_dir, showWarnings = FALSE)

metadata_file = paste0(project_dir, "/data/16S_3batches_metadata.tsv")
metadata<-read_q2metadata(metadata_file)
taxonomy_file="/mnt/raid1/philipos/corona_project/16S_analysis/data/reference/gg_13_8_otus/99_otus_taxonomy.qza"
in_dir="/mnt/raid1/philipos/corona_project/16S_analysis/data/intermediate/06_frequency_filtering"
DA_results_folder="/mnt/raid1/philipos/corona_project/16S_analysis/results/DA_analysis/"

comparisons = data.frame(
  group1=c("severe","mild","asymptomatic","severe","mild","severe"),
  group2=c("negative","negative","negative","asymptomatic","asymptomatic","mild"))

taxonomy_level="genus"
for(taxonomy_level in c("phylum","family","genus")){
  source(paste0(project_dir,"/code/10_Composition/arguments_",taxonomy_level,".R"))
  if (taxonomy_level=="genus"){
    table_name="filtered_table_l6.qza"
    differentials="l6_differentials.tsv"
    substitute1="k__.*g__"
    substitute2=";s__.*"
    substitute3="^k__.*f__"
    substitute4="f__"
  } else if (taxonomy_level=="family"){
    table_name="filtered_table_l5.qza"
    differentials="l5_differentials.tsv"
    substitute1="k__.*f__"
    substitute2=";g__.*"
    substitute3="^k__.*c__"
    substitute4="c__"
  } else if (taxonomy_level=="phylum"){
    table_name="filtered_table_l2.qza"
    differentials="l2_differentials.tsv"
    substitute1="k__.*p__"
    substitute2=";p__.*"
    substitute3="^k__.*f__"
    substitute4="f__"
  }
  significant_taxa_names=c()
  DA_taxa_direction=c()
  for (i in 1:nrow(comparisons)){
    group1=comparisons[i,"group1"]
    group2=comparisons[i,"group2"]
    differential_folder=paste0(DA_results_folder,"DA_",group1,"_vs_",group2,"/differentials")
   
      mytable=read_table(paste0(differential_folder,"/",differentials), col_names =T)[-1,] %>%
                                                   column_to_rownames("featureid") %>%
                                                   mutate(we.eBH=as.numeric(we.eBH)) %>%
                                                   filter(we.eBH<=0.05) %>%
                                                   transmute("in {group1} compared to {group2}" :=ifelse(`diff.btw`>0,"UP","DOWN")) %>%
                                                   rownames_to_column("Feature.ID")        
      
      significant_taxa_names=append(significant_taxa_names, unique(DA_taxa_direction$Feature.ID))
  
    if(is.null(DA_taxa_direction)){
      DA_taxa_direction=mytable
    }     else{
      DA_taxa_direction=merge(DA_taxa_direction,mytable, all=TRUE, by = 'Feature.ID') 
    }
  }
  
  #add Taxon column
  DA_taxa_direction = DA_taxa_direction %>% mutate(Taxon=gsub(substitute1, "", Feature.ID)) %>%
    mutate(Taxon=gsub(substitute2, "", Taxon)) %>%
    mutate(Taxon=gsub(substitute3,substitute4, Taxon)) %>%
    filter(Taxon!="")
    #select(-Feature.ID)
  #DA_taxa_direction %>% rename("Feature.ID"="Row.names")

  #
  OTUs<-read_qza(paste0(in_dir,"/",table_name))$data
  OTUs<-apply(OTUs, 2, function(x) x/sum(x)*100) #convert to percent
  OTUs= OTUs %>%
    as.data.frame() %>%
    rownames_to_column("Feature.ID")
  #we.eBH - Expected Benjamini-Hochberg corrected P value of Welch’s t test
  if(nrow(DA_taxa_direction)==0){next}
  
  #calculate median abundance per group
  OTUs_composition_table=
    OTUs %>%
    as.data.frame() %>%
    gather(-Feature.ID, key="SampleID", value="Abundance") %>%
    filter(Feature.ID %in% significant_taxa_names) %>%
    #mutate(Feature.ID=if_else((Feature.ID %in% significant_taxa_names),  Feature.ID, "_Low abundant")) %>% #flag features to be collapsed
    group_by(SampleID, Feature.ID) %>%
    #summarize(Abundance=sum(Abundance)) %>%   #for "Low abundant group" abundance calculation
    left_join(metadata[,1:25]) %>%
    #mutate(NormAbundance=log10(Abundance+0.01)) %>% # do a log10 transformation after adding a 0.01% pseudocount. Could also add 1 read before transformation to percent
    mutate(Taxon=gsub(substitute1, "", Feature.ID)) %>%
    mutate(Taxon=gsub(substitute2, "", Taxon)) %>%
    mutate(Taxon=gsub(substitute3,substitute4, Taxon))
  #save excel
  median_taxa_abundance=
    OTUs_composition_table %>%
    group_by(Taxon, type) %>%
    summarize(Abundance=median(Abundance)) %>%
    spread(key="type", value = "Abundance") %>%
    mutate_if(is.numeric, round, digits = 2) %>%
    relocate(negative, .before = asymptomatic)  %>%
    #relocate(Change, .before = negative) %>%
    # arrange(desc(Change)) %>%
    filter(Taxon!="")
  FINAL_TABLE=merge(DA_taxa_direction, median_taxa_abundance, all=TRUE, by="Taxon") %>%
    arrange(Taxon) %>%
    arrange(`in severe compared to negative`)
  write.xlsx(as.data.frame(FINAL_TABLE)%>%select(-Feature.ID), paste0(supp_tables_dir,"/Table_S1_",taxonomy_level,"_level_all_DA.xlsx"), row.names=FALSE, showNA=FALSE)
  #select only highly abundant taxa (where median abundance is at least 4% in at least one of the four patient groups:
  median_taxa_abundance= median_taxa_abundance %>% filter(max(severe,mild,asymptomatic,negative)>=4)
  FINAL_TABLE=merge(DA_taxa_direction, median_taxa_abundance, all.y=TRUE, by="Taxon") %>%
    arrange(Taxon) %>%
    arrange(`in severe compared to negative`)
  write.xlsx(as.data.frame(FINAL_TABLE)%>%select(-Feature.ID), paste0(supp_tables_dir,"/Table_S1_",taxonomy_level,"_level.xlsx"), row.names=FALSE, showNA=FALSE)
  top_taxa=unique(FINAL_TABLE$Taxon)
  
#plot some selected taxa
  selected_taxa=c("Actinobacteria","Firmicutes","Proteobacteria", "Fusobacteria","Bacteroidetes")
  if(taxonomy_level=="family"){selected_taxa=c(
    "Weeksellaceae","Alcaligenaceae","Bacillaceae","Corynebacteriaceae","Enterobacteriaceae","Fusobacteriaceae","Halomonadaceae", 
    "Pasteurellaceae","Prevotellaceae","Streptococcaceae", "Veillonellaceae")}
  if(taxonomy_level=="genus"){selected_taxa=c("Corynebacterium","Elizabethkingia","Fusobacterium", 
                                               "Haemophilus", "Halomonas", "Prevotella", "Serratia","Streptococcus", "Veillonella")}
  
  queries <- paste(selected_taxa, collapse = "|")
  OTUs_composition_table_selected = OTUs_composition_table %>% 
    filter(str_detect(Taxon,queries)) %>% 
    filter(Taxon %in% top_taxa) %>%
    mutate(type = factor(type, levels=c("negative","asymptomatic","COVID-19")))

  #create boxplot, all taxa
  myboxplot=qplot(Taxon, Abundance, data = OTUs_composition_table, color = factor(type, levels=c("negative","asymptomatic","COVID-19")), geom = "boxplot", fill=factor(type, levels=c("negative","asymptomatic","COVID-19")))+
    scale_fill_manual(values=c("#56B4E9","#999999","#E69F00")) +
    scale_color_manual(values=c("#56B4E9","#999999","#E69F00"))+
    theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
    labs(fill = "type", color="type")+
    ylab("Relative Abundance (%)")
  myboxplot
  #ggsave(paste0(supp2_plots_dir,"/",taxonomy_level,"_boxplot.png"),myboxplot, "png", width=26.5, height=15.2, units="cm")
  
  
  #create boxplot, selected taxa only
  myboxplot=qplot(Taxon, Abundance, data = OTUs_composition_table_selected, color = factor(type2, levels=c("1.negative","2.asymptomatic","3.mild","4.severe" )), geom = "boxplot", fill=factor(type2, levels=c("1.negative","2.asymptomatic","3.mild","4.severe")))+
    scale_fill_manual(values=c("white","white","white","white")) +
    scale_color_manual(values=c("#56B4E9","#999999","orange","brown"))+
    theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
    labs(fill = "Patient group", color="Patient group")+
    ylab("Relative Abundance (%)")+
    xlab(str_to_sentence(taxonomy_level))
  myboxplot
  ggsave(paste0(supp2_plots_dir,"/selected_",taxonomy_level,"_boxplot.png"),myboxplot, "png", width=26.5, height=15.2, units="cm")
  if(taxonomy_level=="phylum"){p1=myboxplot}
  if(taxonomy_level=="family"){p2=myboxplot}
  if(taxonomy_level=="genus"){p3=myboxplot}
  
}
figure3 = ggpubr::ggarrange(p1,p2,p3,
                     labels=c("A","B"),
                     ncol=1, nrow=3)
ggsave(paste0(plots_dir, "/Figure3.pdf"), width=26.5, height=35, units="cm", device="pdf") # save a PDF 3 inches by 4 inches
ggsave(paste0(plots_dir, "/Figure3.png"), width=26.5, height=35, units="cm", device="png") # save a PNG 3 inches by 4 inches

