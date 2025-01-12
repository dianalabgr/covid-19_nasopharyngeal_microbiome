rm(list=ls())
library("this.path")
library(tidyverse)
library(qiime2R)
# install.packages("ggpubr")
library("ggpubr")
library("vegan")

setwd(dirname2(this.path()))
setwd("../")
project_dir=getwd()
getwd()
plots_dir=paste0(project_dir, "/results/plots")
metadata_file = paste0(project_dir, "/data/16S_3batches_metadata.tsv")
metadata<-read_q2metadata(metadata_file)

in_dir=paste0(project_dir, "/data/intermediate/08_core_metrics/core_metrics_phylogenetic")
shannon<-read_qza(paste0(in_dir,"/shannon_vector.qza"))$data %>% rownames_to_column("SampleID") 

metadata = metadata %>% left_join(shannon) %>% filter(!is.na(shannon_entropy)) 
metadata$type = factor(metadata$type, c("negative","asymptomatic","mild","severe"))

#figure 1b (4 groups)

pvalues=read_csv(paste0(in_dir,"/shannon_4_groups_kruskal-wallis-pairwise-test.csv"))
pvalues$`Group 1` =gsub(" (.*)","",pvalues$`Group 1` )
pvalues$`Group 2` =gsub(" (.*)","",pvalues$`Group 2` )
pvalues=pvalues %>% rename("group1"= "Group 1", "group2"= "Group 2") 
pvalues=pvalues %>% mutate(across(where(is.numeric), ~ round(., digits = 2)))

p2 = ggboxplot(metadata, x = "type", y = "shannon_entropy",
               color = "type", palette = c("#56B4E9", "#999999","orange", "brown"),
               add = "jitter")+
  stat_pvalue_manual(pvalues, label = "q-value",  tip.length = 0.03,
                     y.position = c(8,8.5, 8.5,  9, 9, 9.5), bracket.shorten = 0         ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    ylab("Shannon index") + 
  xlab("Patient group") +labs(colour="patient group")


#figure 1a (2 groups)
pvalues=read_csv(paste0(in_dir,"/shannon_2_groups_kruskal-wallis-pairwise-test.csv"))
pvalues$`Group 1` =gsub(" (.*)","",pvalues$`Group 1` )
pvalues$`Group 2` =gsub(" (.*)","",pvalues$`Group 2` )
pvalues=pvalues %>% rename("group1"= "Group 1", "group2"= "Group 2") 
pvalues=pvalues %>% mutate(across(where(is.numeric), ~ round(., digits = 2)))

p1 = ggboxplot(metadata, x = "test", y = "shannon_entropy",
          color = "test", palette = c("#56B4E9", "#a86d32"),
          add = "jitter")+
  stat_pvalue_manual(pvalues, label = "q-value",  tip.length = 0.03,
                     y.position = c(9), bracket.shorten = 0         ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ylab("Shannon index")+ 
  xlab("COVID-19 test") +labs(colour="COVID-19 test")

figure01 = ggarrange(p1,p2,
                     labels=c("A","B"),
                     ncol=1, nrow=2)
ggsave(paste0(plots_dir, "/Figure1.pdf"), height=8, width=5, device="pdf") # save a PDF 3 inches by 4 inches
ggsave(paste0(plots_dir, "/Figure1.png"), height=8, width=5, device="png") # save a PNG 3 inches by 4 inches
