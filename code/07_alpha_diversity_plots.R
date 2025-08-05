rm(list=ls())
library("this.path")
library(tidyverse)
library(qiime2R)
# install.packages("ggpubr")
library("ggpubr")
library("vegan")

project_dir = this.path::here(..=1)
figures_dir=paste0(project_dir, "/results/figures")
metadata_file = paste0(project_dir, "/data/covid19_study_metadata.tsv")
metadata<-read_q2metadata(metadata_file)

in_dir=paste0(project_dir, "/data/intermediate/08_core_metrics/core_metrics_phylogenetic")
shannon<-read_qza(paste0(in_dir,"/shannon_vector.qza"))$data %>% rownames_to_column("SampleID") 

metadata = metadata %>% left_join(shannon) %>% filter(!is.na(shannon_entropy)) 
metadata$group_1 = factor(metadata$group_1, c("negative","asymptomatic","mild","severe"))

#figure 1b (4 groups)
pvalues=read_csv(paste0(in_dir,"/shannon_4_groups_kruskal-wallis-pairwise-test.csv"))
pvalues$`Group 1` =gsub(" (.*)","",pvalues$`Group 1` )
pvalues$`Group 2` =gsub(" (.*)","",pvalues$`Group 2` )
pvalues=pvalues %>% rename("group1"= "Group 1", "group2"= "Group 2") 
pvalues=pvalues %>% mutate(across(where(is.numeric), ~ round(., digits = 2)))
pvalues$`q-value`<- sprintf("%.2f", pvalues$`q-value`)

p2 = ggboxplot(metadata, x = "group_1", y = "shannon_entropy",
               color = "group_1", palette = c("#56B4E9", "#999999","orange", "brown"),
               add = "jitter")+
  stat_pvalue_manual(pvalues, label = "q-value",  size = 3.5, tip.length = 0.01,
                     y.position = c(8.5, 8, 9.6, 9, 8, 10.1), 
                     bracket.shorten = 0 ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    ylab("Shannon index") + 
  xlab("COVID-19 group") +labs(colour="COVID-19 group")

#figure 1a (2 groups)
pvalues=read_csv(paste0(in_dir,"/shannon_2_groups_kruskal-wallis-pairwise-test.csv"))
pvalues$`Group 1` =gsub(" (.*)","",pvalues$`Group 1` )
pvalues$`Group 2` =gsub(" (.*)","",pvalues$`Group 2` )
pvalues=pvalues %>% rename("group1"= "Group 1", "group2"= "Group 2") 
pvalues=pvalues %>% mutate(across(where(is.numeric), ~ round(., digits = 2)))
pvalues$`q-value`<- sprintf("%.2f", pvalues$`q-value`)

p1 = ggboxplot(metadata, x = "test", y = "shannon_entropy",
          color = "test", palette = c("#56B4E9", "#a86d32"),
          add = "jitter")+
  stat_pvalue_manual(pvalues, label = "q-value", size = 3.5, tip.length = 0.01,
                     y.position = c(9), bracket.shorten = 0         ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ylab("Shannon index")+ 
  xlab("COVID-19 test") +labs(colour="COVID-19 test")

figure01 = ggarrange(p1,p2,
                     labels=c("A","B"),
                     ncol=1, nrow=2)
ggsave(paste0(figures_dir, "/Figure_1.pdf"), height=8, width=5.5, device="pdf") # save a PDF 3 inches by 4 inches
ggsave(paste0(figures_dir, "/Figure_1.png"), height=8, width=5.5, device="png") # save a PNG 3 inches by 4 inches
