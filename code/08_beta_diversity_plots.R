rm(list=ls())
getwd()
library("this.path")
library(tidyverse)
library(qiime2R)
library("ggpubr")

project_dir = this.path::here(..=1)
figures_dir=paste0(project_dir, "/results/figures")
supp_figures_dir=paste0(project_dir, "/results/supplementary_figures")

metadata_file = paste0(project_dir, "/data/covid19_study_metadata.tsv")
metadata<-read_q2metadata(metadata_file)
# metadata$disease_phase=factor(metadata$disease_phase)
in_dir=paste0(project_dir, "/data/intermediate/08_core_metrics/core_metrics_phylogenetic")
shannon<-read_qza(paste0(in_dir,"/shannon_vector.qza"))$data %>% rownames_to_column("SampleID") 
metadata = metadata %>% left_join(shannon) %>% filter(!is.na(shannon_entropy)) 
b_divers_measure=c("unweighted_unifrac",
                   "bray_curtis",
                   "jaccard",
                   "weighted_unifrac")
#figure 2
b="bray_curtis"
pcoa_results=read_qza(paste0(in_dir,"/",b,"_pcoa_results.qza"))
plot_data = pcoa_results$data$Vectors %>%
  select(`SampleID`, PC1, PC2) %>%
  left_join(metadata) %>%
  mutate(group_1 = factor(group_1,c("negative","asymptomatic","mild","severe"))) %>% #, labels = c("NE","AS","MI","SE")
  mutate(group_2 = factor(group_2,c("negative","asymptomatic recovery phase","asymptomatic early phase","mild","severe")))
plot_data %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(alpha=0.7, aes(shape=`group_2`, fill= `group_2`), color ="#000000", size=c(2,4,3,2,2)[factor(plot_data$group_2)])+ theme_q2r() +
  #scale_size_discrete(range = c(0.5,6)) +
  ggtitle("Bray-Curtis dissimilarity")+
  labs(x="PC1 (19.71 %)", y="PC2 (13.32 %)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "right")+
  stat_ellipse(aes(color=`group_1`),level = 0.9, show.legend=FALSE)+ 
  guides(color="none",
         fill=guide_legend(ncol=1),       # <-- two columns for fill legend
         shape=guide_legend(ncol=1))+
  scale_shape_manual(name="COVID-19 group", 
                     values=c(21,4,22, 21, 21),
                     #labels = c("negative", "asymptomatic recovery", "asymptomatic early", "mild", "severe")
                     ) + 
  scale_color_manual(values = c("#56B4E9", "#999999","orange", "brown"))+
  scale_fill_manual(name="COVID-19 group", 
                    values = c("#56B4E9", "#999999","#999999", "orange", "brown"),
                    #labels = c("negative", "asymptomatic recovery", "asymptomatic early", "mild", "severe")
                    )
ggsave(paste0(figures_dir,"/Figure_2.pdf"), height=4, width=4, device="pdf") # save a PDF 3 inches by 4 inches
ggsave(paste0(figures_dir,"/Figure_2.png"), height=4, width=4, device="png") # save a PDF 3 inches by 4 inches

#figure 2 revised
b="bray_curtis"
pcoa_results=read_qza(paste0(in_dir,"/",b,"_pcoa_results.qza"))
plot_data = pcoa_results$data$Vectors %>%
  select(`SampleID`, PC1, PC2) %>%
  left_join(metadata) %>%
  mutate(group_1 = factor(group_1,levels = c("asymptomatic","mild","severe", "negative"))) %>% 
  mutate(group_2 = factor(group_2,levels = c("asymptomatic recovery phase","asymptomatic early phase","mild","severe", "negative")))
plot_data %>%
  ggplot(aes(x=PC1, y=PC2)) +
  geom_point(alpha=0.7, aes(shape=`group_2`, fill= `group_2`), color ="#000000", size=c(4,3,2,2,2)[factor(plot_data$group_2)])+ theme_q2r() +
  #scale_size_discrete(range = c(0.5,6)) +
  ggtitle("Bray-Curtis dissimilarity")+
  labs(x="PC1 (19.71 %)", y="PC2 (13.32 %)") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom",
        legend.justification = "center",
        plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "in"))+ # top, right, bottom, left
  stat_ellipse(aes(color=`group_1`),level = 0.9, show.legend=FALSE)+ 
  guides(color="none",
         fill=guide_legend(ncol=2),       # <-- two columns for fill legend
         shape=guide_legend(ncol=3))+
  scale_shape_manual(name="COVID-19 group", 
                     values=c(4,22, 21, 21,21),
                     #labels = c("negative", "asymptomatic recovery", "asymptomatic early", "mild", "severe")
  ) + 
  scale_color_manual(values = c("#999999","orange", "brown", "#56B4E9"))+
  scale_fill_manual(name="COVID-19 group", 
                    values = c("#999999","#999999", "orange", "brown", "#56B4E9"),
                    #labels = c("negative", "asymptomatic recovery", "asymptomatic early", "mild", "severe")
  )
ggsave(paste0(figures_dir,"/Figure_2r.pdf"), height=4, width=4.1, device="pdf") # save a PDF 3 inches by 4 inches
ggsave(paste0(figures_dir,"/Figure_2r.jpg"), height=4, width=4.1, device="jpeg") # save a PDF 3 inches by 4 inches

#FIGURE S1
b="jaccard"
pcoa_results=read_qza(paste0(in_dir,"/",b,"_pcoa_results.qza"))
plot_s1_a = pcoa_results$data$Vectors %>%
  select(`SampleID`, PC1, PC2) %>%
  left_join(metadata) %>%
  mutate(group_1 = factor(group_1,c("negative","asymptomatic","mild","severe"))) %>%
  ggplot(aes(x=PC1, y=PC2, color=`group_1`, size=`shannon_entropy`)) +
  geom_point(alpha=0.5) + 
  theme_q2r() +
  #scale_shape_manual(values=c(16,1), name="shape") + 
  scale_size_continuous(name="Shannon index", range = c(0.5,6)) +
  scale_color_discrete(name="COVID-19 group", type=c("#56B4E9", "#999999","orange", "brown")) +
  ggtitle(str_to_title(gsub("_"," ",paste0(b, " distance"))))+
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_ellipse(level = 0.9, show.legend=FALSE) 
b="unweighted_unifrac"
pcoa_results=read_qza(paste0(in_dir,"/",b,"_pcoa_results.qza"))
plot_s1_b = pcoa_results$data$Vectors %>%
  select(`SampleID`, PC1, PC2) %>%
  left_join(metadata) %>%
  mutate(group_1 = factor(group_1,c("negative","asymptomatic","mild","severe"))) %>%
  ggplot(aes(x=PC1, y=PC2, color=`group_1`, size=`shannon_entropy`)) +
  geom_point(alpha=0.5) + 
  theme_q2r() +
  #scale_shape_manual(values=c(16,1), name="shape") + 
  scale_size_continuous(name="Shannon index", range = c(0.5,6)) +
  scale_color_discrete(name="COVID-19 group", type=c("#56B4E9", "#999999","orange", "brown")) +
  ggtitle(str_to_title(gsub("_"," ",paste0(b, " distance"))))+
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_ellipse(level = 0.9, show.legend=FALSE)  
b="weighted_unifrac"
pcoa_results=read_qza(paste0(in_dir,"/",b,"_pcoa_results.qza"))
plot_s1_c = pcoa_results$data$Vectors %>%
  select(`SampleID`, PC1, PC2) %>%
  left_join(metadata) %>%
  mutate(group_1 = factor(group_1,c("negative","asymptomatic","mild","severe"))) %>%
  ggplot(aes(x=PC1, y=PC2, color=`group_1`, size=`shannon_entropy`)) +
  geom_point(alpha=0.5) + 
  theme_q2r() +
  #scale_shape_manual(values=c(16,1), name="shape") + 
  scale_size_continuous(name="Shannon index", range = c(0.5,6)) +
  scale_color_discrete(name="COVID-19 group", type=c("#56B4E9", "#999999","orange", "brown")) +
  ggtitle(str_to_title(gsub("_"," ",paste0(b, " distance"))))+
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_ellipse(level = 0.9, show.legend=FALSE)  
figureS1 = ggarrange(plot_s1_a,plot_s1_b,plot_s1_c,
                     labels=c("A","B","c"),
                     ncol=1, nrow=3, common.legend = TRUE, legend = "right") 
figureS1
ggsave(paste0(supp_figures_dir,"/","Figure_S1.pdf"), height=5, width=3, device="pdf", scale = 1.8, bg = "white") # save a PDF 3 inches by 4 inches
ggsave(paste0(supp_figures_dir,"/","Figure_S1.png"), height=5, width=3, device="png", scale = 1.8, bg = "white") # save a PDF 3 inches by 4 inches



