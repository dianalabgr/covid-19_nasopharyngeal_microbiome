rm(list=ls())
getwd()
library("this.path")
library(tidyverse)
library(qiime2R)
library("ggpubr")

setwd(dirname2(this.path()))
setwd("../")
project_dir=getwd()
plots_dir=paste0(project_dir, "/results/plots")
supp_plots_dir=paste0(project_dir, "/results/supplementary_plots")
supp2_plots_dir = paste0(project_dir, "/results/supplementary_plots_2")

metadata_file = paste0(project_dir, "/data/16S_3batches_metadata.tsv")
metadata<-read_q2metadata(metadata_file)
metadata$disease_phase=factor(metadata$disease_phase)



in_dir=paste0(project_dir, "/data/intermediate/08_core_metrics/core_metrics_phylogenetic")

shannon<-read_qza(paste0(in_dir,"/shannon_vector.qza"))$data %>% rownames_to_column("SampleID") 
metadata = metadata %>% left_join(shannon) %>% filter(!is.na(shannon_entropy)) 
b_divers_measure=c("unweighted_unifrac",
                   "bray_curtis",
                   "jaccard",
                   "weighted_unifrac")

#SUPPLEMENTARY FIGURES S2
#point size is Shannon index
for(b in b_divers_measure){
  pcoa_results=read_qza(paste0(in_dir,"/",b,"_pcoa_results.qza"))
  pcoa_results$data$Vectors %>%
    select(`SampleID`, PC1, PC2) %>%
    left_join(metadata) %>%
    mutate(type = factor(type,c("negative","asymptomatic","mild","severe"))) %>%
    ggplot(aes(x=PC1, y=PC2, color=`type`, size=`shannon_entropy`)) +
    geom_point(alpha=0.5) + 
    theme_q2r() +
    #scale_shape_manual(values=c(16,1), name="shape") + 
    scale_size_continuous(name="Shannon index", range = c(0.5,6)) +
    scale_color_discrete(name="Patient group", type=c("#56B4E9", "#999999","orange", "brown")) +
    ggtitle(str_to_title(gsub("_"," ",paste0(b, " distance"))))+
    theme(plot.title = element_text(hjust = 0.5)) +
    stat_ellipse(level = 0.9, show.legend=FALSE)   
  ggsave(paste0(supp2_plots_dir,"/", b,"_PC_plot.pdf"), height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches
  ggsave(paste0(supp2_plots_dir,"/", b,"_PC_plot.png"), height=4, width=5, device="png") # save a PDF 3 inches by 4 inches
}

#figure 2
b="bray_curtis"
pcoa_results=read_qza(paste0(in_dir,"/",b,"_pcoa_results.qza"))
pcoa_results$data$Vectors %>%
  select(`SampleID`, PC1, PC2) %>%
  left_join(metadata) %>%
  mutate(type = factor(type,c("negative","asymptomatic","mild","severe"))) %>%
  ggplot(aes(x=PC1, y=PC2, color=`type`, size=`shannon_entropy`)) +
  geom_point(alpha=0.5) + 
  theme_q2r() +
  #scale_shape_manual(values=c(16,1), name="shape") + 
  scale_size_continuous(name="Shannon index", range = c(0.5,6)) +
  scale_color_discrete(name="Patient group", type=c("#56B4E9", "#999999","orange", "brown")) +
  ggtitle("Bray-Curtis dissimilarity")+
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_ellipse(level = 0.9, show.legend=FALSE)   
ggsave(paste0(plots_dir,"/Figure2.pdf"), height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches
ggsave(paste0(plots_dir,"/Figure2.png"), height=4, width=5, device="png") # save a PDF 3 inches by 4 inches

#figure S1
b="bray_curtis"
pcoa_results=read_qza(paste0(in_dir,"/",b,"_pcoa_results.qza"))
plot_data = pcoa_results$data$Vectors %>%
  select(`SampleID`, PC1, PC2) %>%
  left_join(metadata) %>%
  mutate(type = factor(type,c("negative","asymptomatic","mild","severe"))) 
plot_data %>%
  ggplot(aes(x=PC1, y=PC2, color=`type`)) +
  geom_point(alpha=0.8, aes(shape=`disease_phase`), size=c(1,3,3)[factor(plot_data$disease_phase)])+ theme_q2r() +
  scale_shape_manual(name="Asymptomatic disease phase" , values=c(20,1,4)) + 
  #scale_size_discrete(range = c(0.5,6)) +
  scale_color_discrete(name="Patient group", type=c("#56B4E9", "#999999","orange", "brown"))+
  ggtitle("Bray-Curtis dissimilarity")+
  theme(plot.title = element_text(hjust = 0.5))+
  stat_ellipse(level = 0.9, show.legend=FALSE)
ggsave(paste0(supp_plots_dir,"/FigureS1n.pdf"), height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches
ggsave(paste0(supp_plots_dir,"/FigureS1n.png"), height=4, width=5, device="png") # save a PDF 3 inches by 4 inches

#FIGURE S2
b="jaccard"
pcoa_results=read_qza(paste0(in_dir,"/",b,"_pcoa_results.qza"))
plot_s2_a = pcoa_results$data$Vectors %>%
  select(`SampleID`, PC1, PC2) %>%
  left_join(metadata) %>%
  mutate(type = factor(type,c("negative","asymptomatic","mild","severe"))) %>%
  ggplot(aes(x=PC1, y=PC2, color=`type`, size=`shannon_entropy`)) +
  geom_point(alpha=0.5) + 
  theme_q2r() +
  #scale_shape_manual(values=c(16,1), name="shape") + 
  scale_size_continuous(name="Shannon index", range = c(0.5,6)) +
  scale_color_discrete(name="Patient group", type=c("#56B4E9", "#999999","orange", "brown")) +
  ggtitle(str_to_title(gsub("_"," ",paste0(b, " distance"))))+
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_ellipse(level = 0.9, show.legend=FALSE) 
b="unweighted_unifrac"
pcoa_results=read_qza(paste0(in_dir,"/",b,"_pcoa_results.qza"))
plot_s2_b = pcoa_results$data$Vectors %>%
  select(`SampleID`, PC1, PC2) %>%
  left_join(metadata) %>%
  mutate(type = factor(type,c("negative","asymptomatic","mild","severe"))) %>%
  ggplot(aes(x=PC1, y=PC2, color=`type`, size=`shannon_entropy`)) +
  geom_point(alpha=0.5) + 
  theme_q2r() +
  #scale_shape_manual(values=c(16,1), name="shape") + 
  scale_size_continuous(name="Shannon index", range = c(0.5,6)) +
  scale_color_discrete(name="Patient group", type=c("#56B4E9", "#999999","orange", "brown")) +
  ggtitle(str_to_title(gsub("_"," ",paste0(b, " distance"))))+
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_ellipse(level = 0.9, show.legend=FALSE)  
b="weighted_unifrac"
pcoa_results=read_qza(paste0(in_dir,"/",b,"_pcoa_results.qza"))
plot_s2_c = pcoa_results$data$Vectors %>%
  select(`SampleID`, PC1, PC2) %>%
  left_join(metadata) %>%
  mutate(type = factor(type,c("negative","asymptomatic","mild","severe"))) %>%
  ggplot(aes(x=PC1, y=PC2, color=`type`, size=`shannon_entropy`)) +
  geom_point(alpha=0.5) + 
  theme_q2r() +
  #scale_shape_manual(values=c(16,1), name="shape") + 
  scale_size_continuous(name="Shannon index", range = c(0.5,6)) +
  scale_color_discrete(name="Patient group", type=c("#56B4E9", "#999999","orange", "brown")) +
  ggtitle(str_to_title(gsub("_"," ",paste0(b, " distance"))))+
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_ellipse(level = 0.9, show.legend=FALSE)  
figureS2 = ggarrange(plot_s2_a,plot_s2_b,plot_s2_c,
                     labels=c("A","B","c"),
                     ncol=1, nrow=3, common.legend = TRUE, legend = "right") 
figureS2
ggsave(paste0(supp_plots_dir,"/","Figure_S2.pdf"), height=5, width=3, device="pdf", scale = 1.8, bg = "white") # save a PDF 3 inches by 4 inches
ggsave(paste0(supp_plots_dir,"/","Figure_S2.png"), height=5, width=3, device="png", scale = 1.8, bg = "white") # save a PDF 3 inches by 4 inches

#the point shapes are depicting the asymptomatics disease phase for ALL b metrics
for(b in b_divers_measure){
  pcoa_results=read_qza(paste0(in_dir,"/",b,"_pcoa_results.qza"))
  plot_data = pcoa_results$data$Vectors %>%
    select(`SampleID`, PC1, PC2) %>%
    left_join(metadata) %>%
    mutate(type = factor(type,c("negative","asymptomatic","mild","severe"))) 
  plot_data %>%
    ggplot(aes(x=PC1, y=PC2, color=`type`)) +
    geom_point(alpha=0.8, aes(shape=`disease_phase`), size=c(1,3,3)[factor(plot_data$disease_phase)])+ theme_q2r() +
    scale_shape_manual(name="Asymptomatic disease phase" , values=c(20,1,4)) + 
    #scale_size_discrete(range = c(0.5,6)) +
    scale_color_discrete(name="Patient group", type=c("#56B4E9", "#999999","orange", "brown"))+
    ggtitle(str_to_title(gsub("_"," ",b)))+
    theme(plot.title = element_text(hjust = 0.5))+
    stat_ellipse(level = 0.9, show.legend=FALSE)
  ggsave(paste0(supp2_plots_dir,"/", b,"_PC_plot_AS_disease_phase.pdf"), height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches
  ggsave(paste0(supp2_plots_dir,"/", b,"_PC_plot_AS_disease_phase.png"), height=4, width=5, device="png") # save a PDF 3 inches by 4 inches
}

#POSITIVE VS NEGATIVE
for(b in b_divers_measure){
  pcoa_results=read_qza(paste0(in_dir, "/",b,"_pcoa_results.qza"))
  pcoa_results$data$Vectors %>%
    select(`SampleID`, PC1, PC2) %>%
    left_join(metadata) %>%
    ggplot(aes(x=PC1, y=PC2, color=`test`, size=`shannon_entropy`)) +
    geom_point(alpha=0.5) + 
    theme_q2r() +
    scale_size_continuous(name="Shannon Diversity", range = c(1,4)) +
    scale_color_discrete(name="SARS-CoV-2 Test", type=c("#56B4E9", "#a89051"))+
    ggtitle(gsub("_"," ",b))+
    theme(plot.title = element_text(hjust = 0.5))+
    stat_ellipse(level = 0.9, show.legend=FALSE)
  ggsave(paste0(supp2_plots_dir, "/", b,"_PC_plot_pos_vs_neg.pdf"), height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches
  ggsave(paste0(supp2_plots_dir, "/", b,"_PC_plot_pos_vs_neg.png"), height=4, width=5, device="png") # save a PDF 3 inches by 4 inches
}
# 


#plot Bray Curtis and Age
b="bray_curtis"
metadata = metadata %>% mutate(age_group = case_when(age<18~"child",age<=35 ~"young", age<=55 ~"middle age", age>55 ~"older"))
metadata$age_group = factor(metadata$age_group, levels=c("child","young","middle age","older"))
pcoa_results=read_qza(paste0(in_dir,"/",b,"_pcoa_results.qza"))
pcoa_results$data$Vectors %>%
  select(`SampleID`, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=age_group)) +
  geom_point(alpha=0.5) + 
  theme_q2r() +
  scale_size_continuous(name="Shannon index", range = c(0.5,6)) +
  scale_color_discrete(name="Age group")+ 
  ggtitle("Bray-Curtis dissimilarity")+
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_ellipse(level = 0.9, show.legend=FALSE)   
ggsave(paste0(supp2_plots_dir,"/Figure_bray_curtis_Age.pdf"), height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches
ggsave(paste0(supp2_plots_dir,"/Figure_bray_curtis_Age.png"), height=4, width=5, device="png") # save a PDF 3 inches by 4 inches

#plot Bray Curtis and Sex
b="bray_curtis"
pcoa_results=read_qza(paste0(in_dir,"/",b,"_pcoa_results.qza"))
pcoa_results$data$Vectors %>%
  select(`SampleID`, PC1, PC2) %>%
  left_join(metadata) %>%
  ggplot(aes(x=PC1, y=PC2, color=sex)) +
  geom_point(alpha=0.5) + 
  theme_q2r() +
  scale_size_continuous(name="Shannon index", range = c(0.5,6)) +
  scale_color_discrete(name="Sex")+ 
  ggtitle("Bray-Curtis dissimilarity")+
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_ellipse(level = 0.9, show.legend=FALSE)   
ggsave(paste0(supp2_plots_dir,"/Figure_bray_curtis_Sexnew.pdf"), height=4, width=5, device="pdf") # save a PDF 3 inches by 4 inches
ggsave(paste0(supp2_plots_dir,"/Figure_bray_curtis_Sexnew.png"), height=4, width=5, device="png") # save a PDF 3 inches by 4 inches

