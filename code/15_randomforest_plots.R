# Load necessary libraries
library(readr)
library(pROC)    # for ROC curves
library(ggplot2) # for plotting
library(this.path())
library(qiime2R)
library(dplyr)
library(tidyr)
library(ggrepel)


# Read probability and true label tables
dir = "/mnt/raid1/philipos/corona_project/16S_analysis/data/intermediate/10_randomforest/classifier_1_l6"
project_dir=this.path::here(..=1)
prob <- read_tsv(paste0(dir,"/exported_probabilities/class_probabilities.tsv"))
pred <- read_tsv('exported_predictions/predictions.tsv')
# Suppose you have a metadata TSV with true class labels:
metadata_file=paste0(project_dir,"/data/covid19_study_metadata.tsv")

meta<-read_q2metadata(metadata_file)
meta$group_1 = factor(meta$group_1,c("negative","asymptomatic","mild","severe"))  #, labels = c("NE","AS","MI","SE")
meta$group_2 = gsub(" ","\n",meta$group_2)
meta$group_2 = factor(meta$group_2,c('negative',"asymptomatic\nrecovery\nphase","asymptomatic\nearly\nphase",'mild','severe'), levels=c('negative',"asymptomatic\nrecovery\nphase","asymptomatic\nearly\nphase",'mild','severe'))
meta = metadata %>% mutate(group_2 = factor(group_2, labels = c('negative',"asympt.\nrecovery\nphase","asympt.\nearly\nphase",'mild','severe')))

# Merge by SampleID if necessary
colnames(prob)[1] = "SampleID"
df <- merge(prob, meta, by.x="SampleID", by.y="SampleID")
# Suppose true_label column is named "Class", and positive class is "disease"


# Suppose your probs table is 'prob' (as shown), predictions are in 'pred'
all_classes <- c("asymptomatic", "mild", "negative", "severe")

roc_list <- list()
auc_list <- c()
cls = "severe"
for (cls in all_classes) {
  # 1-vs-rest
  true_cls <- as.integer(pred$prediction == cls)  # 1 if predicted class is 'cls'
  prob_cls <- prob[[cls]]
  roc_cls <- roc(true_cls, prob_cls, quiet=TRUE)
  auc_list[cls] <- auc(roc_cls)
  roc_df <- data.frame(
    fpr = 1 - roc_cls$specificities,
    tpr = roc_cls$sensitivities,
    threshold = roc_cls$thresholds,
    class = cls
  )
  roc_list[[cls]] <- roc_df
}

roc_data <- bind_rows(roc_list, .id = "class")

# Plot all classes' ROC curves & annotate a few threshold values for each
ggplot(roc_data, aes(x=fpr, y=tpr, color=class)) +
  geom_line(size=1) +
  geom_abline(slope=1, intercept=0, linetype='dashed', color='grey') +
  # Annotate threshold cutoffs at select points
  geom_text_repel(
    data = roc_data %>% group_by(class) %>%
      slice(round(seq(1, n(), length.out=5))),
    aes(label=round(threshold,2)),
    size=3, segment.color='grey'
  ) +
  labs(
    title = "ROC Curves for Multiclass Random Forest Classification",
    x = "1 - Specificity", y = "Sensitivity"
  ) +
  scale_color_brewer(palette = "Dark2") +
  theme_minimal()
