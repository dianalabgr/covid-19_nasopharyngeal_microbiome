rm(list=ls())
library(qiime2R)
library(tidyverse)
library(xlsx)
library(tidyr)
library(ggplot2)
library(data.table)
library(wesanderson)
library(stringr)
library(dplyr)

project_dir=this.path::here(..=1)
# setwd(project_dir)
figures_dir=file.path(project_dir, "results/figures")
supp_tables_dir=file.path(project_dir, "results/supplementary_material")
other_figures_dir = file.path(project_dir, "results/intermediate_figures")
dir.create(figures_dir, showWarnings = FALSE)
dir.create(supp_tables_dir, showWarnings = FALSE)
dir.create(other_figures_dir, showWarnings = FALSE)

metadata_file = paste0(project_dir, "/data/covid19_study_metadata.tsv")
metadata<-read_q2metadata(metadata_file) %>% filter(group_2!="asymptomatic recovery phase") %>% filter(!SampleID %in% c("SE-04", "SE-15","NE-27","SE-08"))

taxonomy_file="/mnt/raid1/philipos/corona_project/16S_analysis/data/reference/gg_13_8_otus/99_otus_taxonomy.qza"
in_dir="/mnt/raid1/philipos/corona_project/16S_analysis/data/intermediate/06_frequency_filtering"
DA_results_folder="/mnt/raid1/philipos/corona_project/16S_analysis/results/Differential_Abundance_Analysis/"

comparisons = data.frame(
  group1=c("severe","mild","asymptomatic","severe","mild","severe"),
  group2=c("negative","negative","negative","asymptomatic","asymptomatic","mild"))

taxonomy_level="family"
taxonomy_level="phylum"
taxonomy_level="genus"
for(taxonomy_level in c("phylum","family","genus")){
  if (taxonomy_level=="genus"){
    table_name="filtered_table_l6.qza"
    differentials="l6_ALDEx2_table.tsv"
    substitute1="k__.*g__"
    substitute2=";s__.*"
    substitute3="^k__.*f__"
    substitute4="f__"
  } else if (taxonomy_level=="family"){
    table_name="filtered_table_l5.qza"
    differentials="l5_ALDEx2_table.tsv"
    substitute1="k__.*f__"
    substitute2=";g__.*"
    substitute3="^k__.*c__"
    substitute4="c__"
  } else if (taxonomy_level=="phylum"){
    table_name="filtered_table_l2.qza"
    differentials="l2_ALDEx2_table.tsv"
    substitute1="k__.*p__"
    substitute2=";p__.*"
    substitute3="^k__.*f__"
    substitute4="f__"
  }
  DA_taxa_direction=c()
  for (i in 1:nrow(comparisons)){
    group1=comparisons[i,"group1"]
    group2=comparisons[i,"group2"]
    differential_folder=paste0(DA_results_folder,"DA_",group1,"_vs_",group2)
    mytable=read_table(paste0(differential_folder,"/",differentials), col_names =T)[-1,] %>%
                                                 column_to_rownames("featureid") %>%
                                                 mutate(we.eBH=as.numeric(we.eBH)) %>%  #  #we.eBH - Expected Benjamini-Hochberg corrected P value of Welch’s t test
                                                 filter(we.eBH<=0.05) %>%
                                                 transmute("in {group1} compared to {group2}" :=ifelse(`diff.btw`>0, 
                                                                                                       paste0("UP, p=",sprintf('%.3f',we.eBH)),
                                                                                                       paste0("DOWN, p=",sprintf('%.3f',we.eBH)))) %>%
                                                 rownames_to_column("Feature.ID")       
    if(is.null(DA_taxa_direction)){
      DA_taxa_direction=mytable
    }     else{
      DA_taxa_direction=merge(DA_taxa_direction,mytable, all=TRUE, by = 'Feature.ID') 
    }
  }
  significant_taxa=unique(DA_taxa_direction$Feature.ID)
  #add Taxon column
  DA_taxa_direction = DA_taxa_direction %>% mutate(Taxon=gsub(substitute1, "", Feature.ID)) %>%
    mutate(Taxon=gsub(substitute2, "", Taxon)) %>%
    mutate(Taxon=gsub(substitute3,substitute4, Taxon)) %>%
    filter(Taxon!="")
  OTUs<-read_qza(paste0(in_dir,"/",table_name))$data
  OTUs<-apply(OTUs, 2, function(x) x/sum(x)*100) #convert to percent
  OTUs= OTUs %>%
    as.data.frame() %>%
    rownames_to_column("Feature.ID")
  #we.eBH - Expected Benjamini-Hochberg corrected P value of Welch’s t test
  #if no significant differences detected then skip
  if(nrow(DA_taxa_direction)==0){next}
  #calculate median abundance per group_original
  OTUs_composition_table=
    OTUs %>%
    as.data.frame() %>%
    gather(-Feature.ID, key="SampleID", value="Abundance") %>%
    filter(Feature.ID %in% significant_taxa) %>%
    group_by(SampleID, Feature.ID) %>%
    inner_join(metadata[,1:25]) %>%
    mutate(Taxon=gsub(substitute1, "", Feature.ID)) %>%
    mutate(Taxon=gsub(substitute2, "", Taxon)) %>%
    mutate(Taxon=gsub(substitute3,substitute4, Taxon))
  OTUs_composition_table$group_1=factor(OTUs_composition_table$group_1, c("negative","asymptomatic","mild","severe"))
  #TABLE S1
  #save excel
  median_taxa_abundance=
    OTUs_composition_table %>%
    group_by(Taxon, group_1) %>%
    summarize(Abundance=median(Abundance)) %>%
    spread(key="group_1", value = "Abundance") %>%
    mutate_if(is.numeric, round, digits = 2) %>%
    relocate(negative, .before = asymptomatic)  %>%
    filter(Taxon!="")
  FINAL_TABLE=merge(DA_taxa_direction, median_taxa_abundance, all=TRUE, by="Taxon") %>%
    arrange(Taxon) %>%
    arrange(`in severe compared to negative`)
  write.xlsx(as.data.frame(FINAL_TABLE)%>%select(-Feature.ID), paste0(supp_tables_dir,"/Supplementary_Table_S1_",taxonomy_level,"_level_all_DA_r.xlsx"), row.names=FALSE, showNA=FALSE)
  #select only highly abundant taxa (where median abundance is at least 4% in at least one of the four patient group:
  median_taxa_abundance= median_taxa_abundance %>% filter(max(severe,mild,asymptomatic,negative)>=4)
  FINAL_TABLE=merge(DA_taxa_direction, median_taxa_abundance, all.y=TRUE, by="Taxon") %>%
    arrange(Taxon) %>%
    arrange(`in severe compared to negative`)
  write.xlsx(as.data.frame(FINAL_TABLE)%>%select(-Feature.ID), paste0(supp_tables_dir,"/Supplementary_Table_S1_",taxonomy_level,"_level_r.xlsx"), row.names=FALSE, showNA=FALSE)
  top_significant_taxa=unique(FINAL_TABLE$Taxon)
  
  
  #figures
  increasing_taxa = NULL
  decreasing_taxa = NULL
  # incr_decr_taxa = NULL
  
 
  incr_severe <- grepl("^UP", FINAL_TABLE$`in severe compared to negative`)
  incr_mild   <- grepl("^UP", FINAL_TABLE$`in mild compared to negative`)
  increasing_taxa <- FINAL_TABLE$Taxon[incr_severe | incr_mild]
  decr_severe <- grepl("^DOWN", FINAL_TABLE$`in severe compared to negative`)
  decr_mild   <- grepl("^DOWN", FINAL_TABLE$`in mild compared to negative`)
  decreasing_taxa <- FINAL_TABLE$Taxon[decr_severe | decr_mild]
  
 
  if(!is.null(increasing_taxa)){
    
    # DA_taxa_direction contains: Taxon, and columns 'in X compared to Y' with "UP, p=..." or "DOWN, p=..."
    
    # Pivot to long format for processing
    
    pvalue_long <- FINAL_TABLE %>% filter(Taxon %in% increasing_taxa) %>% arrange(Taxon) %>%
      pivot_longer(cols = starts_with("in "), names_to = "comparison", values_to = "result") %>%
      filter(!is.na(result), str_detect(result, "p=")) %>%
      mutate(
        p_value = as.numeric(str_extract(result, "(?<=p=)[0-9\\.]+")),
        group1 = word(comparison, 2),   # e.g., "severe"
        group2 = word(comparison, 5)    # e.g., "negative"
      ) %>%
      filter(p_value <= 0.05) %>%
      select(Taxon, group1, group2, p_value)
    # # 1. Calculate the Dodge Offsets
    # Get group order and set manual dodge mapping
    group_order <- c("negative", "asymptomatic", "mild", "severe")
    ndodge <- length(group_order)
    dodge_width <- 0.8
    offset <- (-dodge_width/2) + ((2 * (1:ndodge) - 1) * dodge_width / (2 * ndodge))
    names(offset) <- group_order
    
    # Now add geom_vline to your plot_incr/plot_decr
    
    # 2. Convert Taxon to Numeric x-Axis Positions
    taxon_levels <- FINAL_TABLE %>%
      filter(Taxon %in% increasing_taxa) %>%
      pull(Taxon) %>%      # extract as character vector
      unique() %>% 
      sort %>%        # get unique taxon names
      as.list()  
    # convert vector to list
    # taxon_levels <- unique(FINAL_TABLE$Taxon)
    taxon_map <- setNames(seq_along(taxon_levels), taxon_levels)
    guide_positions <- as.vector(sapply(taxon_map, function(x) x + offset))
    
    # 3. Prepare Annotation Data
    
    pvalue_long <- pvalue_long %>%
      mutate(
        base_x     = as.numeric(factor(Taxon, levels = taxon_levels)),
        xstart     = base_x + offset[group1],
        xend       = base_x + offset[group2],
        xmid       = (xstart + xend) / 2   # for centering the p-value label
      )
    
    # 
    y_positions <- OTUs_composition_table %>%
      group_by(Taxon) %>%
      summarize(y = max(Abundance, na.rm = TRUE)) # or use a whisker value to exclude outliers if preferred
    
    # Join y values to pvalue_long
    pvalue_long <- pvalue_long %>%
      left_join(y_positions, by = "Taxon") %>%
      mutate(
        y_line = y + 5,                # line height above boxplot
        y_text = y_line + 2          # adjust space as needed for text
      )
    
    
    pvalue_long <- pvalue_long %>%
      group_by(Taxon) %>%
      mutate(row_index = row_number() - 1,   # 0-based index: 0 for first, 1 for second, etc.
             y      = y + row_index*5,
             y_line = y_line + row_index*5,
             y_text = y_text + row_index*5
      ) %>%
      ungroup()

    # Prepare: match taxa to numeric positions
    taxon_levels <- FINAL_TABLE %>%
      filter(Taxon %in% increasing_taxa) %>%
      pull(Taxon) %>%
      unique() %>%
      sort()
    
    # Add numeric x position to your plotting data
    OTUs_composition_table_i <- OTUs_composition_table %>%
      filter(Taxon %in% increasing_taxa) %>%
      mutate(xpos = as.numeric(factor(Taxon, levels = taxon_levels)))
    
    #
    plot_incr = ggplot(data = OTUs_composition_table_i, aes(x = xpos, y = Abundance)) +
      geom_vline(xintercept = 1:length(increasing_taxa), color = "gray80", size = 0.3) +
      geom_vline(xintercept = guide_positions, linetype = "solid", color = "white", size = 0.3) +
      geom_boxplot(width = 0.8, aes(fill = group_1, group = interaction(group_1, Taxon)), outlier.shape = NA) +
      geom_point(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8, seed = 1234), 
                 aes(group = group_1, fill = group_1), pch = 21, alpha = 0.7) +
      scale_fill_manual(values = c("#56B4E9", "#999999", "orange", "brown")) +
      # Here: set continuous x scale with labels at taxon positions
      scale_x_continuous(
        breaks = seq_along(taxon_levels),
        labels = taxon_levels
      ) +
      labs(fill = "COVID-19 group") +
      theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.text = element_text(size = 12)) +
      ylab("Relative Abundance (%)")+
      xlab(str_to_sentence(taxonomy_level))
    
  
    
    plot_incr = plot_incr+
      # Lines above the relevant group boxes
      geom_segment(
        data = pvalue_long,
        aes(x = xstart, xend = xend, y = y_line, yend = y_line, group = interaction(Taxon, group1, group2)),
        inherit.aes = FALSE,
        color = "black"
      ) +
      geom_text(
        data = pvalue_long,
        aes(x = xmid, y = y_text, label = sprintf("%.2f", p_value), group = interaction(Taxon, group1, group2)),
        inherit.aes = FALSE,
        color = "black",
        size = 2  # or as desired
      ) 
   
      ggsave(paste0(other_figures_dir, "/Figure4_",taxonomy_level,"_increasing_taxa_r2.pdf"), width=26.5, height=12, units="cm", device="pdf") 
      ggsave(paste0(other_figures_dir, "/Figure4_",taxonomy_level,"_increasing_taxa_r2.png"), width=26.5, height=12, units="cm", device="png") 
  }
  if(!is.null(decreasing_taxa)){
   
    
    # DA_taxa_direction contains: Taxon, and columns 'in X compared to Y' with "UP, p=..." or "DOWN, p=..."
    
    # Pivot to long format for processing
    
    pvalue_long <- FINAL_TABLE %>% filter(Taxon %in% decreasing_taxa) %>% arrange(Taxon) %>%
      pivot_longer(cols = starts_with("in "), names_to = "comparison", values_to = "result") %>%
      filter(!is.na(result), str_detect(result, "p=")) %>%
      mutate(
        p_value = as.numeric(str_extract(result, "(?<=p=)[0-9\\.]+")),
        group1 = word(comparison, 2),   # e.g., "severe"
        group2 = word(comparison, 5)    # e.g., "negative"
      ) %>%
      filter(p_value <= 0.05) %>%
      select(Taxon, group1, group2, p_value)
    # # 1. Calculate the Dodge Offsets
    # Get group order and set manual dodge mapping
    group_order <- c("negative", "asymptomatic", "mild", "severe")
    ndodge <- length(group_order)
    dodge_width <- 0.8
    offset <- (-dodge_width/2) + ((2 * (1:ndodge) - 1) * dodge_width / (2 * ndodge))
    names(offset) <- group_order
    # 2. Convert Taxon to Numeric x-Axis Positions
    taxon_levels <- FINAL_TABLE %>%
      filter(Taxon %in% decreasing_taxa) %>%
      pull(Taxon) %>%      # extract as character vector
      unique() %>% 
      sort %>%        # get unique taxon names
      as.list()  
    # convert vector to list
    # taxon_levels <- unique(FINAL_TABLE$Taxon)
    taxon_map <- setNames(seq_along(taxon_levels), taxon_levels)
    guide_positions <- as.vector(sapply(taxon_map, function(x) x + offset))
    
    # 3. Prepare Annotation Data
    
    pvalue_long <- pvalue_long %>%
      mutate(
        base_x     = as.numeric(factor(Taxon, levels = taxon_levels)),
        xstart     = base_x + offset[group1],
        xend       = base_x + offset[group2],
        xmid       = (xstart + xend) / 2   # for centering the p-value label
      )
    
    # 
    y_positions <- OTUs_composition_table %>%
      group_by(Taxon) %>%
      summarize(y = max(Abundance, na.rm = TRUE)) # or use a whisker value to exclude outliers if preferred
    
    # Join y values to pvalue_long
    pvalue_long <- pvalue_long %>%
      left_join(y_positions, by = "Taxon") %>%
      mutate(
        y_line = y + 5,                # line height above boxplot
        y_text = y_line + 2          # adjust space as needed for text
      )
    
    
    pvalue_long <- pvalue_long %>%
      group_by(Taxon) %>%
      mutate(row_index = row_number() - 1,   # 0-based index: 0 for first, 1 for second, etc.
             y      = y + row_index*5,
             y_line = y_line + row_index*5,
             y_text = y_text + row_index*5
      ) %>%
      ungroup()
    
    # Prepare: match taxa to numeric positions
    taxon_levels <- FINAL_TABLE %>%
      filter(Taxon %in% decreasing_taxa) %>%
      pull(Taxon) %>%
      unique() %>%
      sort()
    
    # Add numeric x position to your plotting data
    OTUs_composition_table_d <- OTUs_composition_table %>%
      filter(Taxon %in% decreasing_taxa) %>%
      mutate(xpos = as.numeric(factor(Taxon, levels = taxon_levels)))
    
    #
    plot_decr = ggplot(data = OTUs_composition_table_d, aes(x = xpos, y = Abundance)) +
      geom_vline(xintercept = 1:length(decreasing_taxa), color = "gray80", size = 0.3) +
      geom_vline(xintercept = guide_positions, linetype = "solid", color = "white", size = 0.3) +
      geom_boxplot(width = 0.8, aes(fill = group_1, group = interaction(group_1, Taxon)), outlier.shape = NA) +
      geom_point(position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.8, seed = 1234), 
                 aes(group = group_1, fill = group_1), pch = 21, alpha = 0.7) +
      scale_fill_manual(values = c("#56B4E9", "#999999", "orange", "brown")) +
      # Here: set continuous x scale with labels at taxon positions
      scale_x_continuous(
        breaks = seq_along(taxon_levels),
        labels = taxon_levels
      ) +
      labs(fill = "COVID-19 group") +
      theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.text = element_text(size = 12)) +
      ylab("Relative Abundance (%)")+
      xlab(str_to_sentence(taxonomy_level))
    
    
    plot_decr = plot_decr+
      # Lines above the relevant group boxes
      geom_segment(
        data = pvalue_long,
        aes(x = xstart, xend = xend, y = y_line, yend = y_line, group = interaction(Taxon, group1, group2)),
        inherit.aes = FALSE,
        color = "black"
      ) +
      geom_text(
        data = pvalue_long,
        aes(x = xmid, y = y_text, label = sprintf("%.2f", p_value), group = interaction(Taxon, group1, group2)),
        inherit.aes = FALSE,
        color = "black",
        size = 2  # or as desired
      )
    
    ggsave(paste0(other_figures_dir, "/Figure4_",taxonomy_level,"_decreasing_taxa_r2.pdf"), width=26.5, height=12, units="cm", device="pdf") 
    ggsave(paste0(other_figures_dir, "/Figure4_",taxonomy_level,"_decreasing_taxa_r2.png"), width=26.5, height=12, units="cm", device="png") 
  }
 
  
  if (taxonomy_level=="genus"){
    figure_4e = plot_incr
    figure_4f = plot_decr
  } else if (taxonomy_level=="family"){
    figure_4c = plot_incr
    figure_4d = plot_decr
  } else if (taxonomy_level=="phylum"){
    figure_4a = plot_incr
    figure_4b = plot_decr
  }
  
}

figure4 = ggpubr::ggarrange(figure_4a,figure_4b,figure_4c,figure_4d,figure_4e,figure_4f,
                    labels=c("A","B","C","D","E","F"),
                    ncol=2, nrow=3,common.legend = T)
ggsave(paste0(figures_dir, "/Figure_4r.pdf"), width=26.5, height=35, units="cm", device="pdf") 
ggsave(paste0(figures_dir, "/Figure_4r.jpg"), width=26.5, height=35, units="cm", device="jpg") 

