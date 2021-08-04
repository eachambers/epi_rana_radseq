setwd("~/Box Sync/Epipedobates project/SuppMaterials/5_Data_visualization")

library(tidyverse)
library(cowplot)
library(ggplot2)
library(scales)

theme_set(theme_cowplot())

## This code does the following:
##     Builds figures to determine appropriate clustering threshold value 
##     as it relates to consensus reads and loci in assembly (Fig. S1).

##    FILES REQUIRED:
##          2_Bioinformatics/iPyrad/clustthreshold_ddrad_sumstats.csv

# Data upload and tidying ------------------------------------------------

# Upload files; sumstats files are generated using iPyrad output
sumstats <- read_csv("../2_Bioinformatics/iPyrad/clustthreshold_ddrad_sumstats.csv")

# Remove unneccessary columns
sumstats <-
  sumstats %>%
  select(-X1, -`Unnamed: 0`)

# Add columns with taxon and cluster threshold percent specified
sumstats <-
  sumstats %>%
  mutate(clust_threshold = if_else(grepl("80",sumstats$file_name), "80",
                                 if_else(grepl("81",sumstats$file_name), "81",
                                         if_else(grepl("82",sumstats$file_name), "82",
                                                 if_else(grepl("83",sumstats$file_name), "83",
                                                         if_else(grepl("84",sumstats$file_name), "84",
                                                                 if_else(grepl("85",sumstats$file_name), "85",
                                                                         if_else(grepl("86",sumstats$file_name), "86",
                                                                                 if_else(grepl("87",sumstats$file_name), "87",
                                                                                         if_else(grepl("88",sumstats$file_name), "88",
                                                                                                 if_else(grepl("89",sumstats$file_name), "89",
                                                                                                         if_else(grepl("90",sumstats$file_name), "90",
                                                                                                                 if_else(grepl("91",sumstats$file_name), "91",
                                                                                                                         if_else(grepl("92",sumstats$file_name), "92",
                                                                                                                                 if_else(grepl("93",sumstats$file_name), "93",
                                                                                                                                         if_else(grepl("94",sumstats$file_name), "94",
                                                                                                                                                 if_else(grepl("95",sumstats$file_name), "95",
                                                                                                                                                         if_else(grepl("rana2brad_outfiles",sumstats$file_name), "85", NA_character_)))))))))))))))))) %>%
  mutate(taxon = if_else(grepl("epi", sumstats$file_name), "epi",
                           if_else(grepl("rana", sumstats$file_name), "rana", NA_character_)))

# Custom color palette (partially based off wesanderson package)
palette <- c("#8D529E", "#BA85A2", "#6B2391", "#CA3518", "#C66E1F", 
             "#BA8517", "#D9C339", "#086955", "#7C935E", "#80A7B3", "#386784", "#008CA2")

# Ensure that things are ordered correctly
sumstats$sample <- factor(sumstats$sample, levels=c("Rber_T1113a", "Rber_T1113b", "Rber_T1114", "Rneo_T480", "Rneo_T527", "Rsph_T25870", "Rsph_T26064",
                                                              "Rbla_D2864", "Rbla_D2865", "Rchi_T2034a", "Rchi_T2034b", "Rchi_T2049", "Eant_T6859a", "Eant_T6859b", "Eant_T6857", "Etri_T6836", "Etri_T6842", "Ebou_R0153", "Ebou_R0156",
                                                              "Snub_R0158", "Snub_R0159", "Ahah_R0089a", "Ahah_R0089b", "Ahah_R0090"))


# BUILDING FIGURES -----------------------------------------------


# Consensus reads ---------------------------------------------------------

## To combine these all into a single plot (with loci_assembly), make each figure separately

# Epipedobates
p.epiconsens <-
  sumstats %>%
  filter(taxon=="epi") %>%
  ggplot(aes(as.numeric(clust_threshold), reads_consens, group=sample, color=sample)) +
  geom_line(size=1) +
  geom_point(size=3) +
  geom_vline(xintercept=91, linetype="dashed", size=1) +
  scale_y_continuous(labels=scales::comma, 
                     limits=c(0,200000), 
                     breaks=seq(0,200000, by=50000)) +
  theme(legend.position="none",
        axis.text = element_text(size=14),
        axis.title = element_text(size=17),
        panel.grid.major.x = element_line(size=.75, color="snow2"),
        panel.grid.minor.x = element_line(size=.75, color="snow2")) +
  scale_x_continuous(minor_breaks = seq(80,95,1)) +
  scale_color_manual(values=palette) +
  xlab("Clustering threshold") +
  ylab("Consensus reads")
  

# Rana
p.ranaconsens <-
  sumstats %>%
  filter(taxon=="rana") %>%
  ggplot(aes(as.numeric(clust_threshold), reads_consens, group=sample, color=sample)) +
  geom_line(size=1) +
  geom_point(size=3) +
  geom_vline(xintercept = 91, linetype="dashed", size=1) +
  scale_y_continuous(labels=scales::comma,
                     limits=c(0,200000),
                     breaks=seq(0,200000, by=50000)) +
  theme(legend.position="none",
        axis.text = element_text(size=14),
        axis.title = element_text(size=17),
        panel.grid.major.x = element_line(size=.75, color="snow2"),
        panel.grid.minor.x = element_line(size=.75, color="snow2")) +
  scale_x_continuous(minor_breaks = seq(80,95,1)) +
  scale_color_manual(values=palette) +
  xlab("Clustering threshold") +
  ylab("Consensus reads")

# Loci in assembly --------------------------------------------------------

# Epipedobates
p.epiloci <-
  sumstats %>%
  filter(taxon=="epi") %>%
  ggplot(aes(as.numeric(clust_threshold), loci_assembly, group=sample, color=sample)) +
  geom_line(size=1) +
  geom_point(size=3) +
  geom_vline(xintercept = 91, linetype="dashed", size=1) +
  scale_color_manual(values=palette) +
  scale_y_continuous(labels=scales::comma,
                     limits=c(0,30000),
                     breaks=seq(0,30000, by=5000)) +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=17),
        panel.grid.major.x = element_line(size=.75, color="snow2"),
        panel.grid.major.x = element_line(size=.75, color="snow2")) +
  scale_x_continuous(minor_breaks = seq(80,95,1)) +
  xlab("Clustering threshold") +
  ylab("Loci in assembly")

# Rana
p.ranaloci <-
  sumstats %>%
  filter(taxon=="rana") %>%
  ggplot(aes(as.numeric(clust_threshold), loci_assembly, group=sample, color=sample)) +
  geom_line(size=1) +
  geom_point(size=3) +
  geom_vline(xintercept = 91, linetype="dashed", size=1) +
  scale_color_manual(values=palette) +
  scale_y_continuous(labels=scales::comma,
                     limits=c(0,65000),
                     breaks=seq(0,65000, by=10000)) +
  theme(axis.text = element_text(size=14),
        axis.title = element_text(size=17),
        panel.grid.major.x = element_line(size=.75, color="snow2"),
        panel.grid.major.x = element_line(size=.75, color="snow2")) +
  scale_x_continuous(minor_breaks = seq(80,95,1)) +
  xlab("Clustering threshold") +
  ylab("Loci in assembly")


# Complete the plot figure -------------------------------------------------------

# Pull the legends so that you only have one per taxon in finished plot
legendepi <- get_legend(p.epiloci)
legendrana <- get_legend(p.ranaloci)

# Now plot, adding in legend afterward
p.epi <- plot_grid(p.epiconsens, p.epiloci + theme(legend.position="none"))
p.epi <- plot_grid(p.epi, legendepi, rel_widths = c(1,.15))


p.rana <- plot_grid(p.ranaconsens, p.ranaloci + theme(legend.position="none"))
p.rana <- plot_grid(p.rana, legendrana, rel_widths = c(1,.15))

# Build overall figure
p.clust_threshold <- plot_grid(p.rana, p.epi, nrow=2)

# Combine everything together and export as PDF
save_plot("FigS1_clustthreshold.pdf", p.clust_threshold, base_aspect_ratio = 1.6, base_height = 11)