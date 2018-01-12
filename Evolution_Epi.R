setwd("~/Box Sync/Epipedobates project/BIOINFORMATICS/ANALYSIS/epi_rana_radseq")

library(tidyverse)
library(cowplot)

rana_2brad_snpdist <- read_csv("rana2brad_snpdist.csv")
rana_2brad_sumstats <- read_csv("rana2brad_sumstats.csv")
epi_2brad_snpdist <- read_csv("epi2brad_snpdist.csv")
epi_2brad_sumstats <- read_csv("epi2brad_sumstats.csv")

# Easier labeling of figures
ranaclustlabs = c("80", "81", "84", "85", "86", "87", "88", "89", "90", "91", "93", "95")
epiclustlabs = c("80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91", "92", "93", "94", "95")
ranalabs = c("Rber_a","Rber_b","Rber","Rbla_4","Rbla_5","Rchi_a","Rchi_b","Rchi","Rneo_0","Rneo_7","Rsph_0","Rsph_4")
epilabs = c("Ahah_a","Ahah_b","Ahah","Eant","Eant_a","Eant_b","Ebou_3","Ebou_6","Etri_6","Etri_2","Snub_8","Snub_9")

### 2bRAD_SNPDIST_SUM OF VARIABLE SITES ###
p.rana_2brad_snpdist <-
  rana_2brad_snpdist %>% 
  group_by(clust_threshold) %>%
  filter(sum_var == max(sum_var)) %>%
ggplot(aes(x=clust_threshold, y=sum_var)) + 
  geom_col() +
  scale_x_discrete(labels=ranaclustlabs) +
  scale_y_continuous(labels=scales::comma, limits=c(0,400000)) +
  ggtitle("Rana_2brad")

p.epi_2brad_snpdist <-
  epi_2brad_snpdist %>% 
  group_by(clust_threshold) %>%
  filter(sum_var == max(sum_var)) %>%
  ggplot(aes(x=clust_threshold, y=sum_var)) + 
  geom_col() +
  scale_x_discrete(labels=epiclustlabs) +
  scale_y_continuous(labels=scales::comma, limits=c(0,400000)) +
  ggtitle("Epi_2brad")
  

plot_grid(p.epi_2brad_snpdist, p.rana_2brad_snpdist)

### 2bRAD_SNPDIST_SUM OF PI SITES ###

p.rana_2brad_snpdist_pis <-
  rana_2brad_snpdist %>% 
  group_by(clust_threshold) %>%
  filter(sum_pis == max(sum_pis)) %>%
  ggplot(aes(x=clust_threshold, y=sum_pis)) + 
  geom_col() +
  scale_x_discrete(labels=ranaclustlabs) +
  scale_y_continuous(labels=scales::comma, limits=c(0,1200000)) +
  ggtitle("Rana_2brad")

p.epi_2brad_snpdist_pis <-
  epi_2brad_snpdist %>% 
  group_by(clust_threshold) %>%
  filter(sum_pis == max(sum_pis)) %>%
  ggplot(aes(x=clust_threshold, y=sum_pis)) + 
  geom_col() +
  scale_x_discrete(labels=epiclustlabs) +
  scale_y_continuous(labels=scales::comma, limits=c(0,1200000)) +
  ggtitle("Epi_2brad")

plot_grid(p.epi_2brad_snpdist_pis, p.rana_2brad_snpdist_pis)

### 2bRAD_SUMSTATS HEATMAPS ###
p.rana_2brad_sumstats <-
  rana_2brad_sumstats %>%
  ggplot(aes(clust_threshold, loci_assembly, group=sample, color=sample)) +
  geom_line() +
  scale_x_discrete(labels=ranaclustlabs) +
  scale_y_continuous(labels=scales::comma, limits=c(0,200000)) +
  ggtitle("Rana_2bRAD")
  
p.epi_2brad_sumstats <-
  epi_2brad_sumstats %>%
  ggplot(aes(clust_threshold, loci_assembly, group=sample, color=sample)) +
  geom_line() +
  scale_x_discrete(labels=epiclustlabs) +
  scale_y_continuous(labels=scales::comma, limits=c(0,200000)) +
  ggtitle("Epi_2bRAD")

plot_grid(p.epi_2brad_sumstats, p.rana_2brad_sumstats)

#### COMBINED PLOT
  ggplot(data=rana_2brad_sumstats, aes(clust_threshold, loci_assembly, group=sample, color=sample)) +
  geom_line() +
  geom_point(shape=19) +
  scale_y_continuous(labels=scales::comma, limits=c(0,200000)) +
  geom_line(data=epi_2brad_sumstats, aes(clust_threshold, loci_assembly, group=sample, color=sample)) +
  scale_x_discrete(labels=ranaclustlabs)
  #  geom_point(shape=20)  

p.epi_2brad_sumstats <-
  epi_2brad_sumstats %>%
  ggplot(aes(clust_threshold, loci_assembly, group=sample, color=sample)) +
  geom_line() +
  scale_x_discrete(labels=epiclustlabs) +
  scale_y_continuous(labels=scales::comma, limits=c(0,200000)) +
  ggtitle("Epi_2bRAD")
