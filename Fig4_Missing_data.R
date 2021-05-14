setwd("~/Box Sync/Epipedobates project/SuppMaterials/5_Data_visualization")

library(devtools)
library(tidyverse)
library(ggplot2)
library(cowplot)

theme_set(theme_cowplot())

## The following code generates Fig. 4, which illustrates the amount of missing data
## for each of the sampling depths and for each dataset (only SNP datasets).
## Amounts of missing data were calculated using the PAUP missdata function.

##    FILES REQUIRED:
##          data_files_input_into_scripts/missing_data_snps.txt

# Import data -------------------------------------------------------------

miss <- read_tsv("./data_files_input_into_scripts/missing_data_snps.txt")

# Remove replicate samples
miss <-
  miss %>% 
  filter(sample_ID != "Rber_T1113b", sample_ID != "Rchi_T2034b",
         sample_ID != "Eant_T6859b", sample_ID != "Ahah_R0089b")

# Calculate average percent missing data for each of the sampling depths
# and append new column with that value
new <-
miss %>%
  group_by(.dots=c("dataset", "sampling_depth")) %>%
  summarize(avg_missing = mean(percent_missing))

new_rana <-
new %>%
  filter(dataset=="rana2brad" | dataset=="ranaddrad")
  
new_epi <-
  new %>%
  filter(dataset=="epi2brad" | dataset=="epiddrad")

# Build figures -----------------------------------------------------------

palette <- c("#BA85A2", "#6B2391", "#CA3518", "#C66E1F", 
             "#BA8517", "#D9C339", "#086955", "#7C935E", "#80A7B3", "#386784")

# Missing data - per sample -----------------------------------------------

# Ensure that things are ordered correctly (in phylo order) by changing levels
miss$sample_ID <- factor(miss$sample_ID, levels=c("Rber_T1113a", "Rber_T1114", 
                                                  "Rneo_T480", "Rneo_T527", "Rsph_T25870", 
                                                  "Rsph_T26064", "Rbla_D2864", "Rbla_D2865", 
                                                  "Rchi_T2034a", "Rchi_T2049", 
                                                  "Eant_T6859a", "Eant_T6857", 
                                                  "Etri_T6836", "Etri_T6842", "Ebou_R0153", 
                                                  "Ebou_R0156", "Snub_R0158", "Snub_R0159", 
                                                  "Ahah_R0089a", "Ahah_R0090"))


# Build plots -------------------------------------------------------------

# Epipedobates
p.epiind <-
miss %>%
  filter(dataset=="epi2brad" | dataset=="epiddrad") %>%
  ggplot(aes(x=sampling_depth, y=percent_missing, group=sample_ID, color=sample_ID)) +
  geom_line(size=2) +
  geom_point(size=7) +
  ylab("Missing data (%)") +
  xlab("Sampling depth") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=25),
        axis.line.x = element_line(size=2),
        axis.line.y = element_line(size=2),
        axis.ticks.x = element_line(size=2),
        axis.ticks.y = element_line(size=2),
        strip.background = element_blank()) +
  ggtitle("Epipedobates") +
  facet_grid(~dataset) +
  scale_color_manual(values=palette) +
  scale_y_continuous(limits=c(0,100))

## Add in avg amounts of missing data on top
p.epi_MD <-
p.epiind +
  geom_line(data=new_epi, aes(y=avg_missing, x=sampling_depth, group=dataset), color="black", size=2) +
  geom_point(data=new_epi, aes(y=avg_missing, x=sampling_depth, group=dataset), color="black", size=7)


# ===== do the same for Rana now =====
p.ranaind <-
miss %>%
  filter(dataset=="rana2brad" | dataset=="ranaddrad") %>%
  ggplot(aes(x=sampling_depth, y=percent_missing, group=sample_ID, color=sample_ID)) +
  geom_line(size=2) +
  geom_point(size=7) +
  ylab("Missing data (%)") +
  xlab("Sampling depth") +
  theme(axis.text=element_text(size=22),
        axis.title=element_text(size=25),
        axis.line.x = element_line(size=2),
        axis.line.y = element_line(size=2),
        axis.ticks.x = element_line(size=2),
        axis.ticks.y = element_line(size=2),
        strip.background = element_blank()) +
  ggtitle("Rana") +
  facet_grid(~dataset) +
  scale_color_manual(values=palette) +
  scale_y_continuous(limits=c(0,100))

## Add in avg amounts of missing data on top
p.rana_MD <-
  p.ranaind +
  geom_line(data=new_rana, aes(y=avg_missing, x=sampling_depth, group=dataset), color="black", size=2) +
  geom_point(data=new_rana, aes(y=avg_missing, x=sampling_depth, group=dataset), color="black", size=7)


# Export and save plots ------------------------------------------------------------

plot <-
  plot_grid(p.rana_MD, p.epi_MD, nrow=2)
plot
ggsave("Fig4_missingdata.pdf", width=11.5, height = 17.9)
