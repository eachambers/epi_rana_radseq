library(devtools)
library(tidyverse)
library(ggplot2)
library(cowplot)

theme_set(theme_cowplot())

## The following code generates Figs S2 and S6, which illustrates read depth among reciprocal
## pipelines (Fig S2) and read depth vs missing data (Fig S6). Avg read depth per sample was
## calculated using the calculate_depth_stats.R code

##    FILES REQUIRED:
##          data_files_input_into_scripts/missing_data_snps.txt
##          data_files_input_into_scripts/recip_depths.txt

# Import & process depth data -------------------------------------------------------------

depth <- read_tsv("5_Data_visualization/data_files_input_into_scripts/recip_depths.txt", col_names = TRUE)

# Remove replicate samples
depth <-
  depth %>% 
  filter(sample_ID != "Rber_T1113b", sample_ID != "Rchi_T2034b",
         sample_ID != "Eant_T6859b", sample_ID != "Ahah_R0089b")

# Ensure that things are ordered correctly (in phylo order) by changing levels
depth$sample_ID <- factor(depth$sample_ID, levels=c("Rber_T1113a", "Rber_T1114", 
                                                    "Rneo_T480", "Rneo_T527", "Rsph_T25870", 
                                                    "Rsph_T26064", "Rbla_D2864", "Rbla_D2865", 
                                                    "Rchi_T2034a", "Rchi_T2049", 
                                                    "Eant_T6859a", "Eant_T6857", 
                                                    "Etri_T6836", "Etri_T6842", "Ebou_R0153", 
                                                    "Ebou_R0156", "Snub_R0158", "Snub_R0159", 
                                                    "Ahah_R0089a", "Ahah_R0090"))

# Average read depth per ind per dataset for ipyrad-ddrad and matz-2brad
depth %>% 
  filter(concat_col=="epi_2brad_matz" |
         concat_col=="epi_ddrad_ipyrad" |
         concat_col=="rana_2brad_matz" |
         concat_col=="rana_ddrad_ipyrad") %>% 
  group_by(concat_col) %>% 
  summarize(mean_depth_dataset = mean(mean_depth_ind),
            min_depth = min(mean_depth_ind),
            max_depth = max(mean_depth_ind))

# Import & process missing data -------------------------------------------
miss <- read_tsv("5_Data_visualization/data_files_input_into_scripts/missing_data_snps.txt", col_names = TRUE)
miss <-
  miss %>%
  filter(sample_ID != "Rber_T1113b", sample_ID != "Rchi_T2034b",
         sample_ID != "Eant_T6859b", sample_ID != "Ahah_R0089b")

miss$sample_ID <- factor(miss$sample_ID, levels=c("Rber_T1113a", "Rber_T1114",
                                                  "Rneo_T480", "Rneo_T527", "Rsph_T25870",
                                                  "Rsph_T26064", "Rbla_D2864", "Rbla_D2865",
                                                  "Rchi_T2034a", "Rchi_T2049",
                                                  "Eant_T6859a", "Eant_T6857",
                                                  "Etri_T6836", "Etri_T6842", "Ebou_R0153",
                                                  "Ebou_R0156", "Snub_R0158", "Snub_R0159",
                                                  "Ahah_R0089a", "Ahah_R0090"))

# Define palette -----------------------------------------------------------

palette <- c("#BA85A2", "#6B2391", "#CA3518", "#C66E1F", 
             "#BA8517", "#D9C339", "#086955", "#7C935E", "#80A7B3", "#386784")

# Build Fig S2: read depths -------------------------------------------------------------

p_epi <-
  depth %>%
  filter(taxon == "epi") %>%
  ggplot(aes(x=sample_ID, y=mean_depth_ind, group=sample_ID, fill=sample_ID)) +
  # geom_line(size=2) +
  # geom_point(size=7) +
  geom_bar(stat = "identity") +
  ylab("Mean read depth/ind") +
  xlab("Sample") +
  theme(axis.text=element_text(size=14),
        axis.text.x = element_text(angle = 90),
        axis.title=element_text(size=18),
        axis.line.x = element_line(size=2),
        axis.line.y = element_line(size=2),
        axis.ticks.x = element_line(size=2),
        axis.ticks.y = element_line(size=2),
        strip.background = element_blank(),
        legend.position = "none") +
  ggtitle("Epipedobates") +
  facet_grid(~concat_col) +
  scale_fill_manual(values=palette)

# ===== do the same for Rana now =====
p_rana <-
  depth %>%
  filter(taxon == "rana") %>%
  ggplot(aes(x=sample_ID, y=mean_depth_ind, group=sample_ID, fill=sample_ID)) +
  # geom_line(size=2) +
  # geom_point(size=7) +
  geom_bar(stat = "identity") +
  ylab("Mean read depth/ind") +
  xlab("Sample") +
  theme(axis.text=element_text(size=14),
        axis.text.x = element_text(angle = 90),
        axis.title=element_text(size=18),
        axis.line.x = element_line(size=2),
        axis.line.y = element_line(size=2),
        axis.ticks.x = element_line(size=2),
        axis.ticks.y = element_line(size=2),
        strip.background = element_blank(),
        legend.position = "none") +
  ggtitle("Rana") +
  facet_grid(~concat_col) +
  scale_fill_manual(values=palette)


# Figure S6: read depth vs missing data -----------------------------------------------------------

# Remove irrelevant columns and rows
new_miss <-
  miss %>% 
  filter(sampling_depth == "total") %>% 
  select(-sampling_depth)

new_depth <-
  depth %>% 
  filter(dataset == "2brad" & method == "matz" |
         dataset == "ddrad" & method == "ipyrad") %>% 
  unite(dataset, "taxon", "dataset", sep="") %>% 
  select(-source, -method, -concat_col)
  
joined <- full_join(new_depth, new_miss)

epi_MD_p <-
  joined %>% 
  filter(taxon == "epi") %>% 
  ggplot(aes(x=mean_depth_ind, y=percent_missing, group=sample_ID, color=sample_ID)) +
  geom_point(size=3) +
  facet_grid(~dataset) +
  scale_color_manual(values = palette) +
  xlab("Mean read depth/ind") +
  ylab("Missing data (%)") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18)) +
  scale_y_continuous(limits=c(0,100))

rana_MD_p <-
  joined %>% 
  filter(taxon == "rana") %>% 
  ggplot(aes(x=mean_depth_ind, y=percent_missing, group=sample_ID, color=sample_ID)) +
  geom_point(size=3) +
  facet_grid(~dataset) +
  scale_color_manual(values = palette) +
  xlab("Mean read depth/ind") +
  ylab("Missing data (%)") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=18)) +
  scale_y_continuous(limits=c(0,100))

# Export plots ------------------------------------------------------------

# Export Fig S2 (read depth)
plot_grid(p_rana, p_epi, nrow=2)
ggsave("FigS2_readdepthrecip.pdf", width=11.5, height = 17.9)

# Export Fig S6 (read depth vs missing data)
plot_grid(rana_MD_p, epi_MD_p, nrow=2)
ggsave("FigS2_readdepthrecip.pdf", width=11.5, height = 17.9)

