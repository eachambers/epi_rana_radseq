library(devtools)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(here)

theme_set(theme_cowplot())

## The following code generates Fig. 4, which illustrate read depth among reciprocal
## pipelines and missing data. Avg read depth per sample was calculated using the
## calculate_depth_stats.R code, and missing data was calculated in PAUP*
## using the `missdata` function.

##    FILES REQUIRED:
##          data_files_input_into_scripts/readdepth_missingdata_snps.txt


# Import & process data ---------------------------------------------------

dat <- read_tsv(here("data_files_input_into_scripts", "readdepth_missingdata_snps.txt"), col_names = TRUE)

# Remove replicate samples
dat <-
  dat %>% 
  filter(sample_ID != "Rber_T1113b", sample_ID != "Rchi_T2034b",
         sample_ID != "Eant_T6859b", sample_ID != "Ahah_R0089b")

# Ensure that things are ordered correctly (in phylo order) by changing levels
dat$sample_ID <- factor(dat$sample_ID, levels=c("Rber_T1113a", "Rber_T1114", 
                                                    "Rneo_T480", "Rneo_T527", "Rsph_T25870", 
                                                    "Rsph_T26064", "Rbla_D2864", "Rbla_D2865", 
                                                    "Rchi_T2034a", "Rchi_T2049", 
                                                    "Eant_T6859a", "Eant_T6857", 
                                                    "Etri_T6836", "Etri_T6842", "Ebou_R0153", 
                                                    "Ebou_R0156", "Snub_R0158", "Snub_R0159", 
                                                    "Ahah_R0089a", "Ahah_R0090"))

# Add concatenated columns containing info on taxon, method, and dataset (e.g., epi + 2brad + ipyrad)
# and calculate average read depth per ind per dataset for ipyrad-ddrad and matz-2brad
dat <-
  dat %>% 
  # add concat columns with info on taxon, dataset, and method (e.g., epi_2brad_ipyrad)
  unite(concat_col, c("taxon", "dataset", "method"), remove = FALSE)


# Calculate mean depth per dataset ----------------------------------------

dat %>% 
  # only look at total sampling depth
  filter(sampling_depth == "total") %>% 
  # only consider "correct" (i.e., not reciprocal) pipelines
  filter(concat_col=="epi_2brad_matz" |
         concat_col=="epi_ddrad_ipyrad" |
         concat_col=="rana_2brad_matz" |
         concat_col=="rana_ddrad_ipyrad") %>% 
  # calculate per dataset across all 12 inds
  group_by(concat_col) %>% 
  summarize(mean_depth_dataset = mean(mean_depth_ind),
            min_depth = min(mean_depth_ind),
            max_depth = max(mean_depth_ind))


# Import & process missing data -------------------------------------------
# miss <- read_tsv("5_Data_visualization/data_files_input_into_scripts/missing_data_snps.txt", col_names = TRUE)
# miss <-
#   miss %>%
#   filter(sample_ID != "Rber_T1113b", sample_ID != "Rchi_T2034b",
#          sample_ID != "Eant_T6859b", sample_ID != "Ahah_R0089b")
# 
# miss$sample_ID <- factor(miss$sample_ID, levels=c("Rber_T1113a", "Rber_T1114",
#                                                   "Rneo_T480", "Rneo_T527", "Rsph_T25870",
#                                                   "Rsph_T26064", "Rbla_D2864", "Rbla_D2865",
#                                                   "Rchi_T2034a", "Rchi_T2049",
#                                                   "Eant_T6859a", "Eant_T6857",
#                                                   "Etri_T6836", "Etri_T6842", "Ebou_R0153",
#                                                   "Ebou_R0156", "Snub_R0158", "Snub_R0159",
#                                                   "Ahah_R0089a", "Ahah_R0090"))


# Define palette -----------------------------------------------------------

palette <- c("#BA85A2", "#6B2391", "#CA3518", "#C66E1F", 
             "#BA8517", "#D9C339", "#086955", "#7C935E", "#80A7B3", "#386784")


# Build Fig XX: read depths and missing data ------------------------------

p_epi <-
  dat %>%
  filter(taxon == "epi") %>%
  # only consider total sampling depth
  filter(sampling_depth == "total") %>% 
  ggplot(aes(x = sample_ID, y = mean_depth_ind, group = sample_ID, fill = sample_ID)) +
  geom_bar(stat = "identity") +
  ylab("Mean read depth/ind") +
  xlab("Sample") +
  theme(axis.text = element_text(size = 25),
        # axis.text.x = element_text(angle = 90),
        axis.text.x = element_blank(),
        # axis.title = element_text(size = 18),
        axis.title = element_blank(),
        axis.line = element_line(linewidth = 1),
        axis.ticks.y = element_line(linewidth = 1),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") +
  ggtitle("Epipedobates") +
  facet_grid(~concat_col) +
  scale_fill_manual(values = palette)

# Now, add on missing data proportions as a second axis
# First, figure out how much to scale the second axis units (props in this case)
dat %>% 
  filter(taxon == "epi" & sampling_depth == "total") %>% 
  summarize(max(mean_depth_ind)) # max is 40.2, so we'll need to multiply percent_missing by 0.402 (40.2/100)

p_epi <-
  p_epi +
  geom_line(aes(y = percent_missing*0.402, group = concat_col), color = "darkgrey", linewidth = 1) +
  geom_point(aes(y = percent_missing*0.402), size = 4, color = "darkgrey") +
  facet_grid(~concat_col)

p_epi <-
  p_epi +
  scale_y_continuous(sec.axis = sec_axis(~./0.402))


# ===== do the same for Rana now =====
p_rana <-
  dat %>%
  filter(taxon == "rana") %>%
  # only consider total sampling depth
  filter(sampling_depth == "total") %>% 
  ggplot(aes(x = sample_ID, y = mean_depth_ind, group = sample_ID, fill = sample_ID)) +
  geom_bar(stat = "identity") +
  ylab("Mean read depth/ind") +
  xlab("Sample") +
  theme(axis.text = element_text(size = 25),
        # axis.text.x = element_text(angle = 90),
        axis.text.x = element_blank(),
        # axis.title = element_text(size = 18),
        axis.title = element_blank(),
        axis.line = element_line(linewidth = 1),
        axis.ticks.y = element_line(linewidth = 1),
        axis.ticks.x = element_blank(),
        strip.background = element_blank(),
        legend.position = "none") +
  ggtitle("Rana") +
  facet_grid(~concat_col) +
  scale_fill_manual(values = palette)

# Now, add on missing data proportions as a second axis
# First, figure out how much to scale the second axis units (props in this case)
dat %>% 
  filter(taxon == "rana" & sampling_depth == "total") %>% 
  summarize(max(mean_depth_ind)) # max is 45.9, so we'll need to multiply percent_missing by 0.459 (45.9/100)

p_rana <-
  p_rana +
  geom_line(aes(y = percent_missing*0.459, group = concat_col), color = "darkgrey", linewidth = 1) +
  geom_point(aes(y = percent_missing*0.459), size = 4, color = "darkgrey") +
  facet_grid(~concat_col)

p_rana <-
  p_rana +
  scale_y_continuous(sec.axis = sec_axis(~./0.459))


# Figure S6: read depth vs missing data -----------------------------------------------------------

# Remove irrelevant columns and rows
# new_miss <-
#   miss %>% 
#   filter(sampling_depth == "total") %>% 
#   select(-sampling_depth)
# 
# new_depth <-
#   depth %>% 
#   filter(dataset == "2brad" & method == "matz" |
#          dataset == "ddrad" & method == "ipyrad") %>% 
#   unite(dataset, "taxon", "dataset", sep="") %>% 
#   select(-source, -method, -concat_col)
#   
# joined <- full_join(new_depth, new_miss)
# 
# epi_MD_p <-
#   joined %>% 
#   filter(taxon == "epi") %>% 
#   ggplot(aes(x=mean_depth_ind, y=percent_missing, group=sample_ID, color=sample_ID)) +
#   geom_point(size=3) +
#   facet_grid(~dataset) +
#   scale_color_manual(values = palette) +
#   xlab("Mean read depth/ind") +
#   ylab("Missing data (%)") +
#   theme(axis.text=element_text(size=14),
#         axis.title=element_text(size=18)) +
#   scale_y_continuous(limits=c(0,100))
# 
# rana_MD_p <-
#   joined %>% 
#   filter(taxon == "rana") %>% 
#   ggplot(aes(x=mean_depth_ind, y=percent_missing, group=sample_ID, color=sample_ID)) +
#   geom_point(size=3) +
#   facet_grid(~dataset) +
#   scale_color_manual(values = palette) +
#   xlab("Mean read depth/ind") +
#   ylab("Missing data (%)") +
#   theme(axis.text=element_text(size=14),
#         axis.title=element_text(size=18)) +
#   scale_y_continuous(limits=c(0,100))


# Combine & export plots --------------------------------------------------

plot_grid(p_rana, p_epi, nrow = 2)
ggsave("Fig4_readdepthMD.pdf", width = 17.9, height = 11.5)
