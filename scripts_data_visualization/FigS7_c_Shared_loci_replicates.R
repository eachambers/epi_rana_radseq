setwd("~/Box Sync/Epipedobates project/SuppMaterials/5_Data_visualization")

library(stringr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)

theme_set(theme_cowplot())

## The following code constructs Figure S7, which illustrates the proportion of shared loci
## between replicate samples. There are two replicates in each dataset (Ahah_R0089a&b and
## Eant_T6859a&b for Epipedobates; Rber_T1113a&b and Rchi_T2034a&b for Rana), and these
## analyses were done for all four sampling depths (t1, t2, t3, total).

##    FILES REQUIRED:
##          ../2_Bioinformatics/iPyrad_ddRAD/sampling_depth_summary_data/samplingdepth_ddrad_sumstats.csv [made using extract_data.ipynb]
##          data_files_input_into_scripts/ddRAD_shared_loci_replicates.csv [made using Split_loci_files.ipynb & Shared_loci_replicates.ipynb]
##          data_files_input_into_scripts/2bRAD_shared_loci_replicates.csv [made using 2bRADnative_processdata.R]

# Import data -------------------------------------------------------------

## The following data have info on the number of shared loci between replicate samples:
rep_ddrad <- read_csv("data_files_input_into_scripts/ddRAD_shared_loci_replicates.csv")
rep_2brad <- read_csv("data_files_input_into_scripts/2bRAD_shared_loci_replicates.csv") # already contains total loci


# Tidy the ddRAD data -----------------------------------------------------

# Our 2bRAD data is complete with all required columns, but we need to fill in 
# some info (total loci per sample and prop shared loci) for our ddRAD data

# The number of total loci per sample is in sumstats:
sumstats <- read_csv("../2_Bioinformatics/iPyrad_ddRAD/sampling_depth_summary_data/samplingdepth_ddrad_sumstats.csv")

# Tidy up data ------------------------------------------------------------

# Remove unwanted columns; we're only interested in loci in assembly column and
# also remove unwanted samples; we only want replicate samples
# Finally, append a col that specifies that this is for ddRAD data only
loci_ass <-
sumstats %>%
  dplyr::select(sample, loci_assembly, file_name) %>%
  filter(sample=="Ahah_R0089a" | sample=="Ahah_R0089b" | 
           sample=="Eant_T6859a" | sample=="Eant_T6859b" |
           sample=="Rchi_T2034a" | sample=="Rchi_T2034b" |
           sample=="Rber_T1113a" | sample=="Rber_T1113b") %>%
  mutate(samp_depth = if_else(grepl("t1",file_name), "t1",
                           if_else(grepl("t2",file_name),"t2",
                                   if_else(grepl("t3",file_name),"t3",
                                           if_else(grepl("total", file_name), "total", NA_character_))))) %>%
  mutate(dataset="ddrad") %>%
  dplyr::select(-file_name) %>% 
  mutate(total_loci_sample = loci_assembly) %>% 
  select(-loci_assembly)

# Join datasets together and calculate prop shared loci
rep_ddrad <-
  full_join(rep_ddrad, loci_ass) %>% 
  mutate(prop_loci_within = number_loci_shared/total_loci_sample)
 

# Combine 2bRAD and ddRAD datasets ----------------------------------------

replicates <-
  full_join(rep_2brad, rep_ddrad)
  
# One last calculation!
replicates <-
  replicates %>% 
  group_by(taxon, dataset, samp_depth, species) %>% 
  mutate(union_species=(sum(total_loci_sample))-number_loci_shared) %>% 
  mutate(prop_total=(number_loci_shared/union_species)*100) %>% 
  ungroup()

# Build plots ------------------------------------------

# Change levels so they're ordered consistently
# replicates$sample <- factor(replicates$sample, levels = c("Eant_T6859a", "Eant_T6859b", "Ahah_R0089a", "Ahah_R0089b", 
#                                     "Rber_T1113a", "Rber_T1113b", "Rchi_T2034a", "Rchi_T2034b"))

replicates$species <- factor(replicates$species, levels = c("Eant", "Ahah", "Rber", "Rchi"))

p.epi_reps <-
  replicates %>%
  filter(taxon=="epi") %>%
  ggplot(aes(x=samp_depth, y=prop_total, group=species, color=species)) +
  geom_point(size=2.5) +
  facet_wrap(~dataset, scales="free") +
  geom_line(size=.75) +
  geom_hline(yintercept = 100, color="grey", linetype="dashed", size=.75) +
  xlab("Sampling depth") +
  ylab("Proportion of shared loci (%)") +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=10),
        strip.background = element_blank(),
        legend.position = c(0.8,0.15),
        legend.title = element_blank(),
        legend.text = element_text(face="italic",size=10)) +
  scale_y_continuous(limits=c(28,100)) +
  scale_color_manual(values=c("#BA85A2", "#386784"), # "#8D529E", "#BA85A2", "#80A7B3", "#386784"
                     labels=c("E. anthonyi",
                              "A. hahneli")) +
  ggtitle("Epipedobates")

p.rana_reps <-
  replicates %>%
  filter(taxon=="rana") %>%
  ggplot(aes(x=samp_depth, y=prop_total, group=species, color=species)) +
  geom_point(size=2.5) +
  facet_wrap(~dataset, scales="free") +
  geom_line(size=.75) +
  geom_hline(yintercept = 100, color="grey", linetype="dashed", size=.75) +
  xlab("Sampling depth") +
  ylab("Proportion of shared loci (%)") +
  theme(axis.text=element_text(size=9),
        axis.title=element_text(size=10),
        strip.background = element_blank(),
        legend.position = c(0.8,0.15),
        legend.text = element_text(face="italic",size=10),
        legend.title = element_blank()) +
  scale_y_continuous(limits=c(28,100)) +
  scale_color_manual(values=c("#BA85A2", "#386784"),
                     labels=c("R. berlandieri",
                              "R. chiricahuensis")) +
  ggtitle("Rana")


# Save plots ------------------------------------------------

p.replicates <- plot_grid(p.rana_reps, p.epi_reps, nrow=2, align='v')
save_plot("FigS6_Shared_loci_replicates.pdf", p.replicates, base_width = 6.9, base_height = 6.8)
