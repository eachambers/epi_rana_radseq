library(tidyverse)
library(here)

## This code does the following:
##    The following takes in the sumstats files which are produced from our iPyrad pipeline
##    (they are a portion of the final stats file outputted directly from iPyrad) and processes
##    it such that it can be used downstream to generate Fig 3, the clustering threshold
##    comparison. It also combines ddRAD and 2bRAD results for Fig 3.

##    FILES REQUIRED (from DATA):
##          1_Bioinformatics/iPyrad_ddRAD/clust_threshold_summary_data/clustthreshold_ddrad_sumstats.csv
##          1_Bioinformatics/Matz_2bRAD/epi2brad_Matz/clust_threshold/epi2brad_clust_depthstats.tsv
##          1_Bioinformatics/Matz_2bRAD/epi2brad_Matz/clust_threshold/epi2brad_clust_locistats.tsv
##          1_Bioinformatics/Matz_2bRAD/rana2brad_Matz/clust_threshold/rana2brad_clust_depthstats.tsv
##          1_Bioinformatics/Matz_2bRAD/rana2brad_Matz/clust_threshold/rana2brad_clust_locistats.tsv


# Upload and process ddRAD data -------------------------------------------

# Upload files; sumstats files are generated using iPyrad output
ddrad_sumstats <- read_csv(here("1_Bioinformatics", "iPyrad_ddRAD", "clust_threshold_summary_data", "clustthreshold_ddrad_sumstats.csv"))

# Remove unnecessary columns
ddrad_sumstats <-
  ddrad_sumstats %>%
  select(sample, reads_consens, loci_assembly, file_name)

# Add columns with taxon and cluster threshold percent specified
ddrad_sumstats <-
  ddrad_sumstats %>%
  mutate(clust_threshold = case_when(grepl("80", file_name) ~ 80,
                                     grepl("81", file_name) ~ 81,
                                     grepl("82", file_name) ~ 82,
                                     grepl("83", file_name) ~ 83,
                                     grepl("84", file_name) ~ 84,
                                     grepl("85", file_name) ~ 85,
                                     grepl("86", file_name) ~ 86,
                                     grepl("87", file_name) ~ 87,
                                     grepl("88", file_name) ~ 88,
                                     grepl("89", file_name) ~ 89,
                                     grepl("90", file_name) ~ 90,
                                     grepl("91", file_name) ~ 91,
                                     grepl("92", file_name) ~ 92,
                                     grepl("93", file_name) ~ 93,
                                     grepl("94", file_name) ~ 94,
                                     grepl("95", file_name) ~ 95)) %>% 
  mutate(taxon = case_when(grepl("epi", file_name) ~ "epi",
                           grepl("rana", file_name) ~ "rana")) %>% 
  select(-file_name) %>% 
  mutate(method = "ddRAD")


# Upload and process 2bRAD data -------------------------------------------

epi2b_depth <- read_tsv(here("1_Bioinformatics", "Matz_2bRAD", "epi2brad_Matz", "clust_threshold", "epi2brad_clust_depthstats.tsv"), col_names = TRUE)
epi2b_loci <- read_tsv(here("1_Bioinformatics", "Matz_2bRAD", "epi2brad_Matz", "clust_threshold", "epi2brad_clust_locistats.tsv"), col_names = TRUE)

rana2b_depth <- read_tsv(here("1_Bioinformatics", "Matz_2bRAD", "rana2brad_Matz", "clust_threshold", "rana2brad_clust_depthstats.tsv"))
rana2b_loci <- read_tsv(here("1_Bioinformatics", "Matz_2bRAD", "rana2brad_Matz", "clust_threshold", "rana2brad_clust_locistats.tsv"), col_names = TRUE)

# Join taxa together, joining by sample and clust_threshold cols
epi2b_sumstats <- 
  left_join(epi2b_depth, epi2b_loci) %>% 
  mutate(taxon = "epi")

rana2b_sumstats <- 
  left_join(rana2b_depth, rana2b_loci) %>% 
  mutate(taxon = "rana")

# Join both taxa together and change col names to match ddRAD
twobrad_sumstats <- 
  bind_rows(epi2b_sumstats, rana2b_sumstats) %>% 
  mutate(method = "2bRAD") %>% 
  rename(reads_consens = total_reads,
         loci_assembly = total_loci_sample)


# Join 2bRAD with ddRAD ---------------------------------------------------

dat <- bind_rows(ddrad_sumstats, twobrad_sumstats) # 768 obs of 6 vars
write_tsv(dat, "clust_threshold_data.txt")
