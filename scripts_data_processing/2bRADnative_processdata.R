setwd("~/Box Sync/Epipedobates project/SuppMaterials/2_Bioinformatics/Matz_native")

library(tidyverse)
library(cowplot)
library(ggplot2)
library(phylotools)

## This code does the following:
##      1. Takes output files from the Matz 2bRAD native pipeline (varsites
##      and allsites files) and (a) relabels sequences to correspond with
##      species and sample IDs; and (b) transposes the matrix such that
##      it is in Phylip format; and

##      2. Takes the output from the retabvcf.pl script (which assigns each
##      SNP to a tag) and calculates the number of shared sites between
##      replicate samples for the SNP datasets.

##    FILES REQUIRED:
##          X_varsites
##          X_allsites
##          epi_names
##          rana_names
##          epi_names_retab.txt
##          rana_names_retab.txt

# Before this code is run, make sure you have files that specify the sample
# names in the order they appear in the sequence files. These can be found
# in the X_bam file.

###########################################################################
# Transpose data files from 2bRAD native pipeline to make alignment --------


# For Epi files -----------------------------------------------------------

## Read in SNP files
# Importantly, we don't want to import these using tidyverse functions because
# tibbles don't transpose quickly with baseR
epivart1 <- read.delim("epi2brad_Matz/epi_t1_varsites", header=F)
epiallt1 <- read.delim("epi2brad_Matz/epi_t1_allsites", header=F)

epivart2 <- read.delim("epi2brad_Matz/epi_t2_varsites", header=F)
epiallt2 <- read.delim("epi2brad_Matz/epi_t2_allsites", header=F)

epivart3 <- read.delim("epi2brad_Matz/epi_t3_varsites", header=F)
epiallt3 <- read.delim("epi2brad_Matz/epi_t3_allsites", header=F)

epivartotal <- read.delim("epi2brad_Matz/epi_total_varsites", header=F)
epialltotal <- read.delim("epi2brad_Matz/epi_total_allsites", header=F)
epinames <- read.delim("epi2brad_Matz/epi_names", header=F)
epinames <- as.list(epinames) # Import names; this MUST be in same order as varsites

# Only select relevant columns (i.e. the ones with sequence data)
epivartotal <- epivartotal %>% select(3:14)
epialltotal <- epialltotal %>% select(3:14)

# Replace column names with sample IDs
epivartotal <- setNames(epivartotal, epinames$V1)
epialltotal <- setNames(epialltotal, epinames$V1)

# Transpose the file and concatenate the sequences into a single column
totalEpivar <- epivartotal %>% map(.f = ~paste0(., collapse = "")) %>% 
  bind_cols() %>% 
  t()

totalEpiall <- epialltotal %>% map(.f = ~paste0(., collapse = "")) %>% 
  bind_cols() %>% 
  t()

# Before we can convert the file to phylip, we need to make it into a dataframe (it's
# currently a list of characters), and also make an actual column with the sample IDs
# because they're currently actual row names and not their own column, so there will
# be issues when you go to convert to Phylip format
# Finally, you need to reverse the order of the columns because once you make a new
# column with the sample names it'll place this to the right of the sequence
totalEpivar <- totalEpivar %>% as.data.frame() %>% 
  mutate(sample = row.names(totalEpivar)) %>% 
  select(sample, V1)

totalEpiall <- totalEpiall %>% as.data.frame() %>% 
  mutate(sample = row.names(totalEpiall)) %>% 
  select(sample, V1)

# Verify size of matrix before export
unique(nchar(t1Epivar$V1)) # 27783
unique(nchar(t2Epivar$V1)) # 51073
unique(nchar(t3Epivar$V1)) # 59186
unique(nchar(totalEpivar$V1)) # 63070
unique(nchar(totalEpiall$V1)) # 3208050

# Convert to phylip format and export
dat2phylip(totalEpivar, outfile="Epi_total_SNPs.phylip")
dat2phylip(totalEpiall, outfile="Epi_total_allsites.phylip")

# For Rana files ----------------------------------------------------------

ranavart1 <- read.delim("rana2brad_Matz/rana_t1_varsites", header=F)
ranaallt1 <- read.delim("rana2brad_Matz/rana_t1_allsites", header=F)

ranavart2 <- read.delim("rana2brad_Matz/rana_t2_varsites", header=F)
ranaallt2 <- read.delim("rana2brad_Matz/rana_t2_allsites", header=F)

ranavart3 <- read.delim("rana2brad_Matz/rana_t3_varsites", header=F)
ranaallt3 <- read.delim("rana2brad_Matz/rana_t3_allsites", header=F)

ranavartotal <- read.delim("rana2brad_Matz/rana_total_varsites", header=F)
ranaalltotal <- read.delim("rana2brad_Matz/rana_total_allsites", header=F)

rananames <- read.delim("rana2brad_Matz/rana_names", header=F)
rananames <- as.list(rananames)

# Select only relevant columns
ranavartotal <- ranavartotal %>% select(3:14)
ranaalltotal <- ranaalltotal %>% select(3:14)

# Replace column names with sample IDs
ranavartotal <- setNames(ranavartotal, rananames$V1)
ranaalltotal <- setNames(ranaalltotal, rananames$V1)

# Transpose the file and concatenate the sequences into a single column
totalRanavar <- ranavartotal %>% map(.f = ~paste0(., collapse = "")) %>% 
  bind_cols() %>% 
  t()
totalRanaall <- ranaalltotal %>% map(.f = ~paste0(., collapse = "")) %>% 
  bind_cols() %>% 
  t()

totalRanavar <- totalRanavar %>% as.data.frame() %>% 
  mutate(sample = row.names(totalRanavar)) %>% 
  select(sample, V1)
totalRanaall <- totalRanaall %>% as.data.frame() %>% 
  mutate(sample = row.names(totalRanaall)) %>% 
  select(sample, V1)

# Convert to phylip format and export
dat2phylip(t1ranavar, outfile="Rana_t1_min4_SNPs.phylip")
dat2phylip(totalRanaall, outfile="Rana_total_allsites.phylip")


###########################################################################
###########################################################################

# Replicate analysis ------------------------------------------------------

## We use retab files for SNP datasets for each sampling depth
## Recall that the retab file has samples across the top row and loci ("tag") as each row
## We need to figure out which loci are shared between replicates

# Below is just for Rana t1; this must be done for Epi and all sampling depths
ranat1 <- read.delim("./rana2brad_Matz/rana_t1_snps_retab", header=F)
ranat2 <- read.delim("./rana2brad_Matz/rana_t2_snps_retab", header=F)
ranat3 <- read.delim("./rana2brad_Matz/rana_t3_snps_retab", header=F)
ranatotal <- read.delim("./rana2brad_Matz/rana_total_snps_retab", header=F)

epit1 <- read.delim("./epi2brad_Matz/epi_t1_snps_retab", header=F)
epit2 <- read.delim("./epi2brad_Matz/epi_t2_snps_retab", header=F)
epit3 <- read.delim("./epi2brad_Matz/epi_t3_snps_retab", header=F)
epitotal <- read.delim("./epi2brad_Matz/epi_total_snps_retab", header=F)


# Import names; this MUST be in same order as varsites!!!
rananames <- read.delim("./rana2brad_Matz/rana_names_retab.txt", header=F)
rananames <- as.list(rananames)

epinames <- read.delim("./epi2brad_Matz/epi_names_retab.txt", header=F)
epinames <- as.list(epinames)

# Tidy data ---------------------------------------------------------------

# Dataset you're looking at (RENAME AND RUN for each sampling depth and taxon)
dataset=ranat2
names=rananames

# Rename columns as samples
dataset <-
  dataset %>% 
  setNames(names$V1)

# Tidy data such that we have columns for tag # (locus), site (position on a given locus),
# sample, and base call
dataset <-
  dataset %>%
  gather(3:14, key="sample", value="base")

# Filter for relevant taxa and add species name column
rel_dataset <-
  dataset %>%
  filter(sample=="Rber_T1113a" | sample=="Rber_T1113b" | 
           sample=="Rchi_T2034a" | sample=="Rchi_T2034b" |
           sample=="Eant_T6859a" | sample=="Eant_T6859b" |
           sample=="Ahah_R0089a" | sample=="Ahah_R0089b") %>%
  mutate(species = if_else(grepl("Rber_T1113",sample), "Rber",
                           if_else(grepl("Eant_T6859",sample),"Eant",
                                   if_else(grepl("Ahah_R0089",sample),"Ahah",
                                           if_else(grepl("Rchi_T2034", sample), "Rchi", NA_character_))))) 

# Combines the tag and the site number; this ensures that when we calculate number of shared
# loci, we're looking at the exact same locaiton between samples (rather than just the same
# locus but different sites on that locus)
rel_dataset <-
  rel_dataset %>%
  unite("fullsite", tag:site, remove=FALSE)

# Summarize data ----------------------------------------------------------

# Our first summarize fxn (column 'shared') will give us true/false statements about
# whether loci are shared or need to be removed (because they're not shared); once 
# unshared loci are removed, the number of loci that are shared by replicate samples
# is calculated
shared <-
  rel_dataset %>%
  as_tibble() %>% 
  group_by(species, fullsite) %>%
  summarize(shared = all(!grepl("N", base)),
            to_remove = all(grepl("N", base))) %>% 
  filter(!to_remove) %>%
  group_by(species) %>% 
  summarize(number_loci_shared = sum(shared))


# Now, let's calculate the number of total loci per sample; also append a column
# that specifies the species to ease merging the two together
loci_per_sample <-
  rel_dataset %>%
  as_tibble() %>%
  group_by(sample, fullsite) %>%
  summarize(total_persample = all(!grepl("N", base))) %>%
  group_by(sample) %>%
  summarize(total_loci_sample = sum(total_persample))

loci_per_sample <-
  loci_per_sample %>%
  mutate(species = if_else(grepl("Rber", sample), "Rber",
                           if_else(grepl("Eant", sample), "Eant",
                                   if_else(grepl("Ahah", sample), "Ahah",
                                           if_else(grepl("Rchi", sample), "Rchi", NA_character_)))))

## BE SURE TO RENAME BELOW BASED ON DATASET AND TAXON!
# Combine the two datasets together and calculate the proportion of shared loci per sample
ranat1 <-
  full_join(shared, loci_per_sample) %>%
  mutate(prop_loci_within=number_loci_shared/total_loci_sample,
         samp_depth="t1", # MAKE SURE TO CHANGE THIS
         taxon="rana", # MAKE SURE TO CHANGE THIS
         dataset="2brad")

# Combine all sampling depths and taxa together
all <-
  rbind(ranat1, ranat2, ranat3, ranatotal,
        epit1, epit2, epit3, epitotal)

# Export data file -------------------------------------------------------

write.csv(all, "../Analyses/2bRAD_shared_loci_replicates.csv")
  
