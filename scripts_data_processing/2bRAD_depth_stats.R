setwd("~/Box Sync/Epipedobates project/BIOINFORMATICS/ANALYSIS/epi_rana_radseq")

library(ggplot2)
library(tidyverse)
library(cowplot)

## This code does the following:
##     Calculates the average read depth across samples from 2bRAD data, and
##     plots the distribution of reads with given read depths.

##    FILES REQUIRED:
##          2brad_depth.txt (concatenated data from rana_total_depth.depthSample &
##          epi_total_depth.depthSample; sample IDs given in x_total_bams files)


# Process and tidy data ---------------------------------------------------

depth <- read_tsv("2brad_depth.txt")

# Tidy data into long format
depth <-
  depth %>%
  gather(key="read_depth", value="value", -taxon, -sample)

depth$read_depth <- as.integer(depth$read_depth)

# Calculate avg read depth across reads, per sample and per taxon
sample_avg_epi <-
  depth %>%
  filter(taxon=="epi") %>%
  group_by(sample) %>%
  summarize(avg_depth=sum(read_depth*value/sum(value)))

sample_avg_rana <-
  depth %>%
  filter(taxon=="rana") %>%
  group_by(sample) %>%
  summarize(avg_depth=sum(read_depth*value/sum(value)))

# Now, calculate the avg read depth across all samples per taxon
sample_avg_epi %>%
  summarize(mean(avg_depth)) # 9.37

sample_avg_rana %>%
  summarize(mean(avg_depth)) # 8.81

# Calculate sum of reads at certain depths for plotting purposes
sum_epi <-
  depth %>%
    filter(taxon=="epi") %>%
    group_by(read_depth) %>%
    summarize(total_reads=sum(value))

sum_rana <-
  depth %>%
  filter(taxon=="rana") %>%
  group_by(read_depth) %>%
  summarize(total_reads=sum(value))


# Build figure ------------------------------------------------------------

# To look at the distribution of sample depths
sum_epi %>%
  ggplot(aes(x=read_depth, y=total_reads)) +
  geom_bar(stat="identity") +
  ggtitle("Epi")

sum_rana %>%
  ggplot(aes(x=read_depth, y=total_reads)) +
  geom_bar(stat="identity") +
  ggtitle("Rana")
