library(tidyverse)
library(cowplot)
library(viridis)

theme_set(theme_cowplot())

#     This code contains four functions to process read depth files
#     outputted by ANGSD.
#     1. calc_depth_sample() calculates average read depth per sample
#     2. calc_mean_depth_dataset() calculates the mean read depth for a dataset
#     3. calc_depth_reads() calculates the number of reads at each read depth
#     4. plot_depth() plots depth distribution of all reads in dataset
#
#     Before this code is run, make sure you have files that specify the sample
#     names in the order they appear in the sequence files. These can be found
#     in the X_bams file.


###########################################################################
###########################################################################

#' Calculate average read depth per sample
#'
#' @param input_file_path path to ANGSD .depthSample file
#' @param names_path path to taxon_names file
#'
#' @return avg_depth_sample, the average depth per sample
#'
calc_depth_sample <- function(input_file_path, names_path){
  names = read_tsv(names_path, col_names = FALSE)
  avg_depth_sample <-
    read_tsv(input_file_path, col_names = FALSE) %>% 
    as.data.frame()
  
  col_names = c(1:ncol(avg_depth_sample))
  rownames(avg_depth_sample) = names$X1
  colnames(avg_depth_sample) = col_names
  
  avg_depth_sample <-
    avg_depth_sample %>% 
    rownames_to_column(., "sample") %>% 
    pivot_longer(-sample, names_to = "read_depth", values_to = "reads")
  
  avg_depth_sample$read_depth <- as.integer(avg_depth_sample$read_depth)
  
  avg_depth_sample <-
    avg_depth_sample %>% 
    group_by(sample) %>%
    summarize(corr_reads = sum(read_depth*reads, na.rm=T),
              total_reads = sum(reads, na.rm=T),
              avg_depth_sample = sum(corr_reads/total_reads)) %>% 
    dplyr::select(-corr_reads, -total_reads)
  
  return(avg_depth_sample)
}

#' Calculate mean depth of the dataset (across all samples) and print out value
#'
#' @param avg_depth_sample object returned from calc_depth_sample() function
#'
calc_mean_depth_dataset <- function(avg_depth_sample){
  avg_depth_dataset <-
    avg_depth_sample %>% 
    summarize(avg_depth = mean(avg_depth_sample))
  
  print(paste0("Mean sample depth of dataset: ", avg_depth_dataset, sep=""))
}

#' Calculate number of reads at each read depth of the dataset
#'
#' @param input_file_path 
#' @param names_path path to taxon_names file
#'
#' @return depth_dataset
#'
calc_depth_reads <- function(input_file_path, names_path){
  names = read_tsv(names_path, col_names = FALSE)
  depth_sample <-
    read_tsv(input_file_path, col_names = FALSE) %>% 
    as.data.frame()
  
  col_names = c(1:ncol(depth_sample))
  rownames(depth_sample) = names$X1
  colnames(depth_sample) = col_names
  
  depth_sample <-
    depth_sample %>% 
    rownames_to_column(., "sample") %>% 
    pivot_longer(-sample, names_to = "read_depth", values_to = "reads")
  
  depth_dataset <-
    depth_sample %>% 
    group_by(read_depth) %>% 
    summarize(total_reads = sum(reads)) %>% 
    drop_na()
  
  depth_dataset$read_depth <- as.numeric(depth_dataset$read_depth)
  
  return(depth_dataset)
}

#' Plot depth statistics
#'
#' @param depth_dataset object returned from calc_depth_reads() function
#'
#' @return bar plot with number of reads at each depth
#'
plot_depth <- function(depth_dataset){
  max_depth <- max(depth_dataset$read_depth)
  
  # add buffer for plotting
  max_depth <- max_depth*1.1
  
  depth_dataset %>% 
    ggplot(aes(x=read_depth, y=total_reads)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(angle = 90, size = 10)) +
    scale_x_continuous(expand = c(0,0), limits = c(0, max_depth)) +
    scale_y_continuous(expand = c(0,0)) +
    xlab("Read depth") +
    ylab("Number of reads")
}

###########################################################################
###########################################################################

# Run functions, assuming you're in the main supplementary materials folder

input_file_path = "1_Bioinformatics/Matz_2bRAD/rana2brad_Matz/rana_total_depth.depthSample"
names_path = "1_Bioinformatics/Matz_2bRAD/rana2brad_Matz/rana_names"

avg_depth_sample <- calc_depth_sample(input_file_path, names_path)
avg_depth_sample
calc_mean_depth_dataset(avg_depth_sample)
depth_dataset <- calc_depth_reads(input_file_path, names_path)

plot_depth(depth_dataset)

