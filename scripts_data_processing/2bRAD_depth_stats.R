setwd("~/Box Sync/Epipedobates project/BIOINFORMATICS/ANALYSIS/epi_rana_radseq")

library(ggplot2)
library(tidyverse)
library(cowplot)
library(viridis)

theme_set(theme_cowplot())

#' Calculate depth stats for data generated with Matz native pipeline from 
#' ANGSD output .depthSample file
#' 
#' @param input_file_path is path to .depthSample output file from angsd
#' @param names_path is path to sample names (IN ORDER OF BAMS FILE)
#' 
#' @return calcDepthSample fxn calculates average read depth per sample
#' @return calcMeanDepthDataset fxn calculates the mean read depth for a given dataset
#' @return calcDepthDataset fxn calculates the number of reads at each read depth
#' @return plotDepth fxn plots depth distribution of all reads in dataset
#' 
#' Before this code is run, make sure you have files that specify the sample
#' names in the order they appear in the sequence files. These can be found
#' in the X_bams file.


###########################################################################
###########################################################################

calcDepthSample <- function(input_file_path, names_path){
  names = read_tsv(names_path, col_names = FALSE)
  depthSample <-
    read_tsv(input_file_path, col_names = FALSE) %>% 
    as.data.frame()
  
  col_names = c(1:ncol(depthSample))
  rownames(depthSample) = names$X1
  colnames(depthSample) = col_names
  
  depthSample <-
    depthSample %>% 
    rownames_to_column(., "sample") %>% 
    pivot_longer(-sample, names_to = "read_depth", values_to = "reads")
  
  depthSample$read_depth <- as.integer(depthSample$read_depth)
  
  depthSample <-
    depthSample %>% 
    group_by(sample) %>%
    summarize(corr_reads = sum(read_depth*reads, na.rm=T),
              total_reads = sum(reads, na.rm=T),
              avg_depth_sample = sum(corr_reads/total_reads)) %>% 
    dplyr::select(-corr_reads, -total_reads)
  
  return(depthSample)
}

calcMeanDepthDataset <- function(depthSample){
  depthAvg <-
    depthSample %>% 
    summarize(avg_depth = mean(avg_depth_sample))
  
  return(depthAvg)
}

calcDepthDataset <- function(input_file_path, names_path){
  names = read_tsv(names_path, col_names = FALSE)
  depthSample <-
    read_tsv(input_file_path, col_names = FALSE) %>% 
    as.data.frame()
  
  col_names = c(1:ncol(depthSample))
  rownames(depthSample) = names$X1
  colnames(depthSample) = col_names
  
  depthSample <-
    depthSample %>% 
    rownames_to_column(., "sample") %>% 
    pivot_longer(-sample, names_to = "read_depth", values_to = "reads")
  
  depthDataset <-
    depthSample %>% 
    group_by(read_depth) %>% 
    summarize(total_reads = sum(reads)) %>% 
    drop_na()
  
  return(depthDataset)
}

plotDepth <- function(depthDataset){
  depthDataset %>% 
    ggplot(aes(x=read_depth, y=total_reads)) +
    geom_bar(stat="identity") +
    ggtitle(input_file_path)
}

###########################################################################
###########################################################################

input_file_path = "epi_ddrad_depth.depthSample"
names_path = "epi_ddrad_names"
# input_file_path = "rana_total_depth.depthSample"
# names_path = "rana_names"

depthSample <- calcDepthSample(input_file_path, names_path)
depthAvg <- calcMeanDepthDataset(depthSample)
depthDataset <- calcDepthDataset(input_file_path, names_path)

plotDepth(depthDataset)

