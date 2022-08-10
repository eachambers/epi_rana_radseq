library(tidyverse)
library(cowplot)
library(ggplot2)
library(phylotools)

#' Process ddRAD data generated with Matz native pipeline
#' 
#' Takes output files from ddRAD dataset generated with Matz native pipeline and converts to Phylip
#' 
#' @param input_file_path is path to varsites and allsites files, should contain
#' taxon name (epi/rana), sampling depth, and dataset (2brad)
#' @param names_path is path to sample names (IN ORDER OF BAMS FILE)
#' 
#' @return matz2Phylip fxn saves Phylip file with same name as input_file_path
#' 
#' Before this code is run, make sure you have files that specify the sample
#' names in the order they appear in the sequence files. These can be found
#' in the X_bams file.

###########################################################################
###########################################################################

matz2Phylip <- function(input_file_path, names_path){
  
  # tibbles don't transpose quickly so we'll use baseR
  dat <- read.delim(input_file_path, header=F)
  names <- read.delim(names_path, header=F)
  
  # remove irrelevant first 2 cols
  dat <- 
    dat %>% 
    dplyr::select(3:14)
  
  # set sample IDs to cols
  colnames(dat) = names$V1
  
  # transpose and concatenate sequences
  dat <- dat %>% purrr::map(.f = ~paste0(., collapse = "")) %>% 
    bind_cols() %>% 
    t()
  
  dat <- dat %>% as.data.frame() %>% 
    # make col with sample names
    mutate(sample = row.names(dat)) %>% 
    # ensure cols are in correct order for phylip
    select(sample, V1)
  
  # print size of matrix (i.e., no. sites))
  print(paste("Number of sites in alignment =", unique(nchar(as.character(dat$V1)))), sep="")
  
  # convert and export
  dat2phylip(dat, outfile=(paste(input_file_path, ".phylip", sep="")))
}

###########################################################################
###########################################################################

input_file_path = "epiddrad_Matz/epi_ddrad_varsites"
names_path = "epiddrad_Matz/epi_ddrad_names"

# Run above functions
matz2Phylip(input_file_path, names_path)
