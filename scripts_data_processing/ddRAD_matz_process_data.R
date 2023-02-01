library(tidyverse)
library(cowplot)
library(ggplot2)
library(phylotools)

#     This code contains two functions to process Matz Lab output files from
#     ddRAD data:
#     1. matz2Phylip() takes output files from Matz Lab pipeline and converts
#        to Phylip file
#     2. calc_basic_stats() calculates basic statistics: # sites, loci, SNPs, 
#        and some proportions

#     Before this code is run, make sure you have files that specify the sample
#     names in the order they appear in the sequence files. These can be found
#     in the X_bams file.

###########################################################################
###########################################################################

#' Takes output files from ddRAD dataset generated with Matz native pipeline and converts to Phylip
#' 
#' @param input_file_path is path to varsites and allsites files, should contain
#' @param names_path is path to sample names (IN ORDER OF BAMS FILE)
#' 
#' @return saves Phylip file with same name as input_file_path
#' 
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

#' Calculate and print basic statistics: # sites, loci, SNPs, and some proportions
#'
#' @param retab_path path to retab file
#' @param vcf_varsites path to varsites vcf file
#'
calc_basic_stats <- function(retab_path, vcf_varsites){
  vcf <- read.vcfR(vcf_varsites)
  retab <- read_tsv(retab_path, col_names = FALSE)
  
  nsites <- nrow(retab)
  
  nloci <-
    retab %>% 
    select(X1) %>% 
    unique() %>% 
    nrow()
  
  nsnps = nrow(extract.gt(vcf))
  
  print(paste0("Number of sites in dataset: ", nsites, sep=""))
  print(paste0("Number of SNPs in dataset: ", nsnps, sep=""))
  print(paste0("Number of loci in dataset: ", nloci, sep=""))
  print(paste0("SNPs per locus: ", nsnps/nloci, sep=""))
  print(paste0("SNPs per site: ", nsnps/nsites, sep=""))
}

###########################################################################
###########################################################################

input_file_path = "epiddrad_Matz/epi_ddrad_varsites"
names_path = "epiddrad_Matz/epi_ddrad_names"
vcf_varsites = "2_Bioinformatics/Matz_ddRAD/ranaddrad_Matz/rana_ddrad_varsites.vcf.gz"
retab_path = "2_Bioinformatics/Matz_ddRAD/ranaddrad_Matz/rana_ddrad_retab"

# Run above functions
matz2Phylip(input_file_path, names_path)
