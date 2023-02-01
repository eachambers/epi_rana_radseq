library(tidyverse)
library(cowplot)
library(phylotools)

#     This code contains three functions to process Matz Lab output files from
#     2bRAD data:
#     1. matz2Phylip() takes output files from Matz Lab pipeline and converts
#        to Phylip file
#     2. calcReps() calculates shared sites (and loci) between replicate samples
#     3. summaryReps() summarizes shared sites between replicate samples for plotting
#
#     Before this code is run, make sure you have files that specify the sample
#     names in the order they appear in the sequence files. These can be found
#     in the X_bams file.


###########################################################################
###########################################################################

#' Convert Matz Lab output and convert to Phylip file
#'
#' @param input_file_path path to varsites and allsites files
#' @param names_path path to sample names (IN ORDER OF BAMS FILE)
#'
#' @return saves Phylip file with same name as input_file_path
#' @export
#'
#' @examples
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

#' Calculate shared loci between sets of replicate samples
#'
#' @param retab_file_path path to retab file generated with retab script
#' @param retab_names_path path to retab names (same as names except contains tag and site)
#'
#' @return df with tags & sites between replicate samples
#' @export
#'
#' @examples
calcReps <- function(retab_file_path, retab_names_path){
  retab <- read.delim(retab_file_path, header=F)
  names <- read.delim(names_path, header=F)
  
  # rename cols as tag, site, and sample
  colnames(retab) = names$V1
  
  # tidy data
  retab <- 
    retab %>% 
    pivot_longer(-c(tag, site), names_to = "sample", values_to = "base") %>% 
    # filter for relevant taxa and add species name col
    filter(sample=="Rber_T1113a" | sample=="Rber_T1113b" | 
             sample=="Rchi_T2034a" | sample=="Rchi_T2034b" |
             sample=="Eant_T6859a" | sample=="Eant_T6859b" |
             sample=="Ahah_R0089a" | sample=="Ahah_R0089b") %>%
    mutate(species = if_else(grepl("Rber_T1113",sample), "Rber",
                             if_else(grepl("Eant_T6859",sample),"Eant",
                                     if_else(grepl("Ahah_R0089",sample),"Ahah",
                                             if_else(grepl("Rchi_T2034", sample), "Rchi", NA_character_)))))
  # combine tag and site number
  retab <- retab_new %>% unite("fullsite", tag:site, sep="_")
    
  return(retab)
}

#' Summarize shared loci between replicate samples for repeatability analysis
#'
#' @param retab output retab obj from calcReps function
#' @param retab_file_path path to retab file generated with retab script
#'
#' @return number of shared sites between replicate samples for SNP datasets
#' @export
#'
#' @examples
summaryReps <- function(retab, retab_file_path){
  summary <-
    retab %>% 
    as_tibble() %>% 
    group_by(species, fullsite) %>% 
    # give t/f statements about whether loci are shared or not
    # and remove loci that aren't shared
    summarize(shared = all(!grepl("N", base)),
              to_remove = all(grepl("N", base))) %>% 
    filter(!to_remove) %>% 
    group_by(species) %>% 
    # calculate no. shared loci
    summarize(number_loci_shared = sum(shared))
  
  loci_per_sample <-
    retab %>% 
    as_tibble() %>%
    group_by(sample, fullsite) %>%
    summarize(total_persample = all(!grepl("N", base))) %>%
    group_by(sample) %>%
    summarize(total_loci_sample = sum(total_persample)) %>% 
    mutate(species = if_else(grepl("Rber", sample), "Rber",
                             if_else(grepl("Eant", sample), "Eant",
                                     if_else(grepl("Ahah", sample), "Ahah",
                                             if_else(grepl("Rchi", sample), "Rchi", NA_character_)))))
  
  # join two stats together
  final_summary <-
    full_join(summary, loci_per_sample) %>% 
    mutate(prop_loci_within = number_loci_shared/total_loci_sample,
           samp_depth = ifelse(grepl("t1", retab_file_path), "t1",
                               ifelse(grepl("t2", retab_file_path), "t2",
                                      ifelse(grepl("t3", retab_file_path), "t3",
                                             ifelse(grepl("total", retab_file_path), "total", NA_character_)))),
           taxon = ifelse(grepl("Rana", retab_file_path), "rana",
                          ifelse(grepl("Epi", retab_file_path), "epi", 
                                 ifelse(grepl("rana", retab_file_path), "rana",
                                        ifelse(grepl("epi", retab_file_path), "epi", NA_character_)))),
           dataset = ifelse(grepl("2brad", retab_file_path), "2brad",
                            ifelse(grepl("ddrad", retab_file_path), "ddrad", NA_character_)))
}

###########################################################################
###########################################################################

# Run above functions
# Do for all sampling depths (t1, t2, t3, total) and for each Epi and Rana
input_file_path = "epi2brad_Matz/epi_t1_varsites"
names_path = "epi2brad_Matz/epi_names"
retab_file_path = "epi2brad_Matz/epi_t1_snps_retab"
retab_names_path = "epi2brad_Matz/epi_names_retab.txt"

matz2Phylip(input_file_path, names_path) # saves phylip file to wd
retab <- calcReps(retab_path_path, retab_names_path)
epit1 <- summaryReps(retab, retab_file_path)

###########################################################################
###########################################################################

# If you run summaryReps on each sampling depth and taxon and save as diff datasets, can do
# the following to combine and export
all <-
  rbind(ranat1, ranat2, ranat3, ranatotal,
        epit1, epit2, epit3, epitotal)

write.csv(all, "../Analyses/2bRAD_shared_loci_replicates.csv")

  
