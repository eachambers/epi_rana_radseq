library(vcfR)

#' Calculates mean read depth per sample from vcf and prints avg depth (per sample) across dataset
#'
#' @param vcf_path path to vcf file
#'
#' @return mean_depth_sample, a df with average read depth per sample
#'
mean_sample_depth <- function(vcf_path){
  vcf <- read.vcfR(vcf_path)
  
  depth_stats <- extract.gt(vcf, element = "DP")
  rows_dp <- rownames(depth_stats)
  
  mean_depth_sample <-
    depth_stats %>% 
    as.data.frame(depth_stats, row.names = rows_dp) %>% 
    mutate_if(is.character, as.numeric) %>% 
    pivot_longer(cols = c(1:12), names_to = "sampleID", values_to = "reads") %>% 
    group_by(sampleID) %>% 
    summarize(mean_depth_sample = mean(reads))
  
  mean_depth_sample_ds <-
    mean_depth_sample %>% 
    summarize(mean = mean(mean_depth_sample))
  
  print(paste0("Average read depth per sample across dataset: ", mean_depth_sample_ds$mean, sep=""))
  
  return(mean_depth_sample)
}
