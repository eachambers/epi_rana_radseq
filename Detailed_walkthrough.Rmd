---
title: "Detailed walkthrough"
author: "E. Anne Chambers"
date: "2023-02-01"
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Detailed walkthrough

The following is a detailed walkthrough of all analyses contained in this repository.

A few things to note:
* To investigate the effect of different sequencing depths on phylogenetic inference, we targeted a larger number of reads than is typical and then subsampled these. This resulted in four subsets of our data: *t1*, *t2*, *t3* and *total*. The *total* dataset is the total number of reads actually sequenced, and *t1*, *t2*, and *t3* were subsets of the *total* dataset, in increasing order.

### Parameter selection and basic data characteristics

#### Clustering threshold

#### Reciprocal bioinformatics pipelines

For the majority of our analyses (including all phylogenomic analyses), we used datasets that were generated using typical bioinformatics pipelines for those methods. For 2bRAD data, we used the Matz Lab pipeline (available [here](XXX)), and for ddRAD data, we used [iPyrad](XXX). However, there remained a question about how much using different bioinformatics pipelines might influence our downstream datasets, so we performed a reciprocal bioinformatics pipeline analysis in which we re-ran raw data through reciprocal pipelines: this meant using the Matz Lab pipeline to process ddRAD data, and iPyrad to process 2bRAD data.

The `ddRAD_matz_process_data.R` script processes raw data files produced from running ddRAD data through the Matz Lab pipeline. First, it uses the `matz2Phylip()` function to convert allsites and varsites files to Phylip format, and also contains the `calc_basic_stats()` function which takes in the retab file and varsites vcf file and summarizes basic statistics, including the number of sites, loci, SNPs, SNPs/locus, and SNPs/site.

