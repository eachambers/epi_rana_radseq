## Detailed walkthrough

The following is a detailed walkthrough of all analyses contained in this repository.

*A note on file naming:*

To investigate the effect of different sequencing depths on phylogenetic inference, we targeted a larger number of reads than is typical and then subsampled these. This resulted in four subsets of our data: *t1*, *t2*, *t3* and *total*. The *total* dataset is the total number of reads actually sequenced, and *t1*, *t2*, and *t3* were subsets of the *total* dataset, in increasing order. All data files found on the [Dryad repository](https://doi.org/10.5061/dryad.fbg79cnsp) for this project are named according to these sampling depths. For example, the "epiddrad_t2_4_outfiles/" directory found within the iPyrad ddRAD sampling depth outfiles corresponds to *Epipedobates* output files for the *t2* sampling depth for the ddRAD dataset. The 4 signifies the *min_ind_locus* parameter in iPyrad (all ddRAD datasets were generated with this parameter set to 4).

### General processing of raw output files

For **2bRAD data**, the `2bRADnative_processdata.R` script contains the `matz2Phylip()` function to convert varsites and allsites files to Phylip format. Phylip files were later used to reconstruct ML phylogenies.

For **ddRAD data**, the `extract_data.ipynb` script takes in raw iPyrad stats files produced after each (relevant) step, concatenates them, and generates:
* `X_sumstats.csv` contains read depth information *(will be used later to calculate read depth and perform clustering threshold analysis)*
* `X_coverage.csv` contains the number and sum of locus coverage
* `X_snpdist.csv` contains the total number of SNPs and PIs

All of the above files can be found on Dryad, within the "sampling_depth_summary_data" directories.

### Clustering threshold

We tested the clustering threshold, which is the percent similarity at which two sequences are considered orthologous and assigned to the same cluster (iPyrad parameter 14, *clust_threshold*; cd-hit-est within the Matz Lab pipeline). If this parameter is too high (too stringent), loci may be over-split, meaning that true homologs are interpreted as different loci; however, if the parameter is too low, loci may be under-split, i.e., paralogs incorrectly clustered into a single locus. We tested 16 clustering threshold values from 0.80–0.95 to assess the effect of this parameter on both steps.

The `clust_threshold_processing.R` script takes in a variety of output files from both 2bRAD and ddRAD datasets and combines them into a single dataframe, summarizing the number of consensus reads and loci in assembly for each of 16 clustering thresholds tested. This script generates the `clust_threshold_data.txt` input file for further data visualization.

***Data visualization***: `Fig3_Data_characterization.R` takes in the `clust_threshold_data.txt` input file and produces **Figure 3**.

### Read depth statistics

We also compared read depth between 2bRAD and ddRAD datasets, and its association with missing data.

For **2bRAD datasets**, use the `2bRAD_depth_stats.R` script:
* `calc_depth_sample()` function calculates the average read depth per sample
* `calc_depth_reads()` function calculates the mean depth per dataset
* `plot_depth()` function plots the distribution of read depth

For **ddRAD datasets**, `X_sumstats.csv` contains all information required for this analysis.

`calculate_depth_stats.R` can extract read depth directly from a vcf file.

***Data visualization***: `Fig4_read_depth.R` takes in the `readdepth_missingdata_snps.txt` input file and produces **Figure 4**.

### Phylogenetic signal

To quantify the amount of phylogenetic signal in our datasets, we calculated three different metrics using PAUP*: the number of parsimony-informative sites, the number and proportion of unambiguous synapomorphies on the tree ("shared sites"), and the retention index.

#### Retention index and parsimony-informative sites

***Data visualization***: `Fig6_Retention_index.R` takes in the `Retention_PIs.csv` input file and produces **Figure 6**.

#### Unambiguous synapomorphies or "shared sites"

Plotting these values was a little nuanced; we had four figures in the main manuscript and supplement that needed to be coded identically to ensure that values were mapped to the correct branches of the cladograms.

The `unambig_sums.txt` file contains raw data to be fed into this analysis. You also require `master.nexus`, which has the "master" topology upon which values will be plotted, and `Node_numbering_master_trees.png` to standardize the node numbering when building tree visualizations.

***Data visualization***: `Fig7&S1_PAUP_analysis.R` takes in the `unambig_sums.txt` input file and produces **Figures 7 and S1**.

### Missing data

Missing data was calculated in PAUP* on a per-individual basis and across datasets.

***Data visualization***: `Fig8_Missing_data.R` takes in the `readdepth_missingdata_snps.txt` input file and produces **Figure 8**.

### Allelic dropout

To understand the structure of allelic dropout in our data, we first inferred patterns of gains and losses of loci by analyzing the SNP data under Dollo parsimony. For each assembly, cells with non-missing nucleotide data were recoded as 1, or “present”, and cells with missing data were recoded as 0, or “absent”. Allelic dropout was quantified by counting the unambiguous changes from 1 to 0, using R scripts to parse PAUP* output from the command *describe / apolist chglist diag*.

***Data visualization***: `Fig9_Dollo_analysis.R` takes in the `Plot-Data-for-MS-FigS5.txt` input file and produces **Figure 9**. You will also require the `p.epi` and `p.rana` plot objects from `Fig7&S1_PAUP_analysis.R` which are the cladograms upon which these values will be plotted on to.

#### Significant vs. non-significant allelic dropout

Not all instances of dropout are equally informative about phylogeny. To determine whether an instance of allelic dropout has significant signal, we compared its expected number of changes on the tree for each locus (null expectation) to the observed number of changes with a chi-squared test, using the total datasets (see Supplemental Information for further explanation).

***Data visualization***: `Fig10_Recoded_significance_analysis.R` takes in the `recoded_signonsig.txt` input file and produces **Figure 10**. You will also require the `p.epi` and `p.rana` plot objects from `Fig7&S1_PAUP_analysis.R` which are the cladograms upon which these values will be plotted on to.

### Repeatability

We sequenced two samples twice (i.e., biological replicates; indicated by those taxa with asterisks in cladogram above) to determine how repeatable libraries were across both 2bRAD and ddRAD. To assess repeatability, we calculated the numbers of shared loci (and sites) between replicate samples.

For **2bRAD datasets**, the `2bRADnative_processdata.R` script contains a function, `calcReps()`, which will take in the retab file (generated using the retabvcf.pl script found [here](https://github.com/z0on/2bRAD_denovo/blob/master/retabvcf.pl)) and the path to `retab_names`, and will produce a dataframe with tags and sites between replicate samples. Within this same script, the `summaryReps()` function will summarize shared sites between replicates for plotting, `2bRAD_shared_loci_replicates.csv`.

For **ddRAD datasets**, generating summary data is a two-step process. First, use `Fig11_a_Split_loci_files.ipynb` to take in iPyrad `.loci` files and split them by locus (i.e., split by //). It will then output each locus that is shared by at least two individuals into a separate `file*.txt` file in order to eventually count the number of shared loci between replicate samples. Then, use `Fig11_b_Shared_loci_replicates.ipynb` to take the output from `Fig11_a_Split_loci.ipynb` script (contained within respective output folders for each taxon and sampling depth) and search for occurrences of replicates sharing a locus, counting those. This will produce the `ddRAD_shared_loci_replicates.csv` file.

***Data visualization***: `Fig11_c_Shared_loci_replicates.R` takes in the `X_sumstats.csv`, `2bRAD_shared_loci_replicates.csv`, and `ddRAD_shared_loci_replicates.csv` input files and produces **Figure 11**.

### Reciprocal bioinformatics pipelines

For the majority of our analyses (including all phylogenomic analyses), we used datasets that were generated using typical bioinformatics pipelines for those methods. For 2bRAD data, we used the Matz Lab pipeline (available [here](https://github.com/z0on/2bRAD_denovo)), and for ddRAD data, we used [iPyrad](https://ipyrad.readthedocs.io/en/master/). However, there remained a question about how much using different bioinformatics pipelines might influence our downstream datasets, so we performed a reciprocal bioinformatics pipeline analysis in which we re-ran raw data through reciprocal pipelines: this meant using the Matz Lab pipeline to process ddRAD data, and iPyrad to process 2bRAD data.

For **ddRAD data from Matz Lab pipeline**, the `ddRAD_matz_process_data.R` script processes raw data files produced from running ddRAD data through the Matz Lab pipeline. First, it uses the `matz2Phylip()` function to convert allsites and varsites files to Phylip format, and also contains the `calc_basic_stats()` function which takes in the retab file and varsites vcf file and summarizes basic statistics, including the number of sites, loci, SNPs, SNPs/locus, and SNPs/site.
