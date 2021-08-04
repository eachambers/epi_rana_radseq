# epi_rana_radseq

Scripts for comparing ddRAD and 2bRAD methods in two groups of frogs, from [Chambers et al. (2021)](LINKHERE); raw data files available on Dryad [here](LINKHERE).

## Code

I. Post-processing bioinformatics pipelines
* [Concatenate iPyrad output .stats files](https://github.com/eachambers/epi_rana_radseq/blob/master/scripts_data_processing/extract_data.ipynb)
* [Transpose 2bRAD output data into Phylip format and calculate replicate loci](https://github.com/eachambers/epi_rana_radseq/blob/master/scripts_data_processing/2bRADnative_processdata.R)
* [Calculate the average read depth across samples from 2bRAD data](https://github.com/eachambers/epi_rana_radseq/blob/master/scripts_data_processing/2bRAD_depth_stats.R)

II. Data visualization
* [Figs. 3 & S4: Proportions of missing data](https://github.com/eachambers/epi_rana_radseq/blob/master/scripts_data_visualization/Fig3_PAUP_analysis.R)
* [Fig. 4: Calculating missing data](https://github.com/eachambers/epi_rana_radseq/blob/master/scripts_data_visualization/Fig4_Missing_data.R)
* [Fig. 5: Phylogenetic information along a tree](https://github.com/eachambers/epi_rana_radseq/blob/master/scripts_data_visualization/Fig5_Recoded_significance_analysis.R)
* [Fig. S1: Data characteristization by clustering threshold](https://github.com/eachambers/epi_rana_radseq/blob/master/scripts_data_visualization/FigS1_Data_characterization.R)
* [Fig. S3: Parsimony-informative sites and retention indices](https://github.com/eachambers/epi_rana_radseq/blob/master/scripts_data_visualization/FigS3_Retention_index.R)
* [Fig. S5: Proportions of state changes from recoded datasets](https://github.com/eachambers/epi_rana_radseq/blob/master/scripts_data_visualization/FigS5_Dollo_analysis.R)
* [Fig. S7a: Split iPyrad output .loci files](https://github.com/eachambers/epi_rana_radseq/blob/master/scripts_data_visualization/FigS7_a_Split_loci_files.ipynb)
* [Fig. S7b: Calculate shared loci between replicate samples](https://github.com/eachambers/epi_rana_radseq/blob/master/scripts_data_visualization/FigS7_b_Shared_loci_replicates.ipynb)
* [Fig. S7c: Construct Fig. S7](https://github.com/eachambers/epi_rana_radseq/blob/master/scripts_data_visualization/FigS7_c_Shared_loci_replicates.R)

III. Associated files (for input into above scripts)
* [Amounts and proportions of missing data (Fig. 3)](https://github.com/eachambers/epi_rana_radseq/tree/master/data/missing_data_snps.txt)
* [Master Epipedobates and Rana trees for plotting (Figs. 3, 5, S4, and S5)](https://github.com/eachambers/epi_rana_radseq/tree/master/data/master.nexus)
* [Unambiguous changes (Figs. 3 and S4)](https://github.com/eachambers/epi_rana_radseq/tree/master/data/unambig_sums.txt)
* [Binary-recoded significant unambiguous changes (Fig. 5)](https://github.com/eachambers/epi_rana_radseq/tree/master/data/recoded_signonsig.txt)
* [Retention indices and parsimony-informative sites (Fig. S3)](https://github.com/eachambers/epi_rana_radseq/tree/master/data/Retention_PIs.csv)
* [Standardized node numbering between master tree and PAUP* output (Fig. S5)](https://github.com/eachambers/epi_rana_radseq/tree/master/data/Node_numbering_master_trees.png)
* [Unambiguous changes for binary-recoded datasets (Fig. S5)](https://github.com/eachambers/epi_rana_radseq/tree/master/data/Plot-Data-for-MS-FigS5.txt)
* [Shared loci between 2bRAD replicate samples (Fig. S7)](https://github.com/eachambers/epi_rana_radseq/tree/master/data/2bRAD_shared_loci_replicates.csv)
* [Shared loci between ddRAD replicate samples (Fig. S7)](https://github.com/eachambers/epi_rana_radseq/tree/master/data/ddRAD_shared_loci_replicates.csv)
