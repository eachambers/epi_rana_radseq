library(tidyverse)
library(cowplot)
library(ggplot2)
library(scales)
library(rlang)
library(here)

theme_set(theme_cowplot())

## This code does the following:
##     Builds figures to determine appropriate clustering threshold value 
##     as it relates to consensus reads and loci in assembly (Fig. 3). The input
##     file is generated using the clust_threshold_processing.R script.

##    FILES REQUIRED:
##          data_files_input_into_scripts/clust_threshold_data.txt


# Import and process data -------------------------------------------------

dat <- read_tsv(here("data_files_input_into_scripts", "clust_threshold_data.txt"), col_names = TRUE)

# Ensure that things are ordered correctly
dat$sample <- factor(dat$sample, levels=c("Rber_T1113a", "Rber_T1113b", "Rber_T1114", "Rneo_T480", "Rneo_T527", "Rsph_T25870", "Rsph_T26064",
                                          "Rbla_D2864", "Rbla_D2865", "Rchi_T2034a", "Rchi_T2034b", "Rchi_T2049", "Eant_T6859a", "Eant_T6859b", "Eant_T6857", "Etri_T6836", "Etri_T6842", "Ebou_R0153", "Ebou_R0156",
                                          "Snub_R0158", "Snub_R0159", "Ahah_R0089a", "Ahah_R0089b", "Ahah_R0090"))

# BUILDING FIGURES -----------------------------------------------

# Custom color palette (partially based off wesanderson package)
palette <- c("#8D529E", "#BA85A2", "#6B2391", "#CA3518", "#C66E1F", 
             "#BA8517", "#D9C339", "#086955", "#7C935E", "#80A7B3", "#386784", "#008CA2",
             "#8D529E", "#BA85A2", "#6B2391", "#CA3518", "#C66E1F", 
             "#BA8517", "#D9C339", "#086955", "#7C935E", "#80A7B3", "#386784", "#008CA2")


# Consensus reads ---------------------------------------------------------

## To combine these all into a single plot (with loci_assembly), make each figure separately
# facet on 2brad or ddrad

#' Build two panels (same y axis) of Figure 3, illustrating clustering threshold
#' vs either consensus reads or loci in assembly
#'
#' @param data dataframe (dat from above)
#' @param taxon taxon to be plotted, either "epi" or "rana"
#' @param yaxisval metric on the y axis, either "reads_consens" or "loci_assembly" (must match col name)
#' @param ylabel label of y axis ("Consensus reads" or "Loci in assembly")
#'
#' @return
#' @export
#'
#' @examples
clust_fig <- function(data, taxon, yaxisval, ylabel) {
  data %>% 
    filter(taxon == taxon) %>% 
    ggplot(aes(x=clust_threshold, {{yaxisval}}, group = sample, color = sample)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_vline(xintercept = 91, linetype = "dashed", linewidth = 1) +
    scale_y_continuous(labels = scales::comma) +
    theme(legend.position = "none",
          axis.text = element_text(size = 14),
          axis.title = element_text(size = 17),
          panel.grid.major.x = element_line(linewidth = .75, color = "snow2"),
          panel.grid.minor.x = element_line(linewidth = .75, color = "snow2"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold")) +
    scale_x_continuous(minor_breaks = seq(80,95,1)) +
    scale_color_manual(values = palette) +
    xlab("Clustering threshold") +
    ylab(ylabel) +
    facet_wrap(~method, scales = "free", nrow = 2)
}


# Build plots -------------------------------------------------------------

epi_consens <-
  clust_fig(data = dat, taxon = "epi", yaxisval = reads_consens, ylabel = "Consensus reads")
epi_loci <-
  clust_fig(data = dat, taxon = "epi", yaxisval = loci_assembly, ylabel = "Loci in assembly")

rana_consens <-
  clust_fig(data = dat, taxon = "rana", yaxisval = reads_consens, ylabel = "Consensus reads")
rana_loci <-
  clust_fig(data = dat, taxon = "rana", yaxisval = loci_assembly, ylabel = "Loci in assembly")


# Add in legends ----------------------------------------------------------

# Pull the legends so that you only have one per taxon in finished plot
legendepi <- get_legend(p.epiloci)
legendrana <- get_legend(p.ranaloci)

# Now plot, adding in legend afterward
p_epi <- plot_grid(epi_consens, epi_loci + theme(legend.position = "none"))
p_epi <- plot_grid(p_epi, legendepi, rel_widths = c(1,.15))

p_rana <- plot_grid(rana_consens, rana_loci + theme(legend.position = "none"))
p_rana <- plot_grid(p_rana, legendrana, rel_widths = c(1,.15))

# Build overall figure
p_clust_threshold <- plot_grid(p_rana, p_epi, nrow = 2)


# Export plot -------------------------------------------------------------

# Combine everything together and export as PDF
save_plot("Fig3_clustthreshold.pdf", p_clust_threshold, base_aspect_ratio = 1.6, base_height = 11)
