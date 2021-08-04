setwd("~/Box Sync/Epipedobates project/SuppMaterials/5_Data_visualization")

library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggtree) # had to download directly from Github using devtools: YuLab-SMU/ggtree
library(ape)
library(phytools)
library(ggpmisc)

theme_set(theme_cowplot())

## The following code generates Fig. 5, showing the number of changes that differ significantly & non-sig.
## from random expectations (another way of putting it: detecting phylogenetic signal of state changes along
## the trees). SNP data have been re-coded to 0s and 1s. Values were calculated using PAUP.

##    FILES REQUIRED:
##          data_files_input_into_scripts/recoded_signonsig.txt
##          Use Fig4_PAUP_analysis.R to generate base tree figures (p.epi and p.rana)

##    STRUCTURE OF CODE:
##              (1) Base cladograms built [using Fig4_PAUP_analysis.R], one for each Rana and Epipedobates
##              (2) Place small inset plots along branches of cladograms
##              (3) Place root edge inset plot + blown up plot for y axis scale onto plot


# Import and tidy data ----------------------------------------------------

sig <- read_tsv("data_files_input_into_scripts/recoded_signonsig.txt")

# Separate based on nodes
# Add col for proportion data
prop_sig <-
  sig %>% 
  # replace(is.na(.), 0) %>%
  separate(EdgeName, c("edge_fr", "edge_to"), sep=" --> ") %>% 
  group_by(DataSetName) %>% 
  mutate(sums = sum(n)) %>% 
  mutate(prop = n/sums)

prop_sig %>% 
  group_by(DataSetName, QualityScore) %>% 
  summarize(sum(n))


# Standardize node numbering ----------------------------------------------

epi_data <-
  prop_sig %>%
  filter(Taxon=="epi") %>% 
  mutate(node = case_when(
    edge_to == 11 ~ 15,
    edge_to == 12 ~ 16,
    edge_to == 13 ~ 14,
    edge_to == 14 ~ 17,
    edge_to == 15 ~ 13,
    edge_to == 16 ~ 18,
    edge_to == 18 ~ 11,
    edge_to == "EantT6859a" ~ 2,
    edge_to == "EantT6857" ~ 1,
    edge_to == "EtriT6836" ~ 4,
    edge_to == "EtriT6842" ~ 3,
    edge_to == "EbouR0153" ~ 6,
    edge_to == "EbouR0156" ~ 5,
    edge_to == "SnubR0159" ~ 8,
    edge_to == "SnubR0158" ~ 7,
    edge_to == "AhahR0090" ~ 10,
    edge_to == "AhahR0089a" ~ 9 
  ))

rana_data <-
  prop_sig %>%
  filter(Taxon=="rana") %>% 
  mutate(node = case_when(
    edge_to == 11 ~ 17,
    edge_to == 12 ~ 18,
    edge_to == 13 ~ 16,
    edge_to == 14 ~ 11,
    edge_to == 16 ~ 15,
    edge_to == 17 ~ 13,
    edge_to == 18 ~ 14,
    edge_to == "RberT1114" ~ 2,
    edge_to == "RberT1113a" ~ 1,
    edge_to == "RneoT480" ~ 4,
    edge_to == "RneoT527" ~ 3,
    edge_to == "RblaD2864" ~ 6,
    edge_to == "RblaD2865" ~ 5,
    edge_to == "RsphT25870" ~ 8,
    edge_to == "RsphT26064" ~ 7,
    edge_to == "RchiT2034a" ~ 10,
    edge_to == "RchiT2049" ~ 9
  ))

# Change levels so data are plotted correctly
epi_data$QualityScore = factor(epi_data$QualityScore, levels=c("Sign", "NS"))
rana_data$QualityScore = factor(rana_data$QualityScore, levels=c("Sign", "NS"))

# Get the data into the correct format for the nodebars function

nodebar_format <- function (data, dataset, col) {
  dataset = enquo(dataset)
  col = enquo(col)
  data %>% 
    filter(DataSetName==!!dataset) %>% 
    dplyr::select(!!col, node, StateChange, QualityScore) %>% 
    spread(QualityScore, !!col)
}

## EPI
epi2b_data_sig <- epi_data %>% nodebar_format(dataset = "epi_2brad_total", col = prop)
epidd_data_sig <- epi_data %>% nodebar_format(dataset = "epi_ddrad_total", col = prop)
## RANA
rana2b_data_sig <- rana_data %>% nodebar_format(dataset = "rana_2brad_total", col = prop)
ranadd_data_sig <- rana_data %>% nodebar_format(dataset = "rana_ddrad_total", col = prop)

# (1) Build the base trees ---------------------------------------------------------

# Use PAUP_analysis.R to build p.epi and p.rana

# Adapted nodebar function ---------------------------------------------------

## This is simply altering the existing nodebar function built within ggtree with 
## a custom aesthetic

nodebar_eac <- function (data, cols, ylim, position="stack") {
  if (!"node" %in% colnames(data)) {
    stop("data should have a column 'node'...")
  }
  type <- value <- NULL # creates empty values for type and value, to be filled in next step
  ldf <- gather(data, type, value, !!cols)
  ldf$type = factor(ldf$type, levels=c("Sign", "NS"))
  ldf <- ldf %>% split(., .$node) # creates sep tibbles for each node w type + value cols
  bars <- lapply(ldf, function(df)
    ggplot(df, aes(x = type, y = value, fill = StateChange)) +
      geom_bar(aes(alpha = type), stat = "identity", position = 'stack') +
      theme(strip.background = element_blank(),
            aspect.ratio = 1,
            axis.ticks.x = element_blank(),
            axis.title = element_blank(),
            axis.ticks.y = element_line(color="white", size=0.75),
            axis.line.y = element_line(color="white", size=0.75),
            axis.line.x = element_line(color="white", size=0.75),
            axis.ticks.length = unit(0.2,"cm"),
            legend.position='none',
            axis.text = element_blank(),
            plot.background = element_rect(fill = "#e6e6e6"),
            plot.margin = margin(1,1,-1.5,1,"mm")) + # t, r, b, l
      scale_x_discrete(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0), limits=ylim) +
      scale_fill_manual(values = c("0 <=> 1"="#91d5de", "1 ==> 0"="#eaae37", "1 <=> 0"="#b4674e", "0 ==> 1"="#2e8289")) +
      scale_alpha_manual(values = c(1, 0.4))
  )
  return(bars)
}

# Run function to build plots ---------------------------------------------

## Don't build plot for root branch (node 11); this needs to be enlarged and centred

epi2b_bars_sig <-
  epi2b_data_sig %>% 
  filter(node!=11) %>% 
  nodebar_eac(cols=4:5, ylim=c(0,0.25))

epidd_bars_sig <-
  epidd_data_sig %>% 
  filter(node!=11) %>% 
  nodebar_eac(cols=4:5, ylim=c(0,0.25))

rana2b_bars_sig <-
  rana2b_data_sig %>% 
  filter(node!=11) %>% 
  nodebar_eac(cols=4:5, ylim=c(0,0.15))

ranadd_bars_sig <-
  ranadd_data_sig %>% 
  filter(node!=11) %>% 
  nodebar_eac(cols=4:5, ylim=c(0,0.15))

# (2) Place small plots on trees -------------------------------------------------------

epi2b_plot_sig <-
  inset(p.epi, epi2b_bars_sig, width=0.3, height=0.1, vjust=-.475, hjust=.5)

epidd_plot_sig <-
  inset(p.epi, epidd_bars_sig, width=0.3, height=0.1, vjust=-.475, hjust=.5)

rana2b_plot_sig <-
  inset(p.rana, rana2b_bars_sig, width=0.3, height=0.1, vjust=-.475, hjust=.5)

ranadd_plot_sig <-
  inset(p.rana, ranadd_bars_sig, width=0.3, height=0.1, vjust=-.475, hjust=.5)

# Build blown up single plot for yaxis scale -------------------------------

sig_blownup <- function (data, yaxis, ytitle, ylim) {
  data %>% 
    ggplot(aes(x=QualityScore, y=yaxis, fill=StateChange)) +
    geom_bar(aes(alpha = QualityScore), stat="identity") +
    theme(strip.background = element_blank(),
          aspect.ratio = 1,
          axis.ticks.x = element_blank(),
          axis.title = element_text(size=15, face = "bold"),
          axis.text.x = element_text(size=14, margin = margin(t=-7)),
          axis.text.y = element_text(size=14),
          axis.ticks.y = element_line(color="white", size=1),
          axis.line.y = element_line(color="white", size=1),
          axis.line.x = element_line(color="white", size=1),
          axis.ticks.length = unit(0.5,"cm"),
          legend.position='none',
          plot.background = element_rect(fill = "#e6e6e6", colour="gray48", size=1)) +
    scale_y_continuous(expand = c(0,0), limits=ylim) +
    scale_fill_manual(values=c("0 <=> 1"="#91d5de", "1 ==> 0"="#eaae37", "1 <=> 0"="#b4674e", "0 ==> 1"="#2e8289")) +
    ylab(ytitle) +
    xlab("Significance") +
    scale_alpha_manual(values = c(1, 0.4))
}

##### EPI: node 14 #####
epi2b_inset_sig <-
  epi_data %>% 
  filter(DataSetName == "epi_2brad_total" & node == 14) %>% 
  sig_blownup(yaxis = .$prop, ytitle = "Prop. total changes", ylim=c(0,0.25))

epidd_inset_sig <-
epi_data %>% 
  filter(DataSetName == "epi_ddrad_total" & node == 14) %>% 
  sig_blownup(yaxis = .$prop, ytitle = "Prop. total changes", ylim=c(0,0.25))

##### RANA: node 13 #####
rana2b_inset_sig <-
  rana_data %>% 
  filter(DataSetName == "rana_2brad_total" & node == 13) %>% 
  sig_blownup(yaxis = .$prop, ytitle = "Prop. total changes", ylim=c(0,0.15))

ranadd_inset_sig <-
  rana_data %>% 
  filter(DataSetName == "rana_ddrad_total" & node == 13) %>% 
  sig_blownup(yaxis = .$prop, ytitle = "Prop. total changes", ylim=c(0,0.15))

# Build plot for root edges -----------------------------------------------

# Build root plots; change y to prop for second plots
# These are node 11 for both Epi and Rana

root_plots_sig <- function (data, yaxis, ylim) {
  yaxis = enquo(yaxis)
  ggplot(data, aes(x=QualityScore, y=!!yaxis, fill=StateChange)) +
    geom_bar(aes(alpha = QualityScore), stat="identity", position="stack") +
    theme(strip.background = element_blank(),
          aspect.ratio = 1,
          axis.ticks.x = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(), # size 12
          axis.ticks.y = element_line(color="white", size=1),
          axis.line.y = element_line(color="white", size=1),
          axis.line.x = element_line(color="white", size=1),
          axis.ticks.length = unit(0.5,"cm"),
          legend.position='none',
          plot.background = element_rect(fill = "#e6e6e6")) +
    scale_y_continuous(expand = c(0,0), limits = ylim) +
    scale_fill_manual(values = c("0 <=> 1"="#91d5de", "1 ==> 0"="#eaae37", "1 <=> 0"="#b4674e", "0 ==> 1"="#2e8289")) +
    scale_alpha_manual(values = c(1, 0.4)) 
}

### EPI ###
epi2b_root_sig <-
  epi_data %>% 
  filter(DataSetName=="epi_2brad_total" & node==11) %>% 
  root_plots_sig(yaxis = prop, ylim=c(0,0.25))

epidd_root_sig <-
  epi_data %>% 
  filter(DataSetName=="epi_ddrad_total" & node==11) %>% 
  root_plots_sig(yaxis = prop, ylim=c(0,0.25))

### RANA ###
rana2b_root_sig <-
  rana_data %>% 
  filter(DataSetName=="rana_2brad_total" & node==11) %>% 
  root_plots_sig(yaxis = prop, ylim=c(0,0.15))

ranadd_root_sig <-
  rana_data %>% 
  filter(DataSetName=="rana_ddrad_total" & node==11) %>% 
  root_plots_sig(yaxis = prop, ylim=c(0,0.15))

# (3) Place blow-up and root plots onto main fig --------------------------

epi_final_plot <- function (base_plot, blowup, root, title) {
  base_plot +
    draw_plot(blowup, x=-1, y=8.5, width=3, height=3) +
    draw_plot(root, x=-.9, y=2.5, width=2, height=2) +
    geom_segment(aes(x=2.15, y=8.525, xend=1.575, yend=8.525), size=1, color="gray48") +
    geom_segment(aes(x=2.15, y=9.4, xend=1.8, yend=11.49), size=1, color="gray48") +
    theme(plot.title = element_text(face="bold", hjust = 0.5, size=16)) +
    ggtitle(title)
}

final_epi2b_sig <-
  epi_final_plot(base_plot = epi2b_plot_sig, blowup = epi2b_inset_sig, root = epi2b_root_sig, title = "2bRAD")
final_epidd_sig <-
  epi_final_plot(base_plot = epidd_plot_sig, blowup = epidd_inset_sig, root = epidd_root_sig, title = "ddRAD")

### RANA (dimensions differ slightly)
rana_final_plot <- function (base_plot, blowup, root, title) {
  base_plot +
    draw_plot(blowup, x=-2.25, y=8.5, width=3, height=3) +
    draw_plot(root, x=-.9, y=3, width=1.75, height=1.75) +
    geom_segment(aes(x=1.1, y=8.525, xend=.57, yend=8.525), size=1, color="gray48") +
    geom_segment(aes(x=1.1, y=9.4, xend=.75, yend=11.49), size=1, color="gray48") +
    theme(plot.title = element_text(face="bold", hjust = 0.5, size=16)) +
    ggtitle(title)
}

final_rana2b_sig <-
  rana_final_plot(base_plot = rana2b_plot_sig, blowup = rana2b_inset_sig, root = rana2b_root_sig, title = "2bRAD")
final_ranadd_sig <-
  rana_final_plot(base_plot = ranadd_plot_sig, blowup = ranadd_inset_sig, root = ranadd_root_sig, title = "ddRAD")

# Save four plots ------------------------------------------------------------

final_epi2b_sig
ggsave("epi2b_sig_props.pdf", width=8.36, height = 11.89)
final_epidd_sig
ggsave("epidd_sig_props.pdf", width=8.36, height = 11.89)

final_rana2b_sig
ggsave("rana2b_sig_props.pdf", width=8.9, height = 11.89)
final_ranadd_sig
ggsave("ranadd_sig_props.pdf", width=8.9, height = 11.89)
