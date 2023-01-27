library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggtree) # had to download directly from Github using devtools: YuLab-SMU/ggtree
library(ape)
library(phytools)
library(ggpmisc)
library(here)

theme_set(theme_cowplot())

## The following code generates Figs. 7 and S1, which illustrate the number and proportion
## of unambiguous changes ("shared sites") along the trees for each dataset and sampling depth.
## These values were calculated using PAUP.

##    FILES REQUIRED:
##          data_files_input_into_scripts/unambig_sums.txt
##          data_files_input_into_scripts/master.nexus
##          data_files_input_into_scripts/Node_numbering_master_trees.png for standardizing node numbering

##    STRUCTURE OF CODE:
##          The figure, being slightly complex, requires multiple components be built separately
##          and then combined into a single figure. The figure is built as follows, and the script
##          follows this numbering:
##              (1) Base cladograms built, one for each Rana and Epipedobates
##              (2) Place small inset plots along branches of cladograms
##              (3) Place root edge inset plot + blown up plot for y axis scale onto plot


# Load and tidy data ------------------------------------------------------

# These data were taken from parsed files section B #3
unambig_sums <- read_tsv(here("data_files_input_into_scripts", "unambig_sums.txt"))

unambig_sums <-
  unambig_sums %>%
  separate(edge, c("edge_fr", "edge_to"), sep="->")

## Remove branches leading to terminal tips; don't do this if you want histograms
## along the terminal branches
# sums_notips <-
#   unambig_sums %>%
#   filter(!str_detect(edge_to, 'e|a|R|E'))

## If you want props relative to ONLY internal branches:
#   sums_notips %>% 
#   group_by(dataset, samp_depth) %>% 
#   mutate(total_sums=sum(sum_unambig)) %>% 
#   mutate(prop_unambig_change_int=sum_unambig/total_sums)

###### CHANGE THIS DEPENDING ON THE BRANCHES YOU WANT PLOTTED ###############

data <-
  unambig_sums # change to sums_notips if you only want internal branches plotted


# Import trees ------------------------------------------------------------

## Import tree we'll use to plot with and check node numbering
trees <- read.nexus(here("data_files_input_into_scripts", "master.nexus"))
rana <- trees[1]
epi <- trees[2]

## Plot the tree using ggtree to see the node numbering
## We'll worry about flipping the nodes later
# epi_wlabs <- ggtree(epi) + geom_text2(aes(subset=!isTip, label=node)) + geom_tiplab()
# rana_wlabs <- ggtree(rana) + geom_text2(aes(subset=!isTip, label=node)) + geom_tiplab()

# Standardize node numbering ----------------------------------------------

## PAUP trees are inputted identically, and so node numbering is standardized
## across 2bRAD and ddRAD (see master_trees.png for numbering details)
## But, this node numbering isn't consistent with the base tree plots
## So we're going to add an additional 'node' column to plot onto tree figs
## Also, we'll remove two branches from root that have 0 values

## EPI
epi_data <-
  data %>%
  filter(dataset=="epi2brad" | dataset=="epiddrad") %>% 
  mutate(node = case_when(
    edge_to == 11 ~ 15,
    edge_to == 12 ~ 16,
    edge_to == 13 ~ 14,
    edge_to == 14 ~ 17,
    edge_to == 15 ~ 13,
    edge_to == 16 ~ 18,
    edge_to == 17 ~ 11,
    edge_to == "Eant_T6859a" ~ 2,
    edge_to == "Eant_T6857" ~ 1,
    edge_to == "Etri_T6836" ~ 4,
    edge_to == "Etri_T6842" ~ 3,
    edge_to == "Ebou_R0153" ~ 6,
    edge_to == "Ebou_R0156" ~ 5,
    edge_to == "Snub_R0159" ~ 8,
    edge_to == "Snub_R0158" ~ 7,
    edge_to == "Ahah_R0090" ~ 10,
    edge_to == "Ahah_R0089a" ~ 9
  ))

## RANA
rana_data <-
  data %>%
  filter(dataset=="rana2brad" | dataset=="ranaddrad") %>%
  mutate(node = case_when(
    edge_to == 18 ~ 14,
    edge_to == 16 ~ 15,
    edge_to == 17 ~ 13,
    edge_to == 11 ~ 17,
    edge_to == 12 ~ 18,
    edge_to == 13 ~ 16,
    edge_to == 15 ~ 11,
    edge_to == "Rber_T1114" ~ 2,
    edge_to == "Rber_T1113a" ~ 1,
    edge_to == "Rneo_T480" ~ 4,
    edge_to == "Rneo_T527" ~ 3,
    edge_to == "Rbla_D2864" ~ 6,
    edge_to == "Rbla_D2865" ~ 5,
    edge_to == "Rsph_T25870" ~ 8,
    edge_to == "Rsph_T26064" ~ 7,
    edge_to == "Rchi_T2034a" ~ 10,
    edge_to == "Rchi_T2049" ~ 9
  ))


# Get data into correct format --------------------------------------------

nodebar_format <- function (data, yaxis) {
  yaxis = enquo(yaxis) # nonstandard encoding
  data %>% 
    filter(node!=11) %>% 
    dplyr::select(samp_depth, !!yaxis, node) %>% 
    spread(samp_depth, !!yaxis)
    # %>% filter(t1!=0 & t2!=0 & t3!=0 & total!=0)
}

## EPI
# We'll plot the root branch histogram separately later
epi2b_data_sums <- 
  epi_data %>% filter(dataset=="epi2brad") %>% nodebar_format(yaxis = sum_unambig)
epidd_data_sums <- 
  epi_data %>% filter(dataset=="epiddrad") %>% nodebar_format(yaxis = sum_unambig)

epi2b_data_prop <- 
  epi_data %>% filter(dataset=="epi2brad") %>% nodebar_format(yaxis = prop_unambig_change)
epidd_data_prop <- 
  epi_data %>% filter(dataset=="epiddrad") %>% nodebar_format(yaxis = prop_unambig_change)

## RANA
rana2b_data_sums <- 
  rana_data %>% filter(dataset=="rana2brad") %>% nodebar_format(yaxis = sum_unambig)
ranadd_data_sums <- 
  rana_data %>% filter(dataset=="ranaddrad") %>% nodebar_format(yaxis = sum_unambig)

rana2b_data_prop <- 
  rana_data %>% filter(dataset=="rana2brad") %>% nodebar_format(yaxis = prop_unambig_change)
ranadd_data_prop <- 
  rana_data %>% filter(dataset=="ranaddrad") %>% nodebar_format(yaxis = prop_unambig_change)


# (1) Build the "base" trees ----------------------------------------------------------

rana_labs <- c("berlandieri", "berlandieri", "neovolcanica", "neovolcanica",
               "blairi", "blairi", "sphenocephala", "sphenocephala",
               "chiricahuensis", "chiricahuensis")

epi_labs <- c("anthonyi", "anthonyi",
              "tricolor", "tricolor",
              "boulengeri", "boulengeri",
              "nubicola", "nubicola",
              "hahneli", "hahneli")

epi_genera <- c("E.", "E.", "E.", "E.", "E.", "E.",
                "S.", "S.", "A.", "A.")

rana_genera <- c("R.", "R.", "R.", "R.", "R.", 
                 "R.", "R.", "R.", "R.", "R.")

e_numbers <- c(1, 2, 2, 1, 2, 1, 1, 2, 1, 2)
r_numbers <- c(1, 2, 2, 1, 2, 1, 2, 1, 2, 1)

epilabs <- data.frame(label = epi$epi$tip.label, 
                      species = epi_labs, 
                      sample = e_numbers,
                      genera = epi_genera)

ranalabs <- data.frame(label = rana$rana$tip.label, 
                       species = rana_labs, 
                       sample = r_numbers,
                       genera = rana_genera)

## Build trees upon which we'll put histograms
## In order for ggtree to match phylo object to dataframe, the FIRST COLUMN must match the tip labels
## of the phylo object!!
p.epi <-
  ggtree(epi, size=1.5) %<+% 
  epilabs +
  coord_cartesian(xlim = c(0,5.5), clip="off") +
  geom_tiplab(aes(label=paste0('italic(', genera, ')~italic(', species, ')~', e_numbers)), 
  parse=T, size=5) +
  theme(plot.margin=unit(c(2,3,.75,1.5),"cm")) # top, right, bottom, left

p.epi <-
  p.epi %>% 
  flip(15, 16)

# ====== RANA =======
p.rana <-
  ggtree(rana, size=1.5) %<+% 
  ranalabs +
  coord_cartesian(xlim = c(0,5.5), clip="off") +
  geom_tiplab(aes(label=paste0('italic(', genera, ')~italic(', species, ')~', r_numbers)), 
              parse=T, size=5) +
  theme(plot.margin=unit(c(2,2,.75,5),"cm"))
p.rana <-
  p.rana %>% 
  flip(13, 16) %>% 
  flip(14, 15) %>% 
  flip(17, 18)


# Edit ggtree nodebar function --------------------------------------------

nodebar_eac <- function (data, cols, alpha = 1, position="dodge", ylim) {
  if (!"node" %in% colnames(data)) {
    stop("data should have a column 'node'...")
  }
  type <- value <- NULL # creates empty values for type and value, to be filled in next step
  ldf <- gather(data, type, value, !!cols) %>% split(., .$node) # creates sep tibbles for each node w type + value cols
  bars <- lapply(ldf, function(df)
    ggplot(df, aes(x = 1, y = value, fill = type)) +
      geom_bar(stat = "identity", alpha = alpha, position = 'dodge') +
      theme(strip.background = element_blank(),
            aspect.ratio = 1,
            axis.ticks.x = element_blank(),
            axis.title = element_blank(),
            axis.ticks.y = element_line(color="white", size=0.75),
            axis.line.y = element_line(color="white", size=0.75),
            axis.line.x = element_line(color="white", size=0.75),
            axis.ticks.length = unit(0.2,"cm"),
            legend.position='none',
            plot.background = element_rect(fill = "#e6e6e6"),
            plot.margin = margin(1,1,1,1,"mm"),
            axis.text = element_blank()) +
      scale_x_discrete(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0), limits=ylim) +
      scale_fill_manual(values = c(t1="#a1dab4", t2="#41b6c4", t3="#2c7fb8", total="#253494"))
  )
  return(bars)
}

# Build dataframe for small insets -----------------------------------------------------------

## UNAMBIG SUMS - NO TIPS
epi2b_bars_sums <- nodebar_eac(epi2b_data_sums, cols=2:5, ylim=c(0,900)) #ymax=888
epidd_bars_sums <- nodebar_eac(epidd_data_sums, cols=2:5, ylim=c(0,7100)) #ymax=7093

rana2b_bars_sums <- nodebar_eac(rana2b_data_sums, cols=2:5, ylim=c(0,2550)) #ymax=2506
ranadd_bars_sums <- nodebar_eac(ranadd_data_sums, cols=2:5, ylim=c(0,17000)) #ymax=16673

## INTERNAL BRANCHES ONLY - NO TIPS
# epi2b_bars_prop <- nodebar_eac(epi2b_data_prop, cols=3:6, ylim=c(0,0.41)) # for int ymax=0.303
# epidd_bars_prop <- nodebar_eac(epidd_data_prop, cols=3:6, ylim=c(0,0.41)) # for int ymax=0.401
# rana2b_bars_prop <- nodebar_eac(rana2b_data_prop, cols=3:6, ylim=c(0,0.55)) # for int ymax=0.549
# ranadd_bars_prop <- nodebar_eac(ranadd_data_prop, cols=3:6, ylim=c(0,0.32)) # for int ymax=0.316

## ALL BRANCHES PLOTTED
epi2b_bars_prop <- nodebar_eac(epi2b_data_prop, cols=2:5, ylim=c(0,0.23)) #ymax=0.189
epidd_bars_prop <- nodebar_eac(epidd_data_prop, cols=2:5, ylim=c(0,0.23)) #ymax=0.225

rana2b_bars_prop <- nodebar_eac(rana2b_data_prop, cols=2:5, ylim=c(0,0.31)) #ymax=0.303
ranadd_bars_prop <- nodebar_eac(ranadd_data_prop, cols=2:5, ylim=c(0,0.25)) #ymax=0.241


# (2) Place inset plots on tree -------------------------------------------------------

small_insets <- function (tree, bars) {
  inset(tree, bars, width=0.3, height=0.1, vjust=-.475, hjust=.5)
}

# Sums
epi2b_plot <- small_insets(tree = p.epi, bars = epi2b_bars_sums)
epidd_plot <- small_insets(tree = p.epi, bars = epidd_bars_sums)

rana2b_plot <- small_insets(tree = p.rana, bars = rana2b_bars_sums)
ranadd_plot <- small_insets(tree = p.rana, bars = ranadd_bars_sums)

# Proportions
epi2b_plot_prop <- small_insets(tree = p.epi, bars = epi2b_bars_prop)
epidd_plot_prop <- small_insets(tree = p.epi, bars = epidd_bars_prop)

rana2b_plot_prop <- small_insets(tree = p.rana, bars = rana2b_bars_prop)
ranadd_plot_prop <- small_insets(tree = p.rana, bars = ranadd_bars_prop)


# Build blown up single plot for yaxis scale -------------------------------

blowups <- function (data, yaxis, ylim, ytitle) {
  data %>% 
    ggplot(aes(x=samp_depth, y=yaxis, fill=samp_depth)) +
    geom_bar(stat="identity") +
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
    scale_fill_manual(values=c("#a1dab4", "#41b6c4", "#2c7fb8", "#253494")) +
    ylab(ytitle) +
    xlab("Sampling depth")
}


# For Epi, this is node 14
epi2b_inset <-
  epi_data %>%
  filter(dataset=="epi2brad" & node==14) %>%
  blowups(yaxis=.$sum_unambig, ylim = c(0,900), ytitle = "No. total changes")

epidd_inset <-
  epi_data %>%
  filter(dataset=="epiddrad" & node==14) %>%
  blowups(yaxis=.$sum_unambig, ylim=c(0,7100), ytitle = "No. total changes")

epi2b_inset_prop <-
  epi_data %>% 
  filter(dataset=="epi2brad" & node==14) %>% 
  blowups(yaxis=.$prop_unambig_change, ylim=c(0,0.23), ytitle = "Prop. total changes") # all=0.23; int=0.31

epidd_inset_prop <-
  epi_data %>% 
  filter(dataset=="epiddrad" & node==14) %>% 
  blowups(yaxis=.$prop_unambig_change, ylim=c(0,0.23), ytitle = "Prop. total changes") # all=0.23; int=0.41

########################
# For Rana, node 13
rana2b_inset <-
  rana_data %>%
  filter(dataset=="rana2brad" & node==13) %>%
  blowups(yaxis=.$sum_unambig, ylim=c(0,2550), ytitle = "No. total changes")

ranadd_inset <-
  rana_data %>%
  filter(dataset=="ranaddrad" & node==13) %>%
  blowups(yaxis=.$sum_unambig, ylim=c(0,17000), ytitle = "No. total changes")

rana2b_inset_prop <-
  rana_data %>% 
  filter(dataset=="rana2brad" & node==13) %>% 
  blowups(yaxis=.$prop_unambig_change, ylim=c(0,0.31), ytitle = "Prop. shared sites") # all=0.31; int=0.55

ranadd_inset_prop <-
  rana_data %>% 
  filter(dataset=="ranaddrad" & node==13) %>% 
  blowups(yaxis=.$prop_unambig_change, ylim=c(0,0.25), ytitle = "Prop. shared sites") # all=0.25; int=0.32


# Build plot for root edges -----------------------------------------------

# Build root plots
root_plots <- function (data, yaxis, ylim) {
  ggplot(data, aes(x=samp_depth, y=yaxis, fill=samp_depth)) +
    geom_bar(stat="identity") +
    theme(strip.background = element_blank(),
          aspect.ratio = 1,
          axis.ticks.x = element_blank(),
          axis.title = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_line(color="white", size=1),
          axis.line.y = element_line(color="white", size=1),
          axis.line.x = element_line(color="white", size=1),
          axis.ticks.length = unit(0.5,"cm"),
          legend.position='none',
          plot.background = element_rect(fill = "#e6e6e6")) +
    scale_y_continuous(expand = c(0,0), limits=ylim) +
    scale_fill_manual(values=c("#a1dab4", "#41b6c4", "#2c7fb8", "#253494"))
}

## SUMS
epi2b_root_p <-
  epi_data %>%
  filter(dataset=="epi2brad" & node==11) %>%
  root_plots(yaxis=.$sum_unambig, ylim=c(0,900))

epidd_root_p <-
  epi_data %>%
  filter(dataset=="epiddrad" & node==11) %>%
  root_plots(yaxis=.$sum_unambig, ylim=c(0,7100))

rana2b_root_p <-
  rana_data %>%
  filter(dataset=="rana2brad" & node==11) %>%
  root_plots(yaxis=.$sum_unambig, ylim=c(0,2550))

ranadd_root_p <-
  rana_data %>%
  filter(dataset=="ranaddrad" & node==11) %>%
  root_plots(yaxis=.$sum_unambig, ylim=c(0,24000)) # separate axis for this plot

## PROPORTIONS
epi2b_root_prop_p <-
  epi_data %>%
  filter(dataset=="epi2brad" & node==11) %>%
  root_plots(yaxis=.$prop_unambig_change, ylim=c(0,0.23)) # all=0.23; int=0.31

epidd_root_prop_p <-
  epi_data %>%
  filter(dataset=="epiddrad" & node==11) %>%
  root_plots(yaxis=.$prop_unambig_change, ylim=c(0,0.23)) # all=0.23; int=0.41

rana2b_root_prop_p <-
  rana_data %>%
  filter(dataset=="rana2brad" & node==11) %>%
  root_plots(yaxis=.$prop_unambig_change, ylim=c(0,0.31)) # all=0.31; int=0.55

ranadd_root_prop_p <-
  rana_data %>%
  filter(dataset=="ranaddrad" & node==11) %>%
  root_plots(yaxis=.$prop_unambig_change, ylim=c(0,0.25)) # all=0.25; int=0.32


# (3) Place root edge plot + blow up plot to get final figure ------------------------------------------

epi_final_plot <- function (base_plot, blowup, root, title) {
  base_plot +
    draw_plot(blowup, x=-1, y=8.5, width=3, height=3) +
    draw_plot(root, x=-.9, y=2.5, width=1.75, height=1.75) +
    geom_segment(aes(x=2.15, y=8.525, xend=1.575, yend=8.525), size=1, color="white") +
    geom_segment(aes(x=2.15, y=9.4, xend=1.8, yend=11.49), size=1, color="white") +
    theme(plot.title = element_text(face="bold", hjust = 0.5, size=16)) +
    ggtitle(title)
}

## Sums
final_epi2b_sums <-
  epi_final_plot(base_plot = epi2b_plot, blowup = epi2b_inset, root = epi2b_root_p, title = "2bRAD")
final_epidd_sums <-
  epi_final_plot(base_plot = epidd_plot, blowup = epidd_inset, root = epidd_root_p, title = "ddRAD")

## Proportions
final_epi2b_prop <-
  epi_final_plot(base_plot = epi2b_plot_prop, blowup = epi2b_inset_prop, root = epi2b_root_prop_p, title = "2bRAD")
final_epidd_prop <-
  epi_final_plot(base_plot = epidd_plot_prop, blowup = epidd_inset_prop, root = epidd_root_prop_p, title = "ddRAD")

### RANA (dimensions differ slightly)
rana_final_plot <- function (base_plot, blowup, root, title) {
  base_plot +
    draw_plot(blowup, x=-2.25, y=8.5, width=3, height=3) +
    draw_plot(root, x=-.9, y=3, width=1.75, height=1.75) +
    geom_segment(aes(x=1.1, y=8.525, xend=.57, yend=8.525), size=1, color="white") +
    geom_segment(aes(x=1.1, y=9.4, xend=.75, yend=11.49), size=1, color="white") +
    theme(plot.title = element_text(face="bold", hjust = 0.5, size=16)) +
    ggtitle(title)
}

## Sums
final_rana2b_sums <-
  rana_final_plot(base_plot = rana2b_plot, blowup = rana2b_inset, root = rana2b_root_p, title = "2bRAD")
final_ranadd_sums <-
  rana_final_plot(base_plot = ranadd_plot, blowup = ranadd_inset, root = ranadd_root_p, title = "ddRAD")

## Proportions
final_rana2b_prop <-
  rana_final_plot(base_plot = rana2b_plot_prop, blowup = rana2b_inset_prop, root = rana2b_root_prop_p, title = "2bRAD")
final_ranadd_prop <-
  rana_final_plot(base_plot = ranadd_plot_prop, blowup = ranadd_inset_prop, root = ranadd_root_prop_p, title = "ddRAD")

# Save plots ------------------------------------------------------------

## Sums
final_epi2b_sums
ggsave("epi2b_unambig_sums.pdf", width=8.36, height = 11.89)
final_epidd_sums
ggsave("epidd_unambig_sums.pdf", width=8.36, height = 11.89)

final_rana2b_sums
ggsave("rana2b_unambig_sums.pdf", width=8.9, height = 11.89)
final_ranadd_sums
ggsave("ranadd_unambig_sums.pdf", width=8.9, height = 11.89)

## Proportions
final_epi2b_prop
ggsave("epi2b_unambig_props.pdf", width=8.36, height = 11.89)
final_epidd_prop
ggsave("epidd_unambig_props.pdf", width=8.36, height = 11.89)

final_rana2b_prop
ggsave("rana2b_unambig_props.pdf", width=8.9, height = 11.89)
final_ranadd_prop
ggsave("ranadd_unambig_props.pdf", width=8.9, height = 11.89)
