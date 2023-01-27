library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggtree) # had to download directly from Github using devtools: YuLab-SMU/ggtree
library(ape)
library(phytools)
library(ggpmisc)
library(here)

theme_set(theme_cowplot())

## The following code generates Fig. 9, which illustrates the number of unambiguous
## changes (categorized by state change) along the trees for each dataset and sampling
## depth, when the datasets have been re-coded to 0s and 1s. Values were calculated using PAUP.

##    FILES REQUIRED:
##          data_files_input_into_scripts/Plot-Data-for-MS-FigS5.txt
##          Use Fig7&S1_PAUP_analysis.R to generate base tree figures (p.epi and p.rana)

##    STRUCTURE OF CODE:
##              (1) Base cladograms built [using Fig7&S1_PAUP_analysis.R], one for each Rana and Epipedobates
##              (2) Place small inset plots along branches of cladograms
##              (3) Place root edge inset plot + blown up plot for y axis scale onto plot

# Import and tidy data -------------------------------------------------------------

dollo <- read_tsv(here("data_files_input_into_scripts", "Plot-Data-for-MS-FigS1.txt"))

# Clean up the data file a bit
dollo_data <-
  dollo %>%
  # filter(!str_detect(Edge_to, 'e|a|R|E')) %>% # If you want to remove branches leading to terminal tips
  separate(Edge, c("edge_fr", "edge_to"), sep="->") %>% 
  # filter(!str_detect(Change, '-')) # if ambiguous changes (-) are in data file, they must be removed
  separate(DataSetName, c("taxon", "method", "samp_depth", "parsimony"), sep="-") %>% 
  unite("dataset", taxon:method, sep="") # for ease, merge taxon and method cols together

# Calculate proportions of total numbers of changes
prop_dollo <-
  dollo_data %>% 
  group_by(dataset, samp_depth) %>% 
  mutate(sums = sum(n)) %>% 
  mutate(prop = n/sums)

# Standardize node numbering ----------------------------------------------

# Add additional column with node numbers that correspond with master tree numbering
# from p.epi and p.rana plots
epi_data <-
  prop_dollo %>%
  filter(dataset=="epi2brad" | dataset=="epiddrad") %>% 
  mutate(node = case_when(
    edge_to == "node11" ~ 15,
    edge_to == "node12" ~ 16,
    edge_to == "node13" ~ 14,
    edge_to == "node14" ~ 17,
    edge_to == "node15" ~ 13,
    edge_to == "node16" ~ 18,
    edge_to == "node17" ~ 12,
    edge_to == "node18" ~ 11,
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
  prop_dollo %>%
  filter(dataset=="rana2brad" | dataset=="ranaddrad") %>% 
  mutate(node = case_when(
    edge_to == "node11" ~ 17,
    edge_to == "node12" ~ 18,
    edge_to == "node13" ~ 16,
    edge_to == "node14" ~ 11,
    edge_to == "node15" ~ 12,
    edge_to == "node16" ~ 15,
    edge_to == "node17" ~ 13,
    edge_to == "node18" ~ 14,
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


# Format data for nodebars ------------------------------------------------

# Get the data into the correct format for the nodebars (ggtree) function

nodebar_format <- function (data, dataset, col, treepart) {
  dataset = enquo(dataset)
  treepart = enquo(treepart)
  col = enquo(col)
  data %>% 
    filter(dataset==!!dataset) %>% 
    select(samp_depth, !!col, !!treepart, Change) %>% 
    spread(samp_depth, !!col)
}

## EPI
epi2b_data_dollo <- epi_data %>% nodebar_format(dataset = "epi2brad", col = prop, treepart = node)
epidd_data_dollo <- epi_data %>% nodebar_format(dataset = "epiddrad", col = prop, treepart = node)

## RANA
rana2b_data_dollo <- rana_data %>% nodebar_format(dataset = "rana2brad", col = prop, treepart = node)
ranadd_data_dollo <- rana_data %>% nodebar_format(dataset = "ranaddrad", col = prop, treepart = node)


# (1) Build the base trees ---------------------------------------------------------

# Use Fig7&S1_PAUP_analysis.R to build p.epi and p.rana
  
# Adapt the nodebar function ---------------------------------------------------

## This alters the existing nodebar function (ggtree) with a custom aesthetic

nodebar_eac <- function (data, cols, alpha = 1, position="stack", ylim) {
  if (!"node" %in% colnames(data)) {
    stop("data should have a column 'node'...")
  }
  type <- value <- NULL # creates empty values for type and value, to be filled in next step
  ldf <- gather(data, type, value, !!cols) %>% split(., .$node) # creates sep tibbles for each node w type + value cols
  bars <- lapply(ldf, function(df)
    ggplot(df, aes(x = type, y = value, fill = Change)) +
      geom_bar(stat = "identity", alpha = alpha, position = 'stack') +
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
            plot.margin = margin(1,1,-1.5,1,"mm"), # t, r, b, l
            axis.text = element_blank()) +
      scale_x_discrete(expand = c(0,0)) +
      scale_y_continuous(expand = c(0,0), limits=ylim) +
      scale_fill_manual(values = c("0<=>1"="#91d5de", "1==>0"="#eaae37", "1<=>0"="#b4674e", "0==>1"="#2e8289"))
  )
  return(bars)
}


# Run function to build plots ---------------------------------------------

## Don't build plot for root branch; this needs to be enlarged and centred later
epi2b_bars_dollo <-
  epi2b_data_dollo %>% 
  filter(node!=11) %>% 
  nodebar_eac(cols=4:7, ylim=c(0,0.24)) # ymax=26789; propmax=0.232

epidd_bars_dollo <-
  epidd_data_dollo %>% 
  filter(node!=11) %>% 
  nodebar_eac(cols=4:7, ylim=c(0,0.3)) # ymax=103680; propmax=0.273

######## RANA
## Don't build plot for root branch; this needs to be enlarged and centred
rana2b_bars_dollo <-
  rana2b_data_dollo %>% 
  filter(node!=11) %>% 
  nodebar_eac(cols=4:7, ylim=c(0,0.24)) # ymax=28396; propmax=0.236

ranadd_bars_dollo <-
  ranadd_data_dollo %>% 
  filter(node!=11) %>% 
  nodebar_eac(cols=4:7, ylim=c(0,0.2)) # ymax=120896; propmax=0.198


# (2) Place small plots on trees -------------------------------------------------------

epi2b_plot_dollo <-
  inset(p.epi, epi2b_bars_dollo, width=0.3, height=0.1, vjust=-.475, hjust=.5)

epidd_plot_dollo <-
  inset(p.epi, epidd_bars_dollo, width=0.3, height=0.1, vjust=-.475, hjust=.5)

rana2b_plot_dollo <-
  inset(p.rana, rana2b_bars_dollo, width=0.3, height=0.1, vjust=-.475, hjust=.5)

ranadd_plot_dollo <-
  inset(p.rana, ranadd_bars_dollo, width=0.3, height=0.1, vjust=-.475, hjust=.5)


# Build blown up plot -------------------------------

# This will build a zoomed in plot for one of nodes which
# will be a reference for the yaxis scale

dollo_blownup <- function (data, yaxis, ylim, ytitle) {
  data %>% 
    ggplot(aes(x=samp_depth, y=yaxis, fill=Change)) +
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
    scale_x_discrete(labels=c("t1", "t2", "t3", "total")) +
    scale_fill_manual(values=c("0<=>1"="#91d5de", "1==>0"="#eaae37", "1<=>0"="#b4674e", "0==>1"="#2e8289")) +
    ylab(ytitle) +
    xlab("Sampling depth")
}

# For Epi, this is node 14
epi2b_inset_dollo <-
  epi_data %>% 
  filter(dataset=="epi2brad" & node==14) %>% 
  dollo_blownup(yaxis=.$prop, ylim = c(0,0.24), ytitle = "Prop. total changes")

epidd_inset_dollo <-
  epi_data %>% 
  filter(dataset=="epiddrad" & node==14) %>% 
  dollo_blownup(yaxis=.$prop, ylim = c(0,0.3), ytitle = "Prop. total changes")

########################
# For Rana, node 13
rana2b_inset_dollo <-
  rana_data %>% 
  filter(dataset=="rana2brad" & node==13) %>% 
  dollo_blownup(yaxis=.$prop, ylim=c(0,0.24), ytitle = "Prop. total changes")

ranadd_inset_dollo <-
  rana_data %>% 
  filter(dataset=="ranaddrad" & node==13) %>% 
  dollo_blownup(yaxis=.$prop, ylim=c(0,0.2), ytitle = "Prop. total changes")


# Build plot for root edges -----------------------------------------------

# Build root plots
# These are node 11 for both Epi and Rana

root_plots_dollo <- function (data, yaxis, ylim) {
  yaxis = enquo(yaxis)
  ggplot(data, aes(x=samp_depth, y=!!yaxis, fill=Change)) +
    geom_bar(stat="identity", position="stack") +
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
    scale_fill_manual(values=c("0<=>1"="#91d5de", "1==>0"="#eaae37", "1<=>0"="#b4674e", "0==>1"="#2e8289")) 
}

## change yaxis=prop to yaxis=n if you want raw number plotted

epi2b_root_dollo <-
  epi_data %>% 
  filter(dataset=="epi2brad" & node==11) %>% 
  root_plots_dollo(yaxis=prop, ylim=c(0,0.24))

epidd_root_dollo <-
  epi_data %>% 
  filter(dataset=="epiddrad" & node==11) %>% 
  root_plots_dollo(yaxis=prop, ylim=c(0,0.3))

########################
# For Rana
rana2b_root_dollo <-
  rana_data %>% 
  filter(dataset=="rana2brad" & node==11) %>% 
  root_plots_dollo(yaxis=prop, ylim=c(0,0.24))

ranadd_root_dollo <-
  rana_data %>% 
  filter(dataset=="ranaddrad" & node==11) %>% 
  root_plots_dollo(yaxis=prop, ylim=c(0,0.2))


# (3) Place blown-up and root plots onto main fig --------------------------

### EPIPEDOBATES
epi_final_plot <- function (base_plot, blowup, root, title) {
  base_plot +
    draw_plot(blowup, x=-1, y=8.5, width=3, height=3) +
    draw_plot(root, x=-.9, y=2.5, width=2, height=2) +
    geom_segment(aes(x=2.15, y=8.525, xend=1.575, yend=8.525), size=1, color="gray48") +
    geom_segment(aes(x=2.15, y=9.4, xend=1.8, yend=11.49), size=1, color="gray48") +
    theme(plot.title = element_text(face="bold", hjust = 0.5, size=16)) +
    ggtitle(title)
}

# Build the plots!
final_epi2b_dollo <-
  epi_final_plot(base_plot = epi2b_plot_dollo, blowup = epi2b_inset_dollo, root = epi2b_root_dollo, title = "2bRAD")
final_epidd_dollo <-
  epi_final_plot(base_plot = epidd_plot_dollo, blowup = epidd_inset_dollo, root = epidd_root_dollo, title = "ddRAD")

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

# Build the plots!
final_rana2b_dollo <-
  rana_final_plot(base_plot = rana2b_plot_dollo, blowup = rana2b_inset_dollo, root = rana2b_root_dollo, title = "2bRAD")
final_ranadd_dollo <-
  rana_final_plot(base_plot = ranadd_plot_dollo, blowup = ranadd_inset_dollo, root = ranadd_root_dollo, title = "ddRAD")


# Save four plots ------------------------------------------------------------

## Proportions
final_epi2b_dollo
ggsave("epi2b_dollo_props.pdf", width=8.36, height = 11.89)

final_epidd_dollo
ggsave("epidd_dollo_props.pdf", width=8.36, height = 11.89)

final_rana2b_dollo
ggsave("rana2b_dollo_props.pdf", width=8.9, height = 11.89)

final_ranadd_dollo
ggsave("ranadd_dollo_props.pdf", width=8.9, height = 11.89)

## If you've plotted sums (i.e., Fig S1 and not props as in Fig 9), use following dimensions

# ggsave("epi2b_dollo.pdf", width=8.36, height = 11.89)
# ggsave("epidd_dollo.pdf", width=8.36, height = 11.89)
# ggsave("rana2b_dollo.pdf", width=8.9, height = 11.89)
# ggsave("ranadd_dollo.pdf", width=8.9, height = 11.89)