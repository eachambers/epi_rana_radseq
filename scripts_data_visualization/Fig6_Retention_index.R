library(dplyr)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(here)

theme_set(theme_cowplot())

## The following code generates Fig. 6, which illustrates the retention indices
## and parsimony-informative sites for each sampling depth, taxon, and dataset.

##    FILES REQUIRED:
##          data_files_input_into_scripts/Retention_PIs.csv


#=============================================================================

# Import data -------------------------------------------------------------

# Read in file; these data are first generated from PAUP

data <- read_csv(here("data_files_input_into_scripts", "Retention_PIs.csv"))

# Add a column with the % parsimony inf sites:
data <-
  data %>%
  mutate(percent_PI = prop_pars_inf*100) %>% 
  select(-"...8")

# Build plots -----------------------------------------------

p_epi <-
  data %>%
  filter(taxon == "epi") %>%
  ggplot(aes(x = seq_depth), color = seq_depth) +
  geom_line(aes(y = percent_PI, group = dataset), size = .5) +
  geom_point(aes(y = percent_PI), size = 2.5) +
  facet_grid(~dataset) +
  scale_y_continuous(limits = c(0,70), breaks = c(0,20,40,60)) +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        strip.background = element_blank(),
        legend.position = "none") +
  ylab("Proportion of parsimony-informative characters (%)") +
  xlab("Sampling depth")

p_epi <-
  p_epi +
  geom_line(aes(y = retention_index*60, group = dataset), color = "darkgrey", linewidth = .5) +
  geom_point(aes(y = retention_index*60), size = 2.5, color = "darkgrey") +
  facet_grid(~dataset) +
  scale_y_continuous(limits = c(0,60), breaks = c(0,20,40,60))

## Add on a second axis to the right that represents the retention index
## You need to divide the axis by what the multiplication factor above so
## that it correctly represents those values
p_epi <-
  p_epi +
  scale_y_continuous(sec.axis = sec_axis(~./60))

#### now do the same for Rana ####

p_rana <-
  data %>%
  filter(taxon == "rana") %>%
  ggplot(aes(x = seq_depth), color = seq_depth) +
  geom_line(aes(y = percent_PI, group = dataset), linewidth = .5) +
  geom_point(aes(y = percent_PI), size = 2.5) +
  facet_grid(~dataset) +
  scale_y_continuous(limits = c(0,60), breaks = c(0,20,40,60)) +
  theme(axis.text = element_text(size = 9),
        axis.title = element_text(size = 10),
        strip.background = element_blank(),
        legend.position = "none") +
  ylab("Proportion of parsimony-informative characters (%)") +
  xlab("Sampling depth")

p_rana <-
  p_rana +
  geom_line(aes(y = retention_index*60, group = dataset), color = "darkgrey", linewidth = .5) +
  geom_point(aes(y = retention_index*60), size = 2.5, color = "darkgrey") +
  facet_grid(~dataset) +
  scale_y_continuous(limits = c(0,60), breaks = c(0,20,40,60))

p_rana <-
  p_rana +
  scale_y_continuous(sec.axis = sec_axis(~./60))

# Combine both retention plots together
p_retPIs <- plot_grid(p_rana, p_epi, nrow = 2, align = 'v')


# Export and save plot ---------------------------------------------------

# Save retention plots
# base_height is height of one of plots in inches
# base_aspect_ratio is ratio between width/height of one plot
save_plot("Fig6_RetentionPIs.pdf", p_retPIs, base_width = 6.9, base_height = 6.8)
