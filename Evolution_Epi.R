setwd("~/Box Sync/Epipedobates project/BIOINFORMATICS/ANALYSIS/epi_rana_radseq")

library(tidyverse)
library(cowplot)

# Upload ddrad files; three per dataset
epi_ddrad_snpdist <- read_csv("epiddrad_snpdist.csv")
epi_ddrad_sumstats <- read_csv("epiddrad_sumstats.csv")
epi_ddrad_cov <- read_csv("epiddrad_coverage.csv")
rana_ddrad_snpdist <- read_csv("ranaddrad_snpdist.csv")
rana_ddrad_sumstats <- read_csv("ranaddrad_sumstats.csv")
rana_ddrad_cov <- read_csv("ranaddrad_coverage.csv")

# Upload 2brad files
rana_2brad_snpdist <- read_csv("rana2brad_snpdist.csv")
rana_2brad_sumstats <- read_csv("rana2brad_sumstats.csv")
rana_2brad_cov <- read_csv("rana2brad_coverage.csv")
epi_2brad_snpdist <- read_csv("epi2brad_snpdist.csv")
epi_2brad_sumstats <- read_csv("epi2brad_sumstats.csv")
epi_2brad_cov <- read_csv("epi2brad_coverage.csv")

# Easier labeling of figures
ranaclustlabs = c("80", "81", "84", "85", "86", "87", "88", "89", "90", "91", "93", "95")
epiclustlabs = c("80", "81", "82", "83", "84", "85", "86", "87", "88", "89", "90", "91", "92", "93", "94", "95")
epiclustvals = unique(epi_2brad_cov$clust_threshold)
ranalabs = c("Rber_a","Rber_b","Rber","Rbla_4","Rbla_5","Rchi_a","Rchi_b","Rchi","Rneo_0","Rneo_7","Rsph_0","Rsph_4")
ranaclustvals = unique(rana_2brad_cov$clust_threshold)
epilabs = c("Ahah_a","Ahah_b","Ahah","Eant","Eant_a","Eant_b","Ebou_3","Ebou_6","Etri_6","Etri_2","Snub_8","Snub_9")
epiddclustlabs = c("80", "81", "82", "83", "84", "85", "86", "87", "89", "90", "91")
ranaddclustvals = unique(rana_ddrad_cov$clust_threshold)
ranaddclustlabs = c("80", "81")

################################## PI SITES ################################## 
p.epi_2brad_snpdist_pis <-
  epi_2brad_snpdist %>% 
  group_by(clust_threshold) %>%
  filter(sum_pis == max(sum_pis)) %>%
  ggplot(aes(x=clust_threshold, y=sum_pis)) + 
  geom_col() +
  scale_x_discrete(labels=epiclustlabs) +
  scale_y_continuous(labels=scales::comma, limits=c(0,1200000)) +
  ggtitle("Epi_2brad")

p.epi_ddrad_snpdist_pis <-
  epi_ddrad_snpdist %>% 
  group_by(clust_threshold) %>%
  filter(sum_pis == max(sum_pis)) %>%
  ggplot(aes(x=clust_threshold, y=sum_pis)) + 
  geom_col() +
  scale_x_discrete(labels=epiddclustlabs) +
  scale_y_continuous(labels=scales::comma, limits=c(0,1200000)) +
  ggtitle("Epi_ddrad")

p.rana_2brad_snpdist_pis <-
  rana_2brad_snpdist %>% 
  group_by(clust_threshold) %>%
  filter(sum_pis == max(sum_pis)) %>%
  ggplot(aes(x=clust_threshold, y=sum_pis)) + 
  geom_col() +
  scale_x_discrete(labels=ranaclustlabs) +
  scale_y_continuous(labels=scales::comma, limits=c(0,1200000)) +
  ggtitle("Rana_2brad")

p.rana_ddrad_snpdist_pis <-
  rana_ddrad_snpdist %>% 
  group_by(clust_threshold) %>%
  filter(sum_pis == max(sum_pis)) %>%
  ggplot(aes(x=clust_threshold, y=sum_pis)) + 
  geom_col() +
  scale_x_discrete(labels=ranaddclustlabs) +
  scale_y_continuous(labels=scales::comma, limits=c(0,1200000)) +
  ggtitle("Rana_ddrad")

plot_grid(p.epi_2brad_snpdist_pis, p.epi_ddrad_snpdist_pis, p.rana_2brad_snpdist_pis, p.rana_ddrad_snpdist_pis)

################################## COVERAGE (above 3) ##########################
p.epi_2brad_cov <-
  epi_2brad_cov %>% 
  group_by(clust_threshold) %>%
  filter(number>3) %>%
  ggplot(aes(x=number, y=locus_coverage, group=clust_threshold, color=clust_threshold)) + 
  geom_line() +
  scale_y_continuous(limits=c(0,71000)) +
  ggtitle("Epi_2brad") +
  scale_color_discrete(labels=epiclustlabs)

p.epi_ddrad_cov<-
  epi_ddrad_cov %>% 
  group_by(clust_threshold) %>%
  filter(number>3) %>%
  ggplot(aes(x=number, y=locus_coverage, group=clust_threshold, color=clust_threshold)) + 
  geom_line() +
  scale_y_continuous(limits=c(0,8000)) +
  ggtitle("Epi_ddrad") +
  scale_color_discrete(labels=epiddclustlabs)

p.rana_2brad_cov <-
  rana_2brad_cov %>% 
  group_by(clust_threshold) %>%
  filter(number>3) %>%
  ggplot(aes(x=number, y=locus_coverage, group=clust_threshold, color=clust_threshold)) + 
  geom_line() +
  scale_y_continuous(limits=c(0,71000)) +
  ggtitle("Rana_2brad") +
  scale_color_discrete(labels=ranaclustlabs)

p.rana_ddrad_cov <-
  rana_ddrad_cov %>% 
  group_by(clust_threshold) %>%
  filter(number>3) %>%
  ggplot(aes(x=number, y=locus_coverage, group=clust_threshold, color=clust_threshold)) + 
  geom_line() +
  scale_y_continuous(limits=c(0,8000)) +
  ggtitle("Rana_ddrad") +
  scale_color_discrete(labels=ranaddclustlabs)

plot_grid(p.epi_2brad_cov, p.epi_ddrad_cov, p.rana_2brad_cov, p.rana_ddrad_cov)

################################## LOCI RECOVERED ################################## 
# p.rana_2brad_sumstats <-
#   rana_2brad_sumstats %>%
#   ggplot(aes(clust_threshold, loci_assembly, group=sample, color=sample)) +
#   geom_line() +
#   scale_x_discrete(labels=ranaclustlabs) +
#   scale_y_continuous(labels=scales::comma, limits=c(0,200000)) +
#   ggtitle("Rana_2bRAD")
#   
# p.epi_2brad_sumstats <-
#   epi_2brad_sumstats %>%
#   ggplot(aes(clust_threshold, loci_assembly, group=sample, color=sample)) +
#   geom_line() +
#   scale_x_discrete(labels=epiclustlabs) +
#   scale_y_continuous(labels=scales::comma, limits=c(0,200000)) +
#   ggtitle("Epi_2bRAD")
# 
# plot_grid(p.epi_2brad_sumstats, p.rana_2brad_sumstats)


################################## VARIABLE SITES ####################################
p.rana_2brad_snpdist <-
  rana_2brad_snpdist %>%
  group_by(clust_threshold) %>%
  filter(sum_var == max(sum_var)) %>%
ggplot(aes(x=clust_threshold, y=sum_var)) +
  geom_col() +
  scale_x_discrete(labels=ranaclustlabs) +
  scale_y_continuous(labels=scales::comma, limits=c(0,400000)) +
  ggtitle("Rana_2brad")

p.rana_ddrad_snpdist <-
  rana_ddrad_snpdist %>%
  group_by(clust_threshold) %>%
  filter(sum_var == max(sum_var)) %>%
  ggplot(aes(x=clust_threshold, y=sum_var)) +
  geom_col() +
  scale_x_discrete(labels=ranaddclustlabs) +
  scale_y_continuous(labels=scales::comma, limits=c(0,400000)) +
  ggtitle("Rana_ddrad")

p.epi_2brad_snpdist <-
  epi_2brad_snpdist %>%
  group_by(clust_threshold) %>%
  filter(sum_var == max(sum_var)) %>%
  ggplot(aes(x=clust_threshold, y=sum_var)) +
  geom_col() +
  scale_x_discrete(labels=epiclustlabs) +
  scale_y_continuous(labels=scales::comma, limits=c(0,400000)) +
  ggtitle("Epi_2brad")

p.epi_ddrad_snpdist <-
  epi_ddrad_snpdist %>%
  group_by(clust_threshold) %>%
  filter(sum_var == max(sum_var)) %>%
  ggplot(aes(x=clust_threshold, y=sum_var)) +
  geom_col() +
  scale_x_discrete(labels=epiddclustlabs) +
  scale_y_continuous(labels=scales::comma, limits=c(0,400000)) +
  ggtitle("Epi_ddrad")
# 
plot_grid(p.epi_2brad_snpdist, p.epi_ddrad_snpdist, p.rana_2brad_snpdist, p.rana_ddrad_snpdist)

################################## Calc INF SITES PER LOCUS ###############################
# Epi 2b
epi2b_maxpis <- 
  epi_2brad_snpdist %>%
  group_by(clust_threshold) %>%
  filter(sum_pis == max(sum_pis), sum_var == max(sum_var)) %>%
  select(sum_pis, sum_var, clust_threshold) %>%
  distinct(sum_pis, sum_var, clust_threshold)

epi2b_maxcov <-
  epi_2brad_cov %>%
  group_by(clust_threshold) %>%
  filter(sum_coverage == max(sum_coverage)) %>%
  select(sum_coverage, clust_threshold) %>%
  distinct(sum_coverage, clust_threshold)

epi2b <-
  epi2b_maxpis %>%
  left_join(epi2b_maxcov, by="clust_threshold")

epi2b <-
epi2b %>%
  mutate(stat_pis = sum_pis/sum_coverage) %>%
  mutate(stat_var = sum_var/sum_coverage)

# Epi dd
epidd_maxpis <- 
  epi_ddrad_snpdist %>%
  group_by(clust_threshold) %>%
  filter(sum_pis == max(sum_pis), sum_var == max(sum_var)) %>%
  select(sum_pis, sum_var, clust_threshold) %>%
  distinct(sum_pis, sum_var, clust_threshold)

epidd_maxcov <-
  epi_ddrad_cov %>%
  group_by(clust_threshold) %>%
  filter(sum_coverage == max(sum_coverage)) %>%
  select(sum_coverage, clust_threshold) %>%
  distinct(sum_coverage, clust_threshold)

epidd <-
  epidd_maxpis %>%
  left_join(epidd_maxcov, by="clust_threshold")

epidd <-
  epidd %>%
  mutate(stat_pis = sum_pis/sum_coverage) %>%
  mutate(stat_var = sum_var/sum_coverage)



# Rana 2b
rana2b_maxpis <- 
  rana_2brad_snpdist %>%
  group_by(clust_threshold) %>%
  filter(sum_pis == max(sum_pis), sum_var == max(sum_var)) %>%
  select(sum_pis, sum_var, clust_threshold) %>%
  distinct(sum_pis, sum_var, clust_threshold)

rana2b_maxcov <-
  rana_2brad_cov %>%
  group_by(clust_threshold) %>%
  filter(sum_coverage == max(sum_coverage)) %>%
  select(sum_coverage, clust_threshold) %>%
  distinct(sum_coverage, clust_threshold)

rana2b <-
  rana2b_maxpis %>%
  left_join(rana2b_maxcov, by="clust_threshold")

rana2b <-
  rana2b %>%
  mutate(stat_pis = sum_pis/sum_coverage) %>%
  mutate(stat_var = sum_var/sum_coverage)



# Rana dd
ranadd_maxpis <- 
  rana_ddrad_snpdist %>%
  group_by(clust_threshold) %>%
  filter(sum_pis == max(sum_pis), sum_var == max(sum_var)) %>%
  select(sum_pis, sum_var, clust_threshold) %>%
  distinct(sum_pis, sum_var, clust_threshold)

ranadd_maxcov <-
  rana_ddrad_cov %>%
  group_by(clust_threshold) %>%
  filter(sum_coverage == max(sum_coverage)) %>%
  select(sum_coverage, clust_threshold) %>%
  distinct(sum_coverage, clust_threshold)

ranadd <-
  ranadd_maxpis %>%
  left_join(epidd_maxcov, by="clust_threshold")

ranadd <-
  ranadd %>%
  mutate(stat_pis = sum_pis/sum_coverage) %>%
  mutate(stat_var = sum_var/sum_coverage)
################################## PLOT var&pi SITES PER LOCUS ###############################
p.epidd <-
ggplot(epidd, aes(clust_threshold, stat_var)) +
  geom_col() +
  scale_x_discrete(labels=epiddclustlabs) +
  ggtitle("Epi_ddRAD") +
  ylab("Variable sites / total shared loci")

p.epi2b <- 
  ggplot(epi2b, aes(clust_threshold, stat_var)) +
  geom_col() +
  scale_x_discrete(labels=epiclustlabs) +
  ggtitle("Epi_2bRAD") +
  ylab("Variable sites / total shared loci")

p.ranadd <-
  ggplot(ranadd, aes(clust_threshold, stat_var)) +
  geom_col() +
  scale_x_discrete(labels=ranaddclustlabs) +
  ggtitle("Rana_ddRAD") +
  ylab("Variable sites / total shared loci")

p.rana2b <- 
  ggplot(rana2b, aes(clust_threshold, stat_var)) +
  geom_col() +
  scale_x_discrete(labels=ranaclustlabs) +
  ggtitle("Rana_2bRAD") +
  ylab("Variable sites / total shared loci")

plot_grid(p.epi2b, p.epidd, p.rana2b)  
