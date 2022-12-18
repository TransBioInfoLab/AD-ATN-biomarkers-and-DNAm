#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Article:
# Distinct CSF biomarker-associated DNA methylation in Alzheimer's disease and cognitively normal subjects
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Author: Wei Zhang, Lily Wang
# Date: 07 Nov. 2022
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
dir.base <- "."
dir.results <- file.path(dir.base, "analysis-results/ADNI/LM/")
dir.results.adj <- file.path(dir.results, "adjusted/")
dir.plots.miami <- file.path(dir.base, "analysis-results/plots/miami")
dir.create(dir.plots.miami)
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#Load library
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
library(miamiplot)
library(tidyverse)
#-------------------------------------------------------------------------------
# Load regression results
#-------------------------------------------------------------------------------
DX <- c("AD","CN")
file.list.lm.names <- list.files(
  path = dir.results.adj,
  pattern = paste0(DX, "_results",collapse = "|"),
  full.names = T,
  recursive = T
)

file.list.lm <- file.list.lm.names %>% map(~readr::read_csv(.))

names(file.list.lm) <- unlist(
  paste0(
    "ADNI",
    file.list.lm.names %>% map(~stringr::str_match(.,"_AD|_CN")),
    file.list.lm.names %>% map(~stringr::str_extract(.,"_ABETA|_PTAU|_TAU"))
  )
)
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#Abeta
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
file.list.lm.abeta <- file.list.lm[grep("ABETA", names(file.list.lm))]
single.cpgs.abeta <- file.list.lm.abeta %>% 
  rbindlist(.,idcol = "ID") %>% mutate(ID = gsub("ADNI_","",ID)) %>% 
  separate(ID, into = "Sample", sep = "_",  remove = T)

plot.miami.single.cpgs <- single.cpgs.abeta %>% dplyr::mutate(
  "CpG" = cpg,
  "Sample" = Sample,
  "pVal.final" = pValue.bacon.z,
  "pos" = start,
  "chr" = gsub("chr", "",  seqnames) %>% as.numeric,
  "UCSC_RefGene_Name" = UCSC_RefGene_Name,
  .keep = "none"
)
plot_data <- prep_miami_data(
  data = plot.miami.single.cpgs, split_by = "Sample", split_at = "CN", p = "pVal.final"
)
#-------------------------------------------------------------------------------
# Label
#-------------------------------------------------------------------------------
labels_cn <- plot_data$upper %>%
  filter(pVal.final < 1e-05) %>% 
  dplyr::select(rel_pos, logged_p, UCSC_RefGene_Name) %>%
  dplyr::rename(label = UCSC_RefGene_Name) %>%
  tidyr::separate_rows(label, sep = ";") %>%
  unique() %>%
  arrange(desc(logged_p)) %>% na.omit()
labels_cn <- labels_cn[labels_cn$label != "",]
labels_cn <- labels_cn[!duplicated(labels_cn$label),]

labels_ad <- plot_data$lower %>%
  filter(pVal.final < 1e-05) %>% 
  dplyr::select(rel_pos, logged_p, UCSC_RefGene_Name) %>%
  dplyr::rename(label = UCSC_RefGene_Name) %>%
  tidyr::separate_rows(label, sep = ";") %>%
  unique() %>%
  arrange(desc(logged_p)) %>% na.omit()
labels_ad <- labels_ad[labels_ad$label != "",]
labels_ad <- labels_ad[!duplicated(labels_ad$label),]
#-------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------
plot <- ggmiami(
  data = plot.miami.single.cpgs, 
  split_by = "Sample", 
  split_at = "CN", 
  p = "pVal.final", 
  suggestive_line_color = "red",
  upper_ylab = "CN Sample",
  lower_ylab = "AD Sample",
  genome_line = NULL,
  top_n_hits = 20,
  suggestive_line = 10^-5,
  upper_labels_df = labels_cn,
  lower_labels_df = labels_ad
)

library(gtable)
ggplot2::ggsave(
  plot = gtable_add_padding(plot, unit(1, "cm")),
  filename = file.path(dir.plots.miami,"ADNI_ABETA_miami_plot.pdf"),
  width = 7,
  height = 6
)
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#PTAU
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
file.list.lm.ptau <- file.list.lm[grep("PTAU", names(file.list.lm))]
single.cpgs.ptau <- file.list.lm.ptau %>% 
  rbindlist(.,idcol = "ID") %>% mutate(ID = gsub("ADNI_","",ID)) %>% 
  separate(ID, into = "Sample", sep = "_",  remove = T)

plot.miami.single.cpgs <- single.cpgs.ptau %>% dplyr::mutate(
  "CpG" = cpg,
  "Sample" = Sample,
  "pVal.final" = pValue.bacon.z,
  "pos" = start,
  "chr" = gsub("chr", "",  seqnames) %>% as.numeric,
  "UCSC_RefGene_Name" = UCSC_RefGene_Name,
  .keep = "none"
)
plot_data <- prep_miami_data(
  data = plot.miami.single.cpgs, split_by = "Sample", split_at = "CN", p = "pVal.final"
)
#-------------------------------------------------------------------------------
# Label
#-------------------------------------------------------------------------------
labels_cn <- plot_data$upper %>%
  filter(pVal.final < 1e-05) %>% 
  dplyr::select(rel_pos, logged_p, UCSC_RefGene_Name) %>%
  dplyr::rename(label = UCSC_RefGene_Name) %>%
  tidyr::separate_rows(label, sep = ";") %>%
  unique() %>%
  arrange(desc(logged_p)) %>% na.omit()
labels_cn <- labels_cn[labels_cn$label != "",]
labels_cn <- labels_cn[!duplicated(labels_cn$label),]

labels_ad <- plot_data$lower %>%
  filter(pVal.final < 1e-05) %>% 
  dplyr::select(rel_pos, logged_p, UCSC_RefGene_Name) %>%
  dplyr::rename(label = UCSC_RefGene_Name) %>%
  tidyr::separate_rows(label, sep = ";") %>%
  unique() %>%
  arrange(desc(logged_p)) %>% na.omit()
labels_ad <- labels_ad[labels_ad$label != "",]
labels_ad <- labels_ad[!duplicated(labels_ad$label),]
#-------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------
plot <- ggmiami(
  data = plot.miami.single.cpgs, 
  split_by = "Sample", 
  split_at = "CN", 
  p = "pVal.final", 
  suggestive_line_color = "red",
  upper_ylab = "CN Sample",
  lower_ylab = "AD Sample",
  genome_line = NULL,
  top_n_hits = 20,
  suggestive_line = 10^-5,
  upper_labels_df = labels_cn,
  lower_labels_df = labels_ad
)

ggplot2::ggsave(
  plot = gtable_add_padding(plot, unit(1, "cm")),
  filename = file.path(dir.plots.miami,"ADNI_PTAU_miami_plot.pdf"),
  width = 7,
  height = 6
)
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
#TAU
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
file.list.lm.tau <- file.list.lm[grep("TAU", names(file.list.lm))]
single.cpgs.tau <- file.list.lm.tau %>% 
  rbindlist(.,idcol = "ID") %>% mutate(ID = gsub("ADNI_","",ID)) %>% 
  separate(ID, into = "Sample", sep = "_",  remove = T)

plot.miami.single.cpgs <- single.cpgs.tau %>% dplyr::mutate(
  "CpG" = cpg,
  "Sample" = Sample,
  "pVal.final" = pValue.bacon.z,
  "pos" = start,
  "chr" = gsub("chr", "",  seqnames) %>% as.numeric,
  "UCSC_RefGene_Name" = UCSC_RefGene_Name,
  .keep = "none"
)
plot_data <- prep_miami_data(
  data = plot.miami.single.cpgs, split_by = "Sample", split_at = "CN", p = "pVal.final"
)
#-------------------------------------------------------------------------------
# Label
#-------------------------------------------------------------------------------
labels_cn <- plot_data$upper %>%
  filter(pVal.final < 1e-05) %>% 
  dplyr::select(rel_pos, logged_p, UCSC_RefGene_Name) %>%
  dplyr::rename(label = UCSC_RefGene_Name) %>%
  tidyr::separate_rows(label, sep = ";") %>%
  unique() %>%
  arrange(desc(logged_p)) %>% na.omit()
labels_cn <- labels_cn[labels_cn$label != "",]
labels_cn <- labels_cn[!duplicated(labels_cn$label),]

labels_ad <- plot_data$lower %>%
  filter(pVal.final < 1e-05) %>% 
  dplyr::select(rel_pos, logged_p, UCSC_RefGene_Name) %>%
  dplyr::rename(label = UCSC_RefGene_Name) %>%
  tidyr::separate_rows(label, sep = ";") %>%
  unique() %>%
  arrange(desc(logged_p)) %>% na.omit()
labels_ad <- labels_ad[labels_ad$label != "",]
labels_ad <- labels_ad[!duplicated(labels_ad$label),]
#-------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------
plot <- ggmiami(
  data = plot.miami.single.cpgs, 
  split_by = "Sample", 
  split_at = "CN", 
  p = "pVal.final", 
  suggestive_line_color = "red",
  upper_ylab = "CN Sample",
  lower_ylab = "AD Sample",
  genome_line = NULL,
  top_n_hits = 20,
  suggestive_line = 10^-5,
  upper_labels_df = labels_cn,
  lower_labels_df = labels_ad
)

ggplot2::ggsave(
  plot = gtable_add_padding(plot, unit(1, "cm")),
  filename = file.path(dir.plots.miami,"ADNI_TAU_miami_plot.pdf"),
  width = 7,
  height = 6
)
