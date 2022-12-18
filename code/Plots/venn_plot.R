dir.base <- "."
dir.combp <- file.path(dir.base, "Mike/")
dir.results <- file.path(dir.base, "analysis-results/ADNI/LM/")
dir.results.adj <- file.path(dir.results, "adjusted/")
dir.combp.results <- file.path(dir.base, "analysis-results/comb-p/")

library(tidyverse)
library(ggpubr)
library(ggvenn)

# Load LM results

DX <- c("AD", "CN")
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

# Load comb-p DMR results

load(file.path(dir.combp.results,"ADNI_all_DMR_annotated_filtered_zsidakp_zp_same_direction_results.rda"))

# 5. Venn Diagram

dir.venn <- file.path(dir.base, "analysis-results/plots/venn/")
dir.create(dir.venn)

# 5.1 Separate by CSF

## Single cpgs

### ABETA

venn.list <- list(
  "AD" = file.list.lm$ADNI_AD_ABETA %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "CN" = file.list.lm$ADNI_CN_ABETA %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg)
)

g1 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#0073C2FF","#CD534CFF")
)

### PTAU

venn.list <- list(
  "AD" = file.list.lm$ADNI_AD_PTAU %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "CN" = file.list.lm$ADNI_CN_PTAU %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg)
)

g2 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#0073C2FF","#CD534CFF")
)

### TAU

venn.list <- list(
  "AD" = file.list.lm$ADNI_AD_TAU %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "CN" = file.list.lm$ADNI_CN_TAU %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg)
)

g3 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#0073C2FF","#CD534CFF")
)

plot_single_cpgs <- ggarrange(
  g1,g2,g3,
  ncol = 3,
  labels = list("ABETA","PTAU","TAU"),
  font.label = list(face = "plain"),
  widths = c(2,2,2)
)

## cpgs in DMRs

### ABETA

overlap_two_group <- findOverlapPairs(dmr.filter.list$ADNI_AD_ABETA$region %>% MethReg::make_granges_from_names(),
                                      dmr.filter.list$ADNI_CN_ABETA$region %>% MethReg::make_granges_from_names())

#### Substitude them with the same names
overlap_in_1 <- overlap_two_group@first %>% MethReg::make_names_from_granges()
overlap_in_2 <- overlap_two_group@second %>% MethReg::make_names_from_granges()

region1 <- dmr.filter.list$ADNI_AD_ABETA$region
region2 <- dmr.filter.list$ADNI_CN_ABETA$region
region2[region2 %in% overlap_in_2] <- overlap_in_1

venn.list <- list(
  "AD" = region1,
  "CN" = region2
)

g1 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#EFC000FF","#868686FF")
)

### PTAU

overlap_two_group <- findOverlapPairs(dmr.filter.list$ADNI_AD_PTAU$region %>% MethReg::make_granges_from_names(),
                                      dmr.filter.list$ADNI_CN_PTAU$region %>% MethReg::make_granges_from_names())

#### Substitude them with the same names
overlap_in_1 <- overlap_two_group@first %>% MethReg::make_names_from_granges()
overlap_in_2 <- overlap_two_group@second %>% MethReg::make_names_from_granges()

region1 <- dmr.filter.list$ADNI_AD_PTAU$region
region2 <- dmr.filter.list$ADNI_CN_PTAU$region
region2[region2 %in% overlap_in_2] <- overlap_in_1

venn.list <- list(
  "AD" = region1,
  "CN" = region2
)

g2 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#EFC000FF","#868686FF")
)

### TAU

overlap_two_group <- findOverlapPairs(dmr.filter.list$ADNI_AD_TAU$region %>% MethReg::make_granges_from_names(),
                                      dmr.filter.list$ADNI_CN_TAU$region %>% MethReg::make_granges_from_names())

#### Substitude them with the same names
overlap_in_1 <- overlap_two_group@first %>% MethReg::make_names_from_granges()
overlap_in_2 <- overlap_two_group@second %>% MethReg::make_names_from_granges()

region1 <- dmr.filter.list$ADNI_AD_TAU$region
region2 <- dmr.filter.list$ADNI_CN_TAU$region
region2[region2 %in% overlap_in_2] <- overlap_in_1

venn.list <- list(
  "AD" = region1,
  "CN" = region2
)

g3 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#EFC000FF","#868686FF")
)

plot_cpgs_in_dmr <- ggarrange(
  g1,g2,g3,
  ncol = 3,
  labels = list("ABETA","PTAU","TAU"),
  font.label = list(face = "plain"),
  widths = c(2,2,2)
)

plot_single_and_dmr_cpgs <- ggarrange(
  plot_single_cpgs,
  plot_cpgs_in_dmr,
  labels = c("Sig. CpGs in single cpgs analysis", "Sig. DMRs in comb-p results"),
  nrow = 2,
  hjust = c(0,0),
  vjust = c(-1.5,-1.5),
  align = "v",
  font.label = list(face = "plain")
) + theme(plot.margin = margin(1.5,.5,0.1,.5, "cm")) 

ggsave(
  file.path(dir.venn, "/ADNI_sig_single_cpgs_and_dmrs_by_CSF_venn_plot.pdf"),
  width = 12,
  height = 8
)

## cpgs in DMRs vs single cpgs

### ABETA

venn.list <- list(
  "AD single" = file.list.lm$ADNI_AD_ABETA %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "CN single" = file.list.lm$ADNI_CN_ABETA %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "AD DMR" = str_split(dmr.filter.list$ADNI_AD_ABETA %>% pull(cpgs_in_region), ",") %>% unlist(),
  "CN DMR" = str_split(dmr.filter.list$ADNI_CN_ABETA %>% pull(cpgs_in_region), ",") %>% unlist()
)

g1 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#0073C2FF","#CD534CFF", "#EFC000FF","#868686FF")
)

### PTAU

venn.list <- list(
  "AD single" = file.list.lm$ADNI_AD_PTAU %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "CN single" = file.list.lm$ADNI_CN_PTAU %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "AD DMR" = str_split(dmr.filter.list$ADNI_AD_PTAU %>% pull(cpgs_in_region), ",") %>% unlist(),
  "CN DMR" = str_split(dmr.filter.list$ADNI_CN_PTAU %>% pull(cpgs_in_region), ",") %>% unlist()
)

g2 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#0073C2FF","#CD534CFF", "#EFC000FF","#868686FF")
)

### TAU

venn.list <- list(
  "AD single" = file.list.lm$ADNI_AD_TAU %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "CN single" = file.list.lm$ADNI_CN_TAU %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "AD DMR" = str_split(dmr.filter.list$ADNI_AD_TAU %>% pull(cpgs_in_region), ",") %>% unlist(),
  "CN DMR" = str_split(dmr.filter.list$ADNI_CN_TAU %>% pull(cpgs_in_region), ",") %>% unlist()
)

g3 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#0073C2FF","#CD534CFF", "#EFC000FF","#868686FF")
)

plot_single_vs_dmr <- ggarrange(
  g1,g2,g3,
  ncol = 3,
  labels = list("ABETA","PTAU","TAU"),
  font.label = list(face = "plain"),
  widths = c(2,2,2)
)

ggsave(
  file.path(dir.venn, "/ADNI_sig_single_cpgs_overlap_with_dmrs_by_CSF_venn_plot.pdf"),
  width = 15,
  height = 4
)

### ABETA

venn.list <- list(
  "AD single" = file.list.lm$ADNI_AD_ABETA %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "AD DMR" = str_split(dmr.filter.list$ADNI_AD_ABETA %>% pull(cpgs_in_region), ",") %>% unlist()
)

s <- file.list.lm$ADNI_AD_ABETA %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg)
plyr::l_ply(
  str_split(dmr.filter.list$ADNI_AD_ABETA %>% pull(cpgs_in_region), ","),
  .fun=function(l) print(any(l %in% s))
) #6 #1 #2 #5

g1 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#0073C2FF", "#EFC000FF")
)

### PTAU

venn.list <- list(
  "AD single" = file.list.lm$ADNI_AD_PTAU %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "AD DMR" = str_split(dmr.filter.list$ADNI_AD_PTAU %>% pull(cpgs_in_region), ",") %>% unlist()
)

g2 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#0073C2FF", "#EFC000FF")
)

### TAU

venn.list <- list(
  "AD single" = file.list.lm$ADNI_AD_TAU %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "AD DMR" = str_split(dmr.filter.list$ADNI_AD_TAU %>% pull(cpgs_in_region), ",") %>% unlist()
)

g3 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#0073C2FF", "#EFC000FF")
)

plot_single_vs_dmr_ad <- ggarrange(
  g1,g2,g3,
  ncol = 3,
  labels = list("ABETA","PTAU","TAU"),
  font.label = list(face = "plain"),
  widths = c(2,2,2)
)

### ABETA

venn.list <- list(
  "CN single" = file.list.lm$ADNI_CN_ABETA %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "CN DMR" = str_split(dmr.filter.list$ADNI_CN_ABETA %>% pull(cpgs_in_region), ",") %>% unlist()
)

g1 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#CD534CFF", "#868686FF")
)

### PTAU

venn.list <- list(
  "CN single" = file.list.lm$ADNI_CN_PTAU %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "CN DMR" = str_split(dmr.filter.list$ADNI_CN_PTAU %>% pull(cpgs_in_region), ",") %>% unlist()
)

g2 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#CD534CFF", "#868686FF")
)

### TAU

venn.list <- list(
  "CN single" = file.list.lm$ADNI_CN_TAU %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "CN DMR" = str_split(dmr.filter.list$ADNI_CN_TAU %>% pull(cpgs_in_region), ",") %>% unlist()
)

g3 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#CD534CFF", "#868686FF")
)

plot_single_vs_dmr_cn <- ggarrange(
  g1,g2,g3,
  ncol = 3,
  labels = list("ABETA","PTAU","TAU"),
  font.label = list(face = "plain"),
  widths = c(2,2,2)
)

plot_single_and_dmr_cpgs_v2 <- ggarrange(
  plot_single_vs_dmr_ad,
  plot_single_vs_dmr_cn,
  labels = c("AD group", "CN group"),
  nrow = 2,
  hjust = c(0,0),
  vjust = c(-1.5,-1.5),
  align = "v",
  font.label = list(face = "plain")
) + theme(plot.margin = margin(1.5,.5,0.1,.5, "cm")) 

ggsave(
  file.path(dir.venn, "/ADNI_sig_single_cpgs_overlap_with_dmrs_by_CSF_venn_plot_V2.pdf"),
  width = 12,
  height = 8
)

# 5.1 Separate by group

## Single cpgs

### AD

venn.list <- list(
  "ABETA" = file.list.lm$ADNI_AD_ABETA %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "PTAU" = file.list.lm$ADNI_AD_PTAU %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "TAU" = file.list.lm$ADNI_AD_TAU %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg)
)

g1 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#0073C2FF","#CD534CFF", "#EFC000FF")
)

### CN

venn.list <- list(
  "ABETA" = file.list.lm$ADNI_CN_ABETA %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "PTAU" = file.list.lm$ADNI_CN_PTAU %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg),
  "TAU" = file.list.lm$ADNI_CN_TAU %>% filter(pValue.bacon.z < 1E-05) %>% pull(cpg)
)

g2 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#0073C2FF","#CD534CFF", "#EFC000FF")
)


plot_single_cpgs_by_group <- ggarrange(
  g1,g2,
  ncol = 2,
  labels = list("AD","CN"),
  font.label = list(face = "plain")
) + theme(plot.margin = margin(.5,.5, 1,.5, "cm")) 

## cpgs in DMRs

### AD

overlap_two_group <- findOverlapPairs(dmr.filter.list$ADNI_AD_PTAU$region %>% MethReg::make_granges_from_names(),
                                      dmr.filter.list$ADNI_AD_TAU$region %>% MethReg::make_granges_from_names())

#### Substitude them with the same names
overlap_in_1 <- overlap_two_group@first %>% MethReg::make_names_from_granges()
overlap_in_2 <- overlap_two_group@second %>% MethReg::make_names_from_granges()

region1 <- dmr.filter.list$ADNI_AD_ABETA$region
region2 <- dmr.filter.list$ADNI_AD_PTAU$region
region3 <- dmr.filter.list$ADNI_AD_TAU$region
region3[region3 %in% overlap_in_2] <- overlap_in_1

venn.list <- list(
  "ABETA" = region1,
  "PTAU" = region2,
  "TAU" = region3
)

g1 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#0073C2FF","#CD534CFF", "#EFC000FF")
)

### CN

overlap_two_group <- findOverlapPairs(dmr.filter.list$ADNI_CN_ABETA$region %>% MethReg::make_granges_from_names(),
                                      dmr.filter.list$ADNI_CN_PTAU$region %>% MethReg::make_granges_from_names())
overlap_two_group2 <- findOverlapPairs(dmr.filter.list$ADNI_CN_TAU$region %>% MethReg::make_granges_from_names(),
                                       dmr.filter.list$ADNI_CN_PTAU$region %>% MethReg::make_granges_from_names())

#### Substitude them with the same names
overlap_in_1 <- overlap_two_group@first %>% MethReg::make_names_from_granges()
overlap_in_2 <- overlap_two_group@second %>% MethReg::make_names_from_granges()

overlap_in_3 <- overlap_two_group2@first %>% MethReg::make_names_from_granges()
overlap_in_4 <- overlap_two_group2@second %>% MethReg::make_names_from_granges()

region1 <- dmr.filter.list$ADNI_CN_ABETA$region
region2 <- dmr.filter.list$ADNI_CN_PTAU$region
region3 <- dmr.filter.list$ADNI_CN_TAU$region
region1[region1 %in% overlap_in_1] <- overlap_in_2
region3[region3 %in% overlap_in_3] <- overlap_in_4

venn.list <- list(
  "ABETA" = region1,
  "PTAU" = region2,
  "TAU" = region3
)

g2 <- ggvenn(
  venn.list,
  set_name_size = 4,
  text_size = 3,
  fill_color = c("#0073C2FF","#CD534CFF", "#EFC000FF")
)


plot_cpgs_in_dmr_by_group <- ggarrange(
  g1,g2,
  ncol = 2,
  labels = list("AD","CN"),
  font.label = list(face = "plain")
) + theme(plot.margin = margin(1,.5, .5,.5, "cm")) 

plot_single_and_dmr_cpgs_by_group <- ggarrange(
  plot_single_cpgs_by_group,
  plot_cpgs_in_dmr_by_group,
  labels = c("Sig. CpGs in single cpgs analysis", "Sig. DMRs in comb-p results"),
  nrow = 2,
  hjust = c(0,0),
  vjust = c(-1.5,-1.5),
  label.y = 0.9,
  align = "hv",
  font.label = list(face = "plain")
) + theme(plot.margin = margin(1.5,.5,0.1,.5, "cm")) 

ggsave(
  file.path(dir.venn, "/ADNI_sig_single_cpgs_and_dmrs_by_group_venn_plot.pdf"),
  width = 10,
  height = 10
)


