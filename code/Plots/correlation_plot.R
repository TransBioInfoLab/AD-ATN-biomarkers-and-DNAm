library(tidyverse)
library(ggpubr)
dir.base <- "."
dir.results <- file.path(dir.base, "analysis-results/ADNI/LM/")
dir.results.sig <- file.path(dir.results, "sig_cpgs/")
dir.results.adj <- file.path(dir.results, "adjusted/")
dir.plot <- file.path(dir.base, "analysis-results/plots/correlation_plot")
dir.create(dir.plot)
# Load results
DX <- c("AD", "CN")
file.list <- list.files(
  path = dir.results.adj,
  pattern = paste0(DX, "_results",collapse = "|"),
  full.names = T,
  recursive = T
)
## Merge different group of samples in one dataframe
merge.list <- function(file.list, csf = "ABETA"){
  DX <- c("AD", "CN")
  # read results
  print(grep(paste0("_",csf), file.list, value = T))
  results.list <- grep(paste0("_",csf), file.list, value = T) %>% purrr::map(~readr::read_csv(., show_col_types = F)) %>% 
    purrr::map(~dplyr::select(
      ., 
      cpg, Estimate.bacon, StdErr.bacon, pValue.bacon.z, fdr.bacon.z
    )
    ) 
  # rename dataframe
  dataset <- paste0(csf,"_",stringr::str_extract(
    stringr::str_extract(grep(paste0("_",csf), file.list, value = T), paste0(DX, "_results",collapse = "|")), 
    paste0(DX,collapse = "|")
  ))
  names(results.list) <- dataset
  # select bacon statistics
  results.list <- plyr::alply(
    dataset,
    .margin = 1, 
    .fun = function(d){
      results.list[[d]] %>% rename_with(
        .fn = function(x) {paste0(d,"_",x)},
        contains(c("Estimate.bacon","StdErr.bacon","pValue.bacon.z","fdr.bacon.z"))
      )
    }
  )
  
  results.df <- results.list %>% purrr::reduce(inner_join)
  return(results.df)
}

ABETA.results <- file.list %>% merge.list(csf = "ABETA")
PTAU.results <- file.list %>% merge.list(csf = "PTAU")
TAU.results <- file.list %>% merge.list(csf = "TAU")

results <- left_join(ABETA.results, PTAU.results)
results <- left_join(results, TAU.results)

# Save all significant results of effect size correlation
## Separate by group 
CSF <- c("ABETA", "PTAU", "TAU")
DX <- c("AD", "CN")
CSF.group <- combn(CSF, m = 2, simplify = F)
i <- 1
g <- list()
for(dx in DX){
  for(csf in CSF.group){
    g[[i]] <- ggscatter(
      data = results %>% filter(get(paste0(csf[1], "_", dx, "_pValue.bacon.z")) < 1e-05 | 
                                  get(paste0(csf[2], "_", dx, "_pValue.bacon.z")) < 1e-05 ),
      x = paste0(csf[1], "_", dx,"_Estimate.bacon"),
      y = paste0(csf[2], "_", dx,"_Estimate.bacon"),
      add = "reg.line",
      add.params = list(color = "navy"),
      xlab = paste0(csf[1], " ", dx, " group estimated effect size"),
      ylab = paste0(csf[2], " ", dx, " group estimated effect size"),
      title = paste0("ADNI ", dx, " group ", csf[1], " vs ", csf[2], " \nfor significant CpGs (pValue < 1E-05)"),
      palette = "jco",
      cor.coef = TRUE,
      cor.coeff.args = list(method = "spearman", cor.coef.name = "R", 
                            digits = 3, 
                            label.sep = "\n", color = 'navy')
    )
    i <- i + 1
  }
}
plot <- ggarrange(
  plotlist = g, 
  labels = "AUTO",
  font.label = list(size = 20),
  nrow = 2, ncol = 3
)
ggsave(
  filename = "ADNI_corr_plots_csf_separated_by_ad_and_cn.pdf", 
  plot = plot, 
  path = dir.plot,
  width = 15, height = 10
)

## Separate by CSF 
CSF <- c("ABETA", "PTAU", "TAU")
DX <- c("AD", "CN")
g <- list()
i <- 1
for(csf in CSF){
  g[[i]] <- ggscatter(
    data = results %>% filter(get(paste0(csf,"_", DX[1], "_pValue.bacon.z")) < 1E-05 |
                                get(paste0(csf,"_", DX[2], "_pValue.bacon.z")) < 1E-05),
    x = paste0(csf,"_", DX[1], "_Estimate.bacon"),
    y = paste0(csf,"_", DX[2], "_Estimate.bacon"),
    add = "reg.line",
    add.params = list(color = "navy"),
    xlab = paste0(csf," ", DX[1], " group estimated effect size"),
    ylab = paste0(csf," ", DX[2], " group estimated effect size"),
    title = paste0("ADNI ", csf, " estimated effect size \nfor significant CpGs (pValue < 1E-05)"),
    palette = "jco",
    cor.coef = TRUE,
    cor.coeff.args = list(method = "spearman", cor.coef.name = "R", 
                          digits = 3, 
                          label.sep = "\n", color = 'navy')
  )
  i <- i + 1
}
plot <- ggarrange(
  plotlist = g, 
  labels = "AUTO",
  font.label = list(size = 20),
  nrow = 1, ncol = 3
)
ggsave(
  filename = "ADNI_corr_plots_ad_and_cn_separated_by_csf.pdf", 
  plot = plot, 
  path = dir.plot,
  width = 15, height = 5
)