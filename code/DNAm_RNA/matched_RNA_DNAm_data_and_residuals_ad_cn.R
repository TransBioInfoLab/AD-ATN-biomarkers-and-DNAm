#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Article:
# Distinct CSF biomarker-associated DNA methylation in Alzheimer's disease and cognitively normal subjects
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
# Date: 07 Nov. 2022
#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-==-=-=-=-=-=-=-=-=-=-=-=-=-=
library(tidyverse)
library(data.table)
library(matrixStats)
library(SummarizedExperiment)
library(MethReg)
library(xCell)
library(dorothea)
library(readr)
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
# GET ADNI data
#-----------------------------------------------------------------------------
# DNA methylation
cohort <- "ADNI"
dir.base <- "."
dir.data <- file.path(dir.base,"data/",cohort,"/DNA_methylation") 
dir.data.pca <- file.path(dir.data,"/pca_filtering/") 
dir.data.auc <- file.path(dir.base,"data/Aux_Sync/") 
dir.results <- file.path(dir.base, "analysis-results/ADNI/LM/")
dir.results.adj <- file.path(dir.results, "adjusted/")
dir.combp.results <- file.path(dir.base, "analysis-results/comb-p/")
adni.se <- readRDS(file.path(dir.data.pca, "ADNI_QNBMIQ_PCfiltered_min_age_at_visit_65_with_smoke_probes.RDS"))
for(p in grep("dir",ls(),value = T)) dir.create(get(p),recursive = TRUE,showWarnings = FALSE)

# Gene expression: Affymetrix Human Genome U 219 array 
ADNI_Gene_Expression_Profile <- read_csv(
  file.path(
    dir.base, "..", "AD-meta-analysis-blood-samples/datasets/ADNI/data/gene_expression/ADNI_Gene_Expression_Profile/ADNI_Gene_Expression_Profile.csv"),
    skip = 8
  )
  colnames(ADNI_Gene_Expression_Profile) <- gsub("\\...", "X", colnames(ADNI_Gene_Expression_Profile))
  ADNI_Gene_Expression_Profile$X748 <- NULL
  
  # We have our gene expression as below:
  #   ProbeSet      LocusLink Symbol                  X4    X5    X6    X7    X8    X9   X10
  # 1 11715104_s_at LOC92736  OTOP2                 2.15  2.16  2.52  2.28  2.25  2.24  1.99
  # 2 11715105_at   LOC284099 C17ORF78              2.27  2.13  1.96  2.35  2.15  2.06  2.32
  # 3 11715106_x_at NA        CTAGE6 || CTAGE15     2.43  2.27  2.33  2.26  2.33  2.45  2.17
  # The row (3) must be break in to two
  #   ProbeSet      LocusLink Symbol                  X4    X5    X6    X7    X8    X9   X10
  # 1 11715104_s_at LOC92736  OTOP2                 2.15  2.16  2.52  2.28  2.25  2.24  1.99
  # 2 11715105_at   LOC284099 C17ORF78              2.27  2.13  1.96  2.35  2.15  2.06  2.32
  # 3 11715106_x_at NA        CTAGE6                2.43  2.27  2.33  2.26  2.33  2.45  2.17
  # 4 11715106_x_at NA        CTAGE15               2.43  2.27  2.33  2.26  2.33  2.45  2.17
ADNI_Gene_Expression_Profile <- ADNI_Gene_Expression_Profile %>% tidyr::separate_rows("Symbol")

library(hgu219.db)
x <- hgu219ENSEMBL
# Get the entrez gene IDs that are mapped to an Ensembl ID, this is required by MethReg
probe.to.ensg <- as.data.frame(as.list(x) %>% unlist) %>% na.omit
ADNI_Gene_Expression_Profile$ENGS <- probe.to.ensg[ADNI_Gene_Expression_Profile$ProbeSet,]
ADNI_Gene_Expression_Profile <- ADNI_Gene_Expression_Profile[!is.na(ADNI_Gene_Expression_Profile$ENGS),]

nrow(ADNI_Gene_Expression_Profile) # 43890

# dropping genes in bottom 10 percentile for over 80% of the samples
genes.low.expressed.samples.count <- plyr::aaply(
  ADNI_Gene_Expression_Profile[,grep("^X",colnames(ADNI_Gene_Expression_Profile))] ,
  .margins = 2, # for each sample set get all genes
  .fun = function(sample.genes){
    # for each sample, mark the bottom 10% less expressed genes as 1
    sample.genes[[1]]  <= quantile(sample.genes[[1]] , probs = c(.10),type = 3)
  }) %>% colSums()
genes.idx <- which(genes.low.expressed.samples.count > (length(grep("^X",colnames(ADNI_Gene_Expression_Profile))) * 0.8))
ADNI_Gene_Expression_Profile <- ADNI_Gene_Expression_Profile[-c(genes.idx),]

nrow(ADNI_Gene_Expression_Profile) # 41682

# Since we have multiple probes mapping to the same gene in the array
# we will take the median of the same genes
# as suggest in https://www.nature.com/articles/s41598-020-60595-1
expression.matrix <- plyr::aaply(unique(ADNI_Gene_Expression_Profile$ENGS),.margins = 1,.fun = function(gene){
  dat <- ADNI_Gene_Expression_Profile[ADNI_Gene_Expression_Profile$ENGS == gene,grep("^X",colnames(ADNI_Gene_Expression_Profile))]
  colMedians(dat %>% as.matrix)
},.progress = "time")
rownames(expression.matrix) <- unique(ADNI_Gene_Expression_Profile$ENGS)

ADNI_Gene_Expression_Metadata <- read_csv(
  file.path(dir.base, "../AD-meta-analysis-blood-samples/datasets/ADNI/data/gene_expression/ADNI_Gene_Expression_Profile/ADNI_Gene_Expression_Profile.csv"),
  skip = 0,col_names = FALSE,n_max = 7
)
ADNI_Gene_Expression_Metadata$X748 <- NULL

gene.exp.IDs <- apply(
  ADNI_Gene_Expression_Metadata[,4:ncol(ADNI_Gene_Expression_Metadata)],MARGIN = 2,
  FUN = function(col) {
    paste0(
      stringr::str_extract(pattern = "[0-9]*$",string = col[3]),
      "_", col[1],"_",
      col[2]
    )
  }
)

DXSUM <- readr::read_csv(
  file.path(dir.base,"../AD-meta-analysis-blood-samples/datasets/ADNI/data/Clinical/DXSUM_PDXCONV_ADNIALL_downloaded_2-19-2021.csv")
)
DXSUM <- DXSUM[match(gene.exp.IDs,paste0(stringr::str_extract(pattern = "[0-9]*$",DXSUM$PTID),"_",DXSUM$Phase,"_",DXSUM$VISCODE)),]
gene.exp.IDs <- paste0(DXSUM$RID,"_",DXSUM$Phase,"_",DXSUM$VISCODE2)
colnames(expression.matrix) <- gene.exp.IDs
expression.matrix <- expression.matrix[,colnames(expression.matrix) != "NA_NA_NA"]

#-----------------------------------------------------------------------------
# get matched DNAm and Gene expression
#-----------------------------------------------------------------------------
dnam.IDs <- paste0(adni.se$RID,"_",adni.se$COLPROT,"_",adni.se$VISCODE)
table(gene.exp.IDs %in% dnam.IDs)
# FALSE  TRUE
#  262   482
common.ids <- base::intersect(dnam.IDs %>% as.character,gene.exp.IDs %>% as.character)
adni.se <- adni.se[,paste0(adni.se$RID,"_",adni.se$COLPROT,"_",adni.se$VISCODE) %in% common.ids]
adni.se <- adni.se[,match(common.ids,paste0(adni.se$RID,"_",adni.se$COLPROT,"_",adni.se$VISCODE))]
expression.matrix <- expression.matrix[,colnames(expression.matrix) %in% common.ids]
expression.matrix <- expression.matrix[,match(common.ids,colnames(expression.matrix))]

table(colnames(expression.matrix) == paste0(adni.se$RID,"_",adni.se$COLPROT,"_",adni.se$VISCODE))

save(
  adni.se,
  expression.matrix,
  ADNI_Gene_Expression_Metadata,
  file = file.path(dir.data.auc,"ADNI_all_samples_matched_rna_dnam.rda")
)
#-----------------------------------------------------------------------------
# get AD and CN samples only
#-----------------------------------------------------------------------------
# AD samples
#-----------------------------------------------------------------------------
load(file.path(dir.data.auc,"ADNI_all_samples_matched_rna_dnam.rda"))
pheno <- colData(adni.se)
pheno.ad <- pheno %>% as.data.frame %>% filter(DX %in% c("Dementia"))
adni.se.ad <- adni.se[,rownames(pheno.ad)]
expression.matrix.ad <- expression.matrix[,paste0(adni.se.ad$RID,"_",adni.se.ad$COLPROT,"_",adni.se.ad$VISCODE)]

gene.exp.IDs <- apply(
  ADNI_Gene_Expression_Metadata[,4:ncol(ADNI_Gene_Expression_Metadata)],MARGIN = 2,
  FUN = function(col) {
    paste0(
      stringr::str_extract(pattern = "[0-9]*$",string = col[3]),
      "_", col[1],"_",
      col[2]
    )
  }
)

DXSUM <- readr::read_csv(
  file.path(dir.base,"../AD-meta-analysis-blood-samples/datasets/ADNI/data/Clinical/DXSUM_PDXCONV_ADNIALL_downloaded_2-19-2021.csv")
)
DXSUM <- DXSUM[match(gene.exp.IDs,paste0(stringr::str_extract(pattern = "[0-9]*$",DXSUM$PTID),"_",DXSUM$Phase,"_",DXSUM$VISCODE)),]
gene.exp.IDs <- paste0(DXSUM$RID,"_",DXSUM$Phase,"_",DXSUM$VISCODE2)

colnames(ADNI_Gene_Expression_Metadata)[-c(1:3)] <- gene.exp.IDs
ADNI_Gene_Expression_Metadata.ad <- cbind(ADNI_Gene_Expression_Metadata[,c(1:3)], ADNI_Gene_Expression_Metadata[,colnames(expression.matrix.ad)])
#-----------------------------------------------------------------------------
# CN samples
#-----------------------------------------------------------------------------
pheno.cn <- pheno %>% as.data.frame %>% filter(DX %in% c("CN"))
adni.se.cn <- adni.se[,rownames(pheno.cn)]
expression.matrix.cn <- expression.matrix[,paste0(adni.se.cn$RID,"_",adni.se.cn$COLPROT,"_",adni.se.cn$VISCODE)]

ADNI_Gene_Expression_Metadata.cn <- cbind(ADNI_Gene_Expression_Metadata[,c(1:3)], ADNI_Gene_Expression_Metadata[,colnames(expression.matrix.cn)])

save(
  adni.se.ad,
  adni.se.cn,
  expression.matrix.ad,
  expression.matrix.cn,
  ADNI_Gene_Expression_Metadata.ad,
  ADNI_Gene_Expression_Metadata.cn,
  file = file.path(dir.data.auc,"ADNI_ad_cn_matched_rna_dnam.rda")
)
#-----------------------------------------------------------------------------
# get residuals for sig cpgs
#-----------------------------------------------------------------------------
load(file.path(dir.data.auc,"ADNI_ad_cn_matched_rna_dnam.rda"))
## Load regression results
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

single.sig.cpgs <- file.list.lm %>% map(~filter(., pValue.bacon.z < 1e-05)) %>%
  rbindlist(.,idcol = "ID") %>% mutate(ID = gsub("ADNI_","",ID)) %>%
  separate(ID, into = c("Sample", "CSF"), sep = "_",  remove = T)
## Load comb-p DMR results
load(file.path(dir.combp.results,"ADNI_all_DMR_annotated_filtered_zsidakp_zp_same_direction_results.rda"))
dmr.sig.cpgs <- dmr.filter.list %>%
  map(~separate_rows(.,cpgs_in_region, sep = ",")) %>%
  map(~dplyr::mutate(.,cpg = cpgs_in_region, .keep = "unused")) %>%
  rbindlist(.,idcol = "ID") %>% mutate(ID = gsub("ADNI_","",ID)) %>%
  separate(ID, into = c("Sample", "CSF"), sep = "_",  remove = T)
### Get median
region <- dmr.sig.cpgs$region %>% unique()
beta <- cbind(assay(adni.se.ad), assay(adni.se.cn))
dmr.median <- plyr::ldply(
  region,
  .fun = function(r){
    cpg <- dmr.sig.cpgs %>% filter(region == r) %>% pull(cpg)
    m <- colMedians(beta[cpg,])

    return(m)
  }
)
colnames(dmr.median) <- colnames(beta)
rownames(dmr.median) <- region
#-----------------------------------------------------------------------------
# AD samples
#-----------------------------------------------------------------------------
single.sig.cpgs.ad <- single.sig.cpgs %>% filter(Sample == "AD")
dmr.ad <- dmr.sig.cpgs %>% filter(Sample == "AD") %>% pull(region)
sig.cpgs.dmr.ad <- c(single.sig.cpgs.ad$cpg, dmr.ad) %>% unique
length(sig.cpgs.dmr.ad) #228
se.selected.cpgs <- adni.se.ad[rownames(adni.se.ad) %in% single.sig.cpgs.ad$cpg,]
selected.dmr.median <- dmr.median[rownames(dmr.median) %in% dmr.ad, colnames(adni.se.ad)] %>% as.matrix()
colnames(se.selected.cpgs) <- paste0(se.selected.cpgs$RID,"_",se.selected.cpgs$COLPROT,"_",se.selected.cpgs$VISCODE)

metadata.dnam <- colData(se.selected.cpgs)[,c("RID","DX","CD8T","Mono","Neutro","CD4T","NK","B","Eosino","age_at_visit","PTGENDER","PlateNumber", "Smok_hist", "PTEDUCAT", "APOE4")]
metadata.dnam$granulocytes <- metadata.dnam$Neutro +metadata.dnam$Eosino
metadata.dnam$PlateNumber <- as.factor(metadata.dnam$PlateNumber)

Mvalues.cpgs <- log2(assay(se.selected.cpgs) / (1 - assay(se.selected.cpgs)))
Mvalues.dmrs <- log2(selected.dmr.median/(1-selected.dmr.median))

residuals.matched.met <- MethReg::get_residuals(
  data = rbind(Mvalues.cpgs,Mvalues.dmrs), # m-values
  metadata.samples = metadata.dnam[,c("Mono","CD4T","NK","B","granulocytes", "age_at_visit","PTGENDER","PlateNumber", "Smok_hist", "PTEDUCAT", "APOE4")],
  cores = 4
)

Affy_Plate <- ADNI_Gene_Expression_Metadata.ad[7,-c(1:3)] %>% as.numeric() %>% as.factor()
Affy_Plate <- Affy_Plate[
  match(
    gsub(" ","0",formatC(stringr::str_extract(pattern = "^[0-9]*",colnames(expression.matrix.ad)))),
    stringr::str_extract(pattern = "[0-9]*$",ADNI_Gene_Expression_Metadata.ad[3,-c(1:3)] %>% as.character()))
]

RIN <- ADNI_Gene_Expression_Metadata.ad[6,-c(1:3)] %>% as.numeric()
RIN <- RIN[
  match(
    gsub(" ","0",formatC(stringr::str_extract(pattern = "^[0-9]*",colnames(expression.matrix.ad)))),
    stringr::str_extract(pattern = "[0-9]*$",ADNI_Gene_Expression_Metadata.ad[3,-c(1:3)] %>% as.character()))
]

metadata.exp <- colData(se.selected.cpgs)[,c("age_at_visit", "PTGENDER", "Smok_hist", "PTEDUCAT", "APOE4")]
metadata.exp$Affy_Plate <- Affy_Plate
metadata.exp$RIN <- RIN

aux <- expression.matrix.ad
rownames(aux) <- MethReg:::map_ensg_to_symbol(rownames(aux))
xcell <- xCellAnalysis(aux)
xcell <- xcell[c("B-cells", "NKT", "CD4+ T-cells", "Eosinophils", "Monocytes", "Neutrophils"),]

metadata.exp <- cbind(metadata.exp,t(xcell))
colnames(metadata.exp)[8] <- "B_cells"
colnames(metadata.exp)[10] <- "CD4_T_cells"
metadata.exp$granulocytes <- metadata.exp$Neutrophils + metadata.exp$Eosinophils
metadata.exp$Neutrophils <- NULL
metadata.exp$Eosinophils <- NULL

residuals.matched.exp <- MethReg::get_residuals(
  data = log2(expression.matrix.ad),
  metadata.samples = metadata.exp[,],
  cores = 4
)

all(rownames(metadata.dnam) == rownames(metadata.exp))
all(colnames(residuals.matched.exp) == colnames(residuals.matched.met))

#-----------------------------------------------------------------------------
# get TF activity
#-----------------------------------------------------------------------------
regulons.dorothea <- dorothea::dorothea_hs
rnaseq.tf.es <- MethReg::get_tf_ES(
  exp = residuals.matched.exp,
  regulons = regulons.dorothea
)

#-----------------------------------------------------------------------------
# Save data
#-----------------------------------------------------------------------------

save(
  metadata.dnam,
  se.selected.cpgs,
  selected.dmr.median,
  residuals.matched.met,
  metadata.exp,
  rnaseq.tf.es,
  residuals.matched.exp,
  xcell,
  file = file.path(dir.data.auc,"ADNI_AD_matched_rna_dnam_residuals.rda")
)
#-----------------------------------------------------------------------------
# CN samples
#-----------------------------------------------------------------------------
single.sig.cpgs.cn <- single.sig.cpgs %>% filter(Sample == "CN")
dmr.cn <- dmr.sig.cpgs %>% filter(Sample == "CN") %>% pull(region)
sig.cpgs.dmr.cn <- c(single.sig.cpgs.cn$cpg, dmr.cn) %>% unique
length(sig.cpgs.dmr.cn) #163
se.selected.cpgs <- adni.se.cn[rownames(adni.se.cn) %in% single.sig.cpgs.cn$cpg,]
selected.dmr.median <- dmr.median[rownames(dmr.median) %in% dmr.cn, colnames(adni.se.cn)] %>% as.matrix()

colnames(se.selected.cpgs) <- paste0(se.selected.cpgs$RID,"_",se.selected.cpgs$COLPROT,"_",se.selected.cpgs$VISCODE)
metadata.dnam <- colData(se.selected.cpgs)[,c("RID","DX","CD8T","Mono","Neutro","CD4T","NK","B","Eosino","age_at_visit","PTGENDER","PlateNumber", "Smok_hist", "PTEDUCAT", "APOE4")]
metadata.dnam$granulocytes <- metadata.dnam$Neutro +metadata.dnam$Eosino
metadata.dnam$PlateNumber <- as.factor(metadata.dnam$PlateNumber)

Mvalues.cpgs <- log2(assay(se.selected.cpgs) / (1 - assay(se.selected.cpgs)))
Mvalues.dmrs <- log2(selected.dmr.median/(1-selected.dmr.median))

residuals.matched.met <- MethReg::get_residuals(
  data = rbind(Mvalues.cpgs, Mvalues.dmrs), # m-values
  metadata.samples = metadata.dnam[,c("Mono","CD4T","NK","B","granulocytes", "age_at_visit","PTGENDER","PlateNumber", "Smok_hist", "PTEDUCAT", "APOE4")],
  cores = 4
)

Affy_Plate <- ADNI_Gene_Expression_Metadata.cn[7,-c(1:3)] %>% as.numeric() %>% as.factor()
Affy_Plate <- Affy_Plate[
  match(
    gsub(" ","0",formatC(stringr::str_extract(pattern = "^[0-9]*",colnames(expression.matrix.cn)))),
    stringr::str_extract(pattern = "[0-9]*$",ADNI_Gene_Expression_Metadata.cn[3,-c(1:3)] %>% as.character()))
]

RIN <- ADNI_Gene_Expression_Metadata.cn[6,-c(1:3)] %>% as.numeric()
RIN <- RIN[
  match(
    gsub(" ","0",formatC(stringr::str_extract(pattern = "^[0-9]*",colnames(expression.matrix.cn)))),
    stringr::str_extract(pattern = "[0-9]*$",ADNI_Gene_Expression_Metadata.cn[3,-c(1:3)] %>% as.character()))
]

metadata.exp <- colData(se.selected.cpgs)[,c("age_at_visit", "PTGENDER", "Smok_hist", "PTEDUCAT", "APOE4")]
metadata.exp$Affy_Plate <- Affy_Plate
metadata.exp$RIN <- RIN

aux <- expression.matrix.cn
rownames(aux) <- MethReg:::map_ensg_to_symbol(rownames(aux))
xcell <- xCellAnalysis(aux)
xcell <- xcell[c("B-cells", "NKT", "CD4+ T-cells", "Eosinophils", "Monocytes", "Neutrophils"),]

metadata.exp <- cbind(metadata.exp,t(xcell))
colnames(metadata.exp)[8] <- "B_cells"
colnames(metadata.exp)[10] <- "CD4_T_cells"
metadata.exp$granulocytes <- metadata.exp$Neutrophils + metadata.exp$Eosinophils
metadata.exp$Neutrophils <- NULL
metadata.exp$Eosinophils <- NULL

residuals.matched.exp <- MethReg::get_residuals(
  data = log2(expression.matrix.cn),
  metadata.samples = metadata.exp[,],
  cores = 4
)

all(rownames(metadata.dnam) == rownames(metadata.exp))
all(colnames(residuals.matched.exp) == colnames(residuals.matched.met))

#-----------------------------------------------------------------------------
# get TF activity
#-----------------------------------------------------------------------------
regulons.dorothea <- dorothea::dorothea_hs
rnaseq.tf.es <- MethReg::get_tf_ES(
  exp = residuals.matched.exp,
  regulons = regulons.dorothea
)

#-----------------------------------------------------------------------------
# Save data
#-----------------------------------------------------------------------------

save(
  metadata.dnam,
  se.selected.cpgs,
  selected.dmr.median,
  residuals.matched.met,
  metadata.exp,
  rnaseq.tf.es,
  residuals.matched.exp,
  xcell,
  file = file.path(dir.data.auc,"ADNI_CN_matched_rna_dnam_residuals.rda")
)


