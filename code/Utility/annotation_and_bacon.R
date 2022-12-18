# --------------------------------------------------------------------------------------------------- #
          ############### Function for annotation and inflation adjusted ################
# --------------------------------------------------------------------------------------------------- #

dir.base <- "."
dir.data.aux <- file.path(dir.base,"../AD-meta-analysis-blood-samples/datasets/Aux_Sync/") 

load(file.path(dir.data.aux,"great_EPIC_array_annotation.rda"))
load(file.path(dir.data.aux,"E073_15_coreMarks_segments.rda"))
load(file.path(dir.data.aux,"meta_analysis_cpgs.rda"))

library(SummarizedExperiment)
library(tidyverse)
library(bacon)
library(GWASTools)

# --------------------------------------------------------------------------------------------------- #
# Annotation
# --------------------------------------------------------------------------------------------------- #

annotate_results <- function(result, data.se){
  pval_type <- grep("Pr(>|W|)|Pr(>|t|)|Pr",colnames(result), value = T)
  result$fdr <- p.adjust(result[[pval_type]],method = "fdr")
  result <- cbind(
    result,
    rowRanges(data.se)[result$cpg] %>% as.data.frame() %>% dplyr::select(1:4)
  )
  result$Islands.UCSC.Relation_to_Island <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Islands.UCSC[result$cpg,"Relation_to_Island"]
  result$UCSC_RefGene_Name <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other[result$cpg,"UCSC_RefGene_Name"]       
  result$UCSC_RefGene_Group <- IlluminaHumanMethylationEPICanno.ilm10b4.hg19::Other[result$cpg,"UCSC_RefGene_Group"]     
  
  result <- dplyr::left_join(result, great,by = c("seqnames","start","end","cpg"))
  
  hits <- findOverlaps(rowRanges(data.se)[result$cpg],ChmmModels.gr) %>% as.data.frame()
  hits$state <- ChmmModels.gr$state[hits$subjectHits]
  hits$cpg <- result$cpg[hits$queryHits]
  result <- merge(result, hits[,c("state","cpg")], all.x = TRUE,sort = FALSE, by = "cpg")
  
  # sig.in.brain: 1 if it overlaps with the 3751 CpGs or 119 DMRs from our brain samples meta-analysis
  result$sig.in.brain <- result$cpg %in% brain.meta.analysis.cpgs
  return(result)
}

# --------------------------------------------------------------------------------------------------- #
# Bacon inflation
# --------------------------------------------------------------------------------------------------- #

estimation_of_inflation <- function(data){
  
  ### 1. Compute genomic inflation factor before bacon adjustment
  stats.type <- grep("z value|t value|Wald|t.value",colnames(data), value = T)
  
  data <- data %>% mutate(
    chisq = get(stats.type)^2
  ) 
  
  # inflation factor - last term is median from chisq distrn with 1 df  
  inflationFactor <- median(data$chisq,na.rm = TRUE) / qchisq(0.5, 1)
  print("lambda")
  print(inflationFactor)
  
  se <- data[[grep("Std", colnames(data), value = T)]]
  
  ### 2. bacon analysis
  bc <- bacon(
    teststatistics = NULL,
    effectsizes =  data$Estimate,
    standarderrors = se,
    na.exclude = TRUE,
    verbose = F
  )
  #  posteriors(bc)
  # inflation factor
  print("lambda.bacon")
  print(inflation(bc))
  # bias
  print("estimate bias")
  print(bias(bc))
  print("estimates")
  print(estimates(bc))
  
  ### 3. Create final dataset
  data.with.inflation <- data.frame(
    data,
    Estimate.bacon = bacon::es(bc),
    StdErr.bacon = bacon::se(bc),
    pValue.bacon = pval(bc),
    fdr.bacon = p.adjust(pval(bc), method = "fdr"),
    stringsAsFactors = FALSE
  )
  
  print("o After bacon correction")
  print("Conventional lambda")
  lambda.con <- median((data.with.inflation$Estimate.bacon/data.with.inflation$StdErr.bacon) ^ 2,na.rm = TRUE)/qchisq(0.5, 1)
  print(lambda.con)
  
  # percent_null <- trunc ( estimates(bc)[1]*100, digits = 0)
  # percent_1  <- trunc ( estimates(bc)[2]*100, digits = 0 )
  # percent_2  <- 100 - percent_null - percent_1  
  bc2 <- bacon(
    teststatistics = NULL,
    effectsizes =  data.with.inflation$Estimate.bacon,
    standarderrors = data.with.inflation$StdErr.bacon,
    na.exclude = TRUE,
    priors = list(
      sigma = list(alpha = 1.28,  beta = 0.36), 
      mu = list(lambda = c(0, 3, -3), tau = c(1000, 100, 100)), 
      epsilon = list(gamma = c(99, .5, .5)))
  )
  
  print("inflation")
  print(inflation(bc2))
  print("estimates")
  print(estimates(bc2))
  
  data.with.inflation <- data.with.inflation %>% dplyr::select(-chisq)
  
  inflation.stat <- data.frame(
    "Inflation.org" = inflationFactor,
    "Inflation.bacon" = inflation(bc),
    "Bias.bacon" = bias(bc),
    "Inflation.after.correction" = lambda.con,
    "Inflation.bacon.after.correction" = inflation(bc2),
    "Bias.bacon.after.correction" = bias(bc2)
  )
  return(
    list(
      "data.with.inflation" = data.with.inflation,
      "inflation.stat" = inflation.stat
    )
  )
}

zscore_of_inflation <- function(data){
  
  ### 1. Convert t scores to z scores
  z.scores <- data[["t.value"]]
  
  ### 2. bacon analysis
  bc <- bacon(
    teststatistics = z.scores,
    na.exclude = TRUE,
    verbose = F
  )
  # inflation factor
  print("lambda.bacon")
  print(inflation(bc))
  # bias
  print("estimate bias")
  print(bias(bc))
  print("estimates")
  print(estimates(bc))
  
  ### 3. Create final dataset
  data.with.inflation <- data %>% mutate(
    zScore.bacon = tstat(bc)[,1],
    pValue.bacon.z = pval(bc)[,1],
    fdr.bacon.z = p.adjust(pval(bc), method = "fdr"),
  ) %>% mutate(z.value = z.scores,.after = t.value)
  
  print("o After bacon correction")
  print("Conventional lambda")
  lambda.con <- median((data.with.inflation$zScore.bacon) ^ 2,na.rm = TRUE)/qchisq(0.5, 1)
  print(lambda.con)
  
  percent_null <- trunc ( estimates(bc)[1]*100, digits = 0)
  percent_1  <- trunc ( estimates(bc)[2]*100, digits = 0 )
  percent_2  <- 100 - percent_null - percent_1  
  bc2 <- bacon(
    teststatistics = data.with.inflation$zScore.bacon,
    na.exclude = TRUE,
    priors = list(
      sigma = list(alpha = 1.28,  beta = 0.36), 
      mu = list(lambda = c(0, 3, -3), tau = c(1000, 100, 100)), 
      epsilon = list(gamma = c(90, 5, 5)))
  )
  
  print("inflation")
  print(inflation(bc2))
  print("estimates")
  print(estimates(bc2))
  
  inflation.stat <- data.frame(
    "Inflation.bacon.z" = inflation(bc),
    "Bias.bacon.z" = bias(bc),
    "Inflation.after.correction.z" = lambda.con,
    "Inflation.bacon.after.correction.z" = inflation(bc2),
    "Bias.bacon.after.correction.z" = bias(bc2)
  )
  
  return(
    list(
      "data.with.inflation" = data.with.inflation,
      "bacon.obj" = bc,
      "inflation.stat" = inflation.stat
    )
  )
}