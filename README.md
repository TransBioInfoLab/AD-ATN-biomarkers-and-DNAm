# Distinct CSF biomarker-associated DNA methylation in Alzheimer's disease and cognitively normal subjects

Wei Zhang, Juan I. Young, Lissette Gomez, Michael A. Schmidt, David Lukacsovich, Achintya Varma, X. Steven Chen, Eden R. Martin, Lily Wang

## Description

This github repository includes scripts used for the analyses in the above manuscript.

CSF biomarkers are well-established AD endophenotypes and correlate significantly with neuropathology scores measured on postmortem brain samples. In this work, we performed a comprehensive analysis to identify blood DNA methylation associated with CSF pathological biomarkers for AD in the Alzheimer’s Disease Neuroimaging Initiative (ADNI) cohort. Our study included matched samples of whole blood DNA methylation, CSF Aβ42, phosphorylated tau181 (pTau181), and total tau (tTau) biomarkers data, measured on the same subjects and at the same clinical visits from a total of 202 subjects (123 cognitively normal, 79 AD cases).

We identified a number of novel associations between blood DNAm and CSF Aβ42, pTau181, and total tau biomarkers, demonstrating that changes in various pathological processes in the CSF are also reflected in the blood epigenome. Overall, the CSF biomarker-associated DNAm is relatively distinct in cognitively normal (CN) and AD subjects, highlighting the importance of analyzing omics data measured on cognitively normal (including preclinical AD) subjects to identify diagnostic biomarkers, and considering disease stages in the development and testing of AD treatment strategies. Moreover, our pathway analysis of CSF biomarker associated-DNAm from CN subjects revealed that biological processes associated with early brain impairment relevant to AD are marked by DNAm in the blood.


### 1. Preprocessing of ADNI DNA methylation data

The DNA methylation samples were measured with the Illumina HumanMethylation EPIC beadchip. Quality control of both probes and samples were performed. 

| File                 | Link |
|---------------------|-------------|
| ADNI_preprocessing/ADNI.Rmd  | [Link to the script](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/code/ADNI_preprocessing/ADNI.Rmd) |

### 2. Identification of CSF biomarker-associated CpGs

We studied the association between DNA methylation and CSF biomarkers for AD in the cognitive normal (CN) and AD subjects separately. The bacon method was used to correct genomic inflation and bias. In addition, interaction models were used to compare the methylation-to-CSF biomarker associations in CN samples and AD samples.

|File | Link | Link to the script |
|---------|--------------------|-------------|
|AD| [Link to results](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/results/DNAm-to-CSF-biomarkers/AD/) | [Link to the script](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/code/Linear_model/ADNI_last_visit_CSF_LM_AD_Samples.Rmd)  |
|CN| [Link to results](https://github.com/TransBioInfoLab/AAD-ATN-biomarkers-and-DNAm/blob/main/results/DNAm-to-CSF-biomarkers/CN/) |[Link to the script](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/code/Linear_model/ADNI_last_visit_CSF_LM_CN_Samples.Rmd)  |
|Interaction||[Link to the script](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/code/Linear_model/ADNI_last_visit_CSF_Interaction_model.Rmd)  |


### 3. Pathway analysis

We used the function `methylRRA` in *methylGSA* R package for this analysis, to identify biological pathways enriched with CSF biomarker-associated DNA methylation. 

| File | Link |
|---------|--------------------|
| Pathway_analysis/methylGSA.Rmd| [Link to the script](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/code/Pathway_analysis/methylGSA.Rmd)  |



### 4. Integrative methylation-to-gene expression analysis

To prioritize significant DNAm with downstream functional effects, we correlated DNA methylation levels of the significant DMRs or CpGs with the expression levels of genes found in their vicinity, using matched DNAm and gene expression samples generated from 263 independent subjects (84 AD cases and 179 CN) in the ADNI cohort.
To associate genes with DNA methylation sites, we used the *MethReg* R package and considered CpGs located in the promoter and distal regions separately.

| File |  Link |
|---------------------|-------------|
| DNAm_RNA/RNA_vs_DNAm.Rmd  | [Link to the script](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/code/DNAm_RNA/RNA_vs_DNAm.Rmd) |


### 5. Correlation and overlap with genetic susceptibility loci

To identify methylation quantitative trait loci (mQTLs) for the significant DMRs and CpGs, we performed look-up analyses using the GoDMC database for mQTLs. 

| File |  Link to the script |
|---------------------|-------------|
| mQTLs_analysis/ADbiomarkers_mQTLs_analysis.R | [Link to the script](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/code/mQTLs_analysis/ADbiomarkers_mQTLs_analysis.R) |

### 6. Whole genome analysis results (CN samples)

| File |  Link to the script |
|---------------------|-------------|
| DNAm associated with abeta | [Link to results](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/results/DNAm-to-CSF-biomarkers/CN/ADNI_Last_Visit_CSF_ABETA_with_covariates_LM_CN_results_with_inflation.csv) |
| DNAm associated with pTau | [Link to results](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/results/DNAm-to-CSF-biomarkers/CN/ADNI_Last_Visit_CSF_PTAU_with_covariates_LM_CN_results_with_inflation.csv) |
| DNAm associated with tTau | [Link to results](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/results/DNAm-to-CSF-biomarkers/CN/ADNI_Last_Visit_CSF_TAU_with_covariates_LM_CN_results_with_inflation.csv) |

### 7. Whole genome analysis results (AD samples)

| File |  Link to the script |
|---------------------|-------------|
| DNAm associated with abeta | [Link to results](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/results/DNAm-to-CSF-biomarkers/AD/ADNI_Last_Visit_CSF_ABETA_with_covariates_LM_AD_results_with_inflation.csv) |
| DNAm associated with pTau | [Link to results](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/results/DNAm-to-CSF-biomarkers/AD/ADNI_Last_Visit_CSF_PTAU_with_covariates_LM_AD_results_with_inflation.csv) |
| DNAm associated with tTau | [Link to results](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/results/DNAm-to-CSF-biomarkers/AD/ADNI_Last_Visit_CSF_TAU_with_covariates_LM_AD_results_with_inflation.csv) |


## For reproducible research

The following R packages are required:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install(version = "3.15")

list.of.packages <- c(
  "bacon",
  "coMethDMR",
  "devtools",
  "DMRcate",                                      
  "doParallel",  
  "dorothea",
  "dplyr",                                        
  "DT",                                           
  "EpiDISH",                                      
  "ExperimentHub",                                
  "fgsea",                                        
  "GenomicRanges",                                
  "GEOquery",                                     
  "ggpubr",                                       
  "ggrepel",
  "ggvenn",
  "gridExtra",                                    
  "gt",
  "GWASTools",  
  "hgu219.db",
  "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
  "lubridate",                                    
  "lumi",                                         
  "meta",                                         
  "metap",
  "methylGSA",
  "MethReg",                                      
  "minfi",                                        
  "missMethyl",                                   
  "mygene",                                       
  "plyr",                                         
  "readr",                                        
  "readxl",                                       
  "ReMapEnrich",
  "rGREAT",
  "RPMM",                                         
  "sm",                                           
  "stats",                                        
  "SummarizedExperiment",                         
  "tidyverse",                                        
  "wateRmelon",                                   
  "writexl"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)

devtools::install_github("igordot/msigdbr")
devtools::install_github('dviraran/xCell')
```

The script for installing the packages can be found at: Utility/load_package.R ([Link to the script](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/code/Utility/load_package.R))


The ADNIMERGE R package can be downloaded from https://ida.loni.usc.edu/: Merged ADNI 1/GO/2 Packages for R

The CSF biomarker data can be accessed from https://ida.loni.usc.edu/: file "UPENNBIOMK9.CSV"

```r
install.packages("/path/to/ADNIMERGE_0.0.1.tar.gz", repos = NULL, type = "source")
```

The platform information is:

```r
version  R version 4.2.1 (2022-06-23)
os       macOS Ventura 13.0
system   x86_64, darwin17.0
ui       RStudio
language (EN)
collate  en_US.UTF-8
ctype    en_US.UTF-8
tz       America/New_York
date     2022-12-16
```

## Acknowledgement
Data used in the preparation of this article were obtained from the Alzheimer’s Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu). As such, the investigators within the ADNI contributed to the design and implementation of ADNI and/or provided data but did not participate in the analysis or writing of this report. A complete listing of ADNI investigators can be found at: http://adni.loni.usc.edu/wp-content/uploads/how_to_apply/ADNI_Acknowledgement_List.pdf


