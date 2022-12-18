# Distinct CSF biomarker-associated DNA methylation in Alzheimer's disease and cognitively normal subjects

Wei Zhang, Juan I. Young, Lissette Gomez, Michael A. Schmidt, David Lukacsovich, Achintya Varma, X. Steven Chen, Eden R. Martin, Lily Wang

### Description

This github repository includes scripts used for the analyses in the above manuscript.

**Method** CSF biomarkers are well-established AD endophenotypes and correlate significantly with neuropathology scores measured on postmortem brain samples. In this work, we performed a comprehensive analysis to identify blood DNA methylation associated with CSF pathological biomarkers for AD in the Alzheimer’s Disease Neuroimaging Initiative (ADNI) cohort. Our study included matched samples of whole blood DNA methylation, CSF Aβ42, phosphorylated tau181 (pTau181), and total tau (tTau) biomarkers data, measured on the same subjects and at the same clinical visits from a total of 202 subjects (123 cognitively normal, 79 AD cases).

**Results** We identified a number of novel associations between blood DNAm and CSF Aβ42, pTau181, and total tau biomarkers, demonstrating that changes in various pathological processes in the CSF are also reflected in the blood epigenome. Overall, the CSF biomarker-associated DNAm is relatively distinct in cognitively normal (CN) and AD subjects, highlighting the importance of analyzing omics data measured on cognitively normal (including preclinical AD) subjects to identify diagnostic biomarkers, and considering disease stages in the development and testing of AD treatment strategies. Moreover, our pathway analysis of CSF biomarker associated-DNAm from CN subjects revealed that biological processes associated with early brain impairment relevant to AD are marked by DNAm in the blood.


### 1. Preprocessing of ADNI DNA methylation data

The ADNI DNA methylation samples were measured with the Illumina HumanMethylation EPIC beadchip. Probes QC and samples QC were conducted. The final sample set was selected by matched CSF measures and clinical information. Outliers were removed by PCA method.

| File                 | Link |
|---------------------|-------------|
| ADNI_preprocessing/ADNI.Rmd  | [Link to the script](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/code/ADNI_preprocessing/ADNI.Rmd) |

### 2. Identification of CSF biomarker-associated CpGs

In this part we analyzed the association between methylation and each CSF biomarker in the cognitive normal (CN) and AD samples separately. A bacon method of adjusting inflation and bias was applied to the final results. In addition, interaction models were fitted to compare the effects of methylation-to-CSF biomarker associations in CN samples and AD samples.

|Result | Link | Link to the script |
|---------|--------------------|-------------|
|AD| [Link](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/results/DNAm-to-CSF-biomarkers/AD/) | [Link to the script](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/code/Linear_model/ADNI_last_visit_CSF_LM_AD_Samples.Rmd)  |
|CN| [Link](https://github.com/TransBioInfoLab/AAD-ATN-biomarkers-and-DNAm/blob/main/results/DNAm-to-CSF-biomarkers/CN/) |[Link to the script](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/code/Linear_model/ADNI_last_visit_CSF_LM_CN_Samples.Rmd)  |
|Interaction||[Link to the script](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/code/Linear_model/ADNI_last_visit_CSF_Interaction_model.Rmd)  |


### 3. Pathway analysis

For pathway analysis, we used a robust rank aggregation method to perform the analysis. This method is implemented in the function `methylRRA` in *methylGSA* R package.

| File | Link |
|---------|--------------------|
| Pathway_analysis/methylGSA.Rmd| [Link to the script](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/code/Pathway_analysis/methylGSA.Rmd)  |



### 4. Integrative methylation-to-gene expression analysis

To evaluate the DNA methylation effect on the gene expression of nearby genes, we analyzed matched gene expression (Affymetrix Human GenomeU 219 array) and DNA methylation (EPIC array) in ADNI data. To associate genes with DNA methylation sites, we used the MethReg R package and considered CpGs located in the promoter regions and distal regions separate.

| File |  Link |
|---------------------|-------------|
| DNAm_RNA/RNA_vs_DNAm.Rmd  | [Link to the script](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/code/DNAm_RNA/RNA_vs_DNAm.Rmd) |


### 5. Correlation and overlap with genetic susceptibility loci

| File |  Link to the script |
|---------------------|-------------|
| mQTLs_analysis/ADbiomarkers_mQTLs_analysis.R | [Link to the script](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/code/mQTLs_analysis/ADbiomarkers_mQTLs_analysis.R) |


# For reproducible research

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

The script for install and library the packages can be found in: Utility/load_package.R ([Link to the script](https://github.com/TransBioInfoLab/AD-ATN-biomarkers-and-DNAm/blob/main/code/Utility/load_package.R))


For ADNIMERGE, download it from https://ida.loni.usc.edu/: Merged ADNI 1/GO/2 Packages for R

For CSF biomarker information, download it from https://ida.loni.usc.edu/: UPENNBIOMK9.CSV

```r
install.packages("/path/to/ADNIMERGE_0.0.1.tar.gz", repos = NULL, type = "source")
```

The platform information are:

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

# Acknowledgement
Data used in the preparation of this article were obtained from the Alzheimer’s Disease Neuroimaging Initiative (ADNI) database (adni.loni.usc.edu). As such, the investigators within the ADNI contributed to the design and implementation of ADNI and/or provided data but did not participate in the analysis or writing of this report. A complete listing of ADNI investigators can be found at: http://adni.loni.usc.edu/wp-content/uploads/how_to_apply/ADNI_Acknowledgement_List.pdf

# References

1. Vasanthakumar, A. et al. Harnessing peripheral DNA methylation differences in the Alzheimer's Disease Neuroimaging Initiative (ADNI) to reveal novel biomarkers of disease. Clin Epigenetics 12, 84 (2020).
