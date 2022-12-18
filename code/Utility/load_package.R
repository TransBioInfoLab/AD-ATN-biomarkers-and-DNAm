### For reproducible research, the following R packages are required:

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

