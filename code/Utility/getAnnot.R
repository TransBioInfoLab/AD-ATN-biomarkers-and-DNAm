# ----------------------------------------------------------------------------------
# Source code from methylGSA package
# ----------------------------------------------------------------------------------
getAnnot = function(array.type, group = "all"){
  if(array.type=="450K"){
    FullAnnot = tryCatch({
      getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
    },
    error = function(e){
      stop("IlluminaHumanMethylation450kanno.ilmn12.hg19 needs to
be installed and loaded before running methylglm/methylRRA")
    })
  }else{
    FullAnnot = tryCatch({
      getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
    },
    error = function(e){
      stop("IlluminaHumanMethylationEPICanno.ilm10b4.hg19 needs to
be installed and loaded before running methylglm/methylRRA")
    })
  }
  
  FullAnnot = FullAnnot[,c("Name","UCSC_RefGene_Name","UCSC_RefGene_Group")]
  FullAnnot = FullAnnot[str_length(rownames(FullAnnot))==10,]
  FullAnnot = FullAnnot[!FullAnnot$UCSC_RefGene_Name=="",]
  ## get the first gene in each USCS_RefGene_Name
  temp = vapply(strsplit(FullAnnot$UCSC_RefGene_Name,split=";"),
                '[', 1, FUN.VALUE=character(1))
  FullAnnot$UCSC_RefGene_Name = temp
  ## get the first gene group in each UCSC_RefGene_Group
  temp = vapply(strsplit(FullAnnot$UCSC_RefGene_Group,split=";"),
                '[', 1, FUN.VALUE=character(1))
  FullAnnot$UCSC_RefGene_Group = temp
  
  if(group == "body"){
    FullAnnot = 
      FullAnnot[FullAnnot$UCSC_RefGene_Group%in%c("Body", "1stExon"),]
  }
  
  if(group == "promoter1"){
    FullAnnot = FullAnnot[grepl("TSS",FullAnnot$UCSC_RefGene_Group),]
  }
  
  if(group == "promoter2"){
    FullAnnot = 
      FullAnnot[FullAnnot$UCSC_RefGene_Group%in%c("TSS200", "TSS1500", 
                                                  "1stExon", "5'UTR"),]
  }
  
  return(FullAnnot)
}
