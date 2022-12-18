setwd("M:/AD/analysis/projects/Sex_strat/Aim3/Meta-analysis/Lissette/bloodMeta_mQTLs/ADbiomarkers_mqtls")


cnCpGs <- read.csv("cnCpGs.csv")
colnames(cnCpGs) <- "cpg"

adCpGs <- read.csv("adCpGs.csv")
colnames(adCpGs) <- "cpg"

#blood <- read.csv("../assoc_meta_all.csv.gz")
#m450k <- read.csv("../450kmethyArray.csv")
#blood$sig <- ifelse(blood$cistrans == "TRUE" & blood$pval_are < 1e-8, 1, 
#                    ifelse(blood$cistrans == "FALSE" & blood$pval_are < 1e-14, 1, 0))
#bloodsig <- blood[which(blood$sig == 1),]
#bloodsig_cpg450k <- bloodsig[which(bloodsig$cpg%in%m450k$Name),]
#write.csv(bloodsig, "assoc_meta_all_sigBloodCpGsALL.csv")
#write.csv(bloodsig_cpg450k, "bloodsig_cpg450k.csv")
bloodsig <- read.csv("../sexStrat/sexStrat_mqtls_071522/assoc_meta_all_sigBloodCpGsALL.csv")


cnCpGs_bloodmQTLs <- merge(cnCpGs, bloodsig, by = "cpg")
write.csv(cnCpGs_bloodmQTLs, "cnCpGs_bloodmQTLs.csv")

adCpGs_bloodmQTLs <- merge(adCpGs, bloodsig, by = "cpg")
write.csv(adCpGs_bloodmQTLs, "adCpGs_bloodmQTLs.csv")

#cnCpGs_bloodmQTLs <- read.csv("cnCpGs_bloodmQTLs.csv")
#adCpGs_bloodmQTLs <- read.csv("adCpGs_bloodmQTLs.csv")


#Annotation

library(dplyr)
library(ELMER.data)
library(sesameData)

probes.info <- sesameDataGet("HM450.hg19.manifest") %>% as.data.frame()
probes.info <- probes.info[, c("seqnames", "start", "end")]
colnames(probes.info)[1] <- "chr"
probes.info$chr <- as.character(probes.info$chr)


UCSC_annot <- read.csv("../infinium-methylationepic-v-1-0-b5-manifest-file.csv")
UCSC_annot <- UCSC_annot[
  , c("IlmnID", "UCSC_RefGene_Group", "UCSC_RefGene_Accession",
      "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island")]
colnames(UCSC_annot) <- c(
  "cpg", "UCSC_RefGene_Group", "UCSC_RefGene_Accession",
  "UCSC_RefGene_Name", "Relation_to_Island")
UCSC_annot$cpg <- as.character(UCSC_annot$cpg)
UCSC_annot$UCSC_RefGene_Group <- as.character(UCSC_annot$UCSC_RefGene_Group)
UCSC_annot$UCSC_RefGene_Accession <- as.character(UCSC_annot$UCSC_RefGene_Accession)
UCSC_annot$UCSC_RefGene_Name <- as.character(UCSC_annot$UCSC_RefGene_Name)
UCSC_annot$Relation_to_Island <- as.character(UCSC_annot$Relation_to_Island)


annotateCpG <- function(data){
  
  ### add chrom and position
  dat <- merge(
    probes.info, data, by.x = "row.names", by.y = "cpg", sort = FALSE)
  colnames(dat)[1] <- "cpg"
  dat$cpg <- as.character(dat$cpg)
  
  ### add UCSC annotation
  dat <- merge(UCSC_annot, dat, by = "cpg", sort = FALSE)
  
}

cnCpGs_bloodmQTLs_annot <- annotateCpG(data = cnCpGs_bloodmQTLs)
write.csv(cnCpGs_bloodmQTLs_annot, "cnCpGs_bloodmQTLs_annot.csv")

adCpGs_bloodmQTLs_anoot <- annotateCpG(data = adCpGs_bloodmQTLs)
write.csv(adCpGs_bloodmQTLs_anoot, "adCpGs_bloodmQTLs_annot.csv")


#Kunkle results
gwas <- read.table("../Kunkle_etal_Stage1_results.txt", header = TRUE)
gwas$chrPos <- paste0("chr", gwas$Chromosome, ":", gwas$Position)

cnCpGs_bloodmQTLs_annot$chrPos1 <- gsub(":SNP", "", cnCpGs_bloodmQTLs_annot$snp)
cnCpGs_bloodmQTLs_annot$chrPos <- gsub(":INDEL", "", cnCpGs_bloodmQTLs_annot$chrPos1)
cnBlood_gwas <- gwas[which(gwas$chrPos %in% cnCpGs_bloodmQTLs_annot$chrPos),]
cnBlood_gwas$sig <- ifelse(cnBlood_gwas$Pvalue < 5e-8, 1, 0)
write.csv(cnBlood_gwas, "cnBlood_gwas.csv", row.names = FALSE)

adCpGs_bloodmQTLs_anoot$chrPos1 <- gsub(":SNP", "", adCpGs_bloodmQTLs_anoot$snp)
adCpGs_bloodmQTLs_anoot$chrPos <- gsub(":INDEL", "", adCpGs_bloodmQTLs_anoot$chrPos1)
adBlood_gwas <- gwas[which(gwas$chrPos %in% adCpGs_bloodmQTLs_anoot$chrPos),]
adBlood_gwas$sig <- ifelse(adBlood_gwas$Pvalue < 5e-8, 1, 0)
write.csv(adBlood_gwas, "adBlood_gwas.csv", row.names = FALSE)


#Overlap of mQTLs with 24 AD regions
library(GenomicRanges)

#cn
GWASloci <- read.csv("../AD_LDblocks.csv")
cnCpGs_bloodmQTLs_annot$Chr1 <- gsub(":.*", "", cnCpGs_bloodmQTLs_annot$chrPos)
cnCpGs_bloodmQTLs_annot$Chr <- as.numeric(gsub("chr", "", cnCpGs_bloodmQTLs_annot$Chr1))

cnCpGs_bloodmQTLs_annot$Pos <- as.numeric(gsub(".*:", "", cnCpGs_bloodmQTLs_annot$chrPos))

cnCpGs_bloodmQTLs_annot <- cnCpGs_bloodmQTLs_annot[order(cnCpGs_bloodmQTLs_annot$Chr, cnCpGs_bloodmQTLs_annot$Pos),]

GWASloci_ranges <- GRanges(
  seqnames = GWASloci$chr,
  ranges = IRanges(GWASloci$start, GWASloci$end)
)

cnBloodSNPs_ranges <- GRanges(
  seqnames = cnCpGs_bloodmQTLs_annot$Chr,
  ranges = IRanges(cnCpGs_bloodmQTLs_annot$Pos, cnCpGs_bloodmQTLs_annot$Pos)
)

overlap <- findOverlaps(GWASloci_ranges, cnBloodSNPs_ranges)
overlap_df <- as.data.frame(overlap)

GWASlociOverlap <- GWASloci[overlap_df$queryHits,]
cnBloodsnpsOverlap <- cnCpGs_bloodmQTLs_annot[overlap_df$subjectHits,]

GWASloci_cnBloodOverlap <- cbind(GWASlociOverlap, cnBloodsnpsOverlap)

write.csv(GWASloci_cnBloodOverlap, "GWAS_loci_cnBlood_overlap.csv", row.names = FALSE)


#ad
adCpGs_bloodmQTLs_anoot$Chr1 <- gsub(":.*", "", adCpGs_bloodmQTLs_anoot$chrPos)
adCpGs_bloodmQTLs_anoot$Chr <- as.numeric(gsub("chr", "", adCpGs_bloodmQTLs_anoot$Chr1))

adCpGs_bloodmQTLs_anoot$Pos <- as.numeric(gsub(".*:", "", adCpGs_bloodmQTLs_anoot$chrPos))

adBloodSNPs_ranges <- GRanges(
  seqnames = adCpGs_bloodmQTLs_anoot$Chr,
  ranges = IRanges(adCpGs_bloodmQTLs_anoot$Pos, adCpGs_bloodmQTLs_anoot$Pos)
)

overlap <- findOverlaps(GWASloci_ranges, adBloodSNPs_ranges)
overlap_df <- as.data.frame(overlap)

GWASlociOverlap <- GWASloci[overlap_df$queryHits,]
adBloodsnpsOverlap <- adCpGs_bloodmQTLs_anoot[overlap_df$subjectHits,]

GWASloci_adBloodOverlap <- cbind(GWASlociOverlap, adBloodsnpsOverlap)

write.csv(GWASloci_adBloodOverlap, "GWAS_loci_adBlood_overlap.csv", row.names = FALSE)


#Overlap of CpGs with 24 AD regions
#cn

cnCpGs_annot <- annotateCpG(data = cnCpGs)
adCpGs_annot <- annotateCpG(data = adCpGs)

GWASloci_ranges <- GRanges(
  seqnames = paste0("chr", GWASloci$chr),
  ranges = IRanges(GWASloci$start, GWASloci$end)
)

cnBloodCpGs_ranges <- GRanges(
  seqnames = cnCpGs_annot$chr,
  ranges = IRanges(cnCpGs_annot$start, cnCpGs_annot$end)
)

overlap <- findOverlaps(GWASloci_ranges, cnBloodCpGs_ranges)
overlap_df <- as.data.frame(overlap)

# -> no overlap

# GWASlociOverlap <- GWASloci[overlap_df$queryHits,]
# cnBloodCpGoverlap <- cnCpGs[overlap_df$subjectHits,]
# 
# GWASloci_cnBloodCpGoverlap <- cbind(GWASlociOverlap, cnBloodCpGoverlap)
# 
# write.csv(GWASloci_cnBloodCpGoverlap, "GWAS_loci_cnBloodCpGs_overlap.csv", row.names = FALSE)


#ad
adBloodCpGs_ranges <- GRanges(
  seqnames = adCpGs_annot$chr,
  ranges = IRanges(adCpGs_annot$start, adCpGs_annot$end)
)

overlap <- findOverlaps(GWASloci_ranges, adBloodCpGs_ranges)
overlap_df <- as.data.frame(overlap)

# -> no overlap

# GWASlociOverlap <- GWASloci[overlap_df$queryHits,]
# adBloodCpGoverlap <- adCpGs[overlap_df$subjectHits,]
# 
# GWASloci_adBloodCpGoverlap <- cbind(GWASlociOverlap, adBloodCpGoverlap)
# 
# write.csv(GWASloci_adBloodCpGoverlap, "GWAS_loci_adBloodCpGs_overlap.csv", row.names = FALSE)


#Overlap with Deming et al., 2017 - Supp tables 2-4
demsnps <- read.csv("Deming2017_snps_chrPos.csv")
demsnpss<-demsnps[order(demsnps$chr, demsnps$pos),]

demsnps_ranges <- GRanges(
  seqnames = demsnpss$chr,
  ranges = IRanges(demsnpss$pos, demsnpss$pos)
)

#cn
overlap <- findOverlaps(demsnps_ranges, cnBloodSNPs_ranges)
overlap_df <- as.data.frame(overlap)

# -> no overlap

#ad
overlap <- findOverlaps(demsnps_ranges, adBloodSNPs_ranges)
overlap_df <- as.data.frame(overlap)

