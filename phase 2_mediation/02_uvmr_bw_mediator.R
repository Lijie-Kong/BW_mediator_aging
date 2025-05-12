# UVMR Analysis: Birthweight to Mediators

# -------------------------------
# 1. Configuration Settings
# -------------------------------
library(MRInstruments)
library(data.table)
library(MRPRESSO)
library(readr)
library(vroom)
library(tidyr)
library(tibble)
library(dplyr)
library(MVMR)
library(ieugwasr)
library(MendelianRandomization)

# --------------------------
# Section 2: Data Loading
# --------------------------
# Load exposure data
exposure_dat <- read_exposure_data(
  filename = "data/exposure/BW_IV.csv",
  sep = ",",
  phenotype_col = "exposure",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure",
  chr_col = "chr.exposure",
  pos_col = "pos.exposure"
)
exposure_dat <- subset(exposure_dat, (exposure %in% c("Overall birthweight")))

# Load mediator data
setwd("data/mediator/")
filename <- dir()
outname <- gsub("\\.txt$", "", filename)
for( i in 1:length(filename))
{
  out_gwas <-fread(filename[i],header=T,fill=T)
  out_gwas$pval <- as.numeric(out_gwas$pval)
  out_gwas <- as.data.frame(out_gwas)
  out_extra <- format_data(
    out_gwas,
    type = "outcome",
    snps = exposure_dat$SNP,
    header = TRUE,
    phenotype_col = "phenotype",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "EAF",
    effect_allele_col = "EA",
    other_allele_col = "OA",
    pval_col = "pval",
    chr_col = "Chr",
    pos_col = "BP",
    samplesize_col = "samplesize"
  )
  out_extra$outcome <- outname[i]
  
  if ("chr.outcome" %in% colnames(out_extra)) {
    out_extra$chr.outcome <- as.character(out_extra$chr.outcome)
  }
  if ("pos.outcome" %in% colnames(out_extra)) {
    out_extra$pos.outcome <- as.character(out_extra$pos.outcome)
  }
  
  H_data <- harmonise_data(exposure_dat  = exposure_dat, 
                           outcome_dat =  out_extra)
  mr_results<-mr(H_data)
  res <- generate_odds_ratios(mr_results)
  # sensitivity analysis
  pleio <- mr_pleiotropy_test(H_data)
  heter <- mr_heterogeneity(H_data)
  # F-statistics 
  H2 <- H_data
  H2$fval <- (H2$beta.exposure/H2$se.exposure)^2
  f <- aggregate(H2$fval,by = list(H2$outcome),mean)
  
  if (i == 1) {
    hamdata <- H_data
    result <- res
    ple <- pleio
    het <- heter
    fstat <-f 
  }else{
    hamdata_new<- H_data
    result_new <- res
    ple_new <- pleio
    het_new <- heter
    fstat_new <-f 
    hamdata <- bind_rows(hamdata,hamdata_new)
    result <- bind_rows(result, result_new)
    ple <- bind_rows(ple, ple_new)
    het <- bind_rows(het, het_new)
    fstat <- bind_rows(fstat,fstat_new)
  }
}



