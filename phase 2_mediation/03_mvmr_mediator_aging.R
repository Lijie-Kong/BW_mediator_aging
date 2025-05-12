# MVMR Analysis: Mediators to Healthy Aging Traits (with adjustment for birthweight)

# -------------------------------
# 1. Environment Setup
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


# -------------------------------
# 2: Extraction of MVMR IV
# -------------------------------

# Load BW exposure data
BW_IV <- read_exposure_data(
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

exposure_types <- unique(BW_IV$exposure)
gwas_path_map <- list(
  "Overall birthweight" = "data/gwas/BW.txt"
)

setwd("data/mediator/")
filename <- dir()
outname <- gsub("\\.txt$", "", filename)

for (exposure_type in exposure_types) {
  BW_IV2 <- BW_IV %>% filter(exposure == exposure_type)
  
  gwas_path <- ifelse(
    exposure_type %in% names(gwas_path_map),
    gwas_path_map[[exposure_type]],
    stop(paste("Undefined exposure type path:", exposure_type))
  )
  
  for(i in 1:length(filename))
  {
    exp <-fread(filename[i],header=T,fill=T) 
    exp <- as.data.frame(exp)
    exp$pval <- as.numeric(exp$pval)
    exp_IV<-exp[exp$pval<5*10^-8,]
    exp_IV <- format_data(
      exp_IV,
      type = "exposure",
      snps = NULL,
      header = TRUE,
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
    exp_IV2<-clump_data(exp_IV,clump_r2 = 0.001)
    
    exp2_exp <- exp_IV2
    exp1_exp2 <- subset(BW_IV2,select = c(effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure,se.exposure,pval.exposure,SNP,exposure,id.exposure))
    exp2_exp2 <- subset(exp2_exp,select = c(effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure,se.exposure,pval.exposure,SNP,exposure,id.exposure))
    exposure_dat <- rbind(exp1_exp2,exp2_exp2)
    
    temp <- exposure_dat
    temp$id.exposure <- 1
    temp <- temp[order(temp$pval.exposure, decreasing = FALSE), ]
    temp <- subset(temp, !duplicated(SNP))
    temp2 <- clump_data(temp,clump_r2 = 0.001)
    exposure_dat <- subset(exposure_dat, SNP %in% temp2$SNP)
    ####MVMR IVs
    
    #BW GWAS
    raw1 <- read_outcome_data(
      filename  = gwas_path,  
      sep = " ",
      snps = exposure_dat$SNP,
      snp_col = "RSID",
      beta_col = "beta",
      se_col = "se",
      eaf_col = "eaf",
      effect_allele_col = "ea",
      other_allele_col = "nea",
      pval_col = "p",
      # samplesize_col = "samplesize",
      chr_col = "chr",
      pos_col = "pos"
    )
    data1 <- convert_outcome_to_exposure(raw1)
    
    #exp2 GWAS
    raw2 <- format_data(
      exp,
      type = "outcome",
      snps = data1$SNP,
      header = TRUE,
      snp_col = "SNP",
      beta_col = "BETA",
      se_col = "SE",
      eaf_col = "EAF",
      effect_allele_col = "EA",
      other_allele_col = "OA",
      pval_col = "pval",
      units_col = "units",
      ncase_col = "af_alt_cases",
      ncontrol_col = "af_alt_controls",
      samplesize_col = "samplesize"
    )
    rawdata<- harmonise_data(exposure_dat  = data1, 
                             outcome_dat = raw2)
    
    rawdata1 <- subset(rawdata,select = c(SNP, exposure, id.exposure, effect_allele.exposure , 
                                          other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, 
                                          pval.exposure))
    rawdata1$exposure <- exposure_type
    
    rawdata2 <- subset(rawdata,select = c(SNP, outcome, id.outcome, effect_allele.outcome, 
                                          other_allele.outcome, eaf.outcome, beta.outcome, se.outcome, 
                                          pval.outcome))
    names(rawdata2) <- gsub("outcome", "exposure", names(rawdata2))
    rawdata2$exposure <- outname[i]
    rawdatafinal <- rbind(rawdata1, rawdata2)
    write.csv(rawdatafinal,output_filename)
  }
}

# -------------------------------
# 3: MVMR analysis
# -------------------------------

process_outcome <- function(rawdata, outcome_data, outcome_name, snp_col, beta_col, se_col, eaf_col, effect_allele_col, other_allele_col, pval_col, chr_col, pos_col, samplesize_col) {
  outcome_mvout <- format_data(
    outcome_data,
    type = "outcome",
    snps = rawdata$SNP,
    header = TRUE,
    phenotype_col = "Phenotype",
    snp_col = snp_col,
    beta_col = beta_col,
    se_col = se_col,
    eaf_col = eaf_col,
    effect_allele_col = effect_allele_col,
    other_allele_col = other_allele_col,
    pval_col = pval_col,
    chr_col = chr_col,
    pos_col = pos_col,
    samplesize_col = samplesize_col
  )
  outcome_mvout$outcome <- outcome_name
  mvdat <- mv_harmonise_data(exposure_dat=rawdata, outcome_dat = outcome_mvout)
  res <- mv_multiple(mvdat)
  MRMVObject <- mr_mvinput(bx = cbind(mvdat$exposure_beta), bxse = cbind(mvdat$exposure_se), 
                           by = mvdat$outcome_beta, byse=mvdat$outcome_se)
  MVMR1 <- mr_mvivw(MRMVObject)
  MVMR2 <- mr_mvegger(MRMVObject)
  MVMR3 <- mr_mvmedian(MRMVObject)
  MVMR4 <- mr_mvlasso(MRMVObject)
  MVMRres <- bind_cols(MVMR1@Estimate, MVMR1@StdError, MVMR1@Pvalue, MVMR2@Estimate, MVMR2@StdError.Est, MVMR2@Pvalue.Est, MVMR3@Estimate, MVMR3@StdError, MVMR3@Pvalue, MVMR4@Estimate, MVMR4@StdError, MVMR4@Pvalue)
  colnames(MVMRres) <- c("IVW_beta", "IVW_SE", "IVW_P", "Egger_beta", "Egger_SE", "Egger_P", "median_beta", "median_SE", "median_P", "lasso_beta", "lasso_SE", "lasso_P")
  MVMRres <- cbind(MVMRres, phenotype = mvdat$expname)
  return(list(result = res$result, sensitivity = MVMRres))
}

setwd("data/MVIV/")
filename <- dir()
all_results <- list()

for (i in 1:length(filename)) {
  exp <- read.csv(filename[i])
  rawdatafinal <- subset(exp, select = c(effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure, SNP, exposure, id.exposure))
  
  outcomes <- list(
    list(data = GIP, name = "GIP", snp_col = "rsid", beta_col = "beta1", se_col = "se", eaf_col = "freq1", effect_allele_col = "a1", other_allele_col = "a0", pval_col = "p", chr_col = "chr", pos_col = "pos", samplesize_col = "n"),
    list(data = mvAge, name = "mvAge", snp_col = "SNP", beta_col = "beta", se_col = "se", effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "Pvalue", chr_col = "CHR", pos_col = "BP", samplesize_col = "n")
  )
  
  for (outcome in outcomes) {
    res <- process_outcome(rawdatafinal, outcome$data, outcome$name, outcome$snp_col, outcome$beta_col, outcome$se_col, outcome$eaf_col, outcome$effect_allele_col, outcome$other_allele_col, outcome$pval_col, outcome$chr_col, outcome$pos_col, outcome$samplesize_col)
    
    result_key <- paste0(outcome$name, "_result_mvmr")
    sensitivity_key <- paste0(outcome$name, "_result_mvegger")
    
    if (i == 1) {
      all_results[[result_key]] <- res$result
      all_results[[sensitivity_key]] <- res$sensitivity
    } else {
      all_results[[result_key]] <- rbind(all_results[[result_key]], res$result)
      all_results[[sensitivity_key]] <- rbind(all_results[[sensitivity_key]], res$sensitivity)
    }
  }
}


