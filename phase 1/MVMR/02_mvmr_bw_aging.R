########################################################
# MVMR Analysis: Birthweight to healthy aging phenotypes (with adjustment for obesity indicators)
########################################################

library(data.table)
library(MRPRESSO)
library(readr)
library(vroom)
library(tidyr)
library(tibble)
library(dplyr)
library(MVMR)
library(MendelianRandomization)
library(TwoSampleMR)
library(writexl)


# Read outcome GWAS 
GIP1_1 <-fread("GWAS_outcome/aging_gip.tsv.gz",header=T) 
GIP1_1 <- as.data.frame(GIP1_1)

mvAge <-fread("GWAS_outcome/mvAge.txt",header=T)
mvAge <- as.data.frame(mvAge)

healthspan <-fread("GWAS_outcome/healthspan.txt",header=T) 
healthspan <- as.data.frame(healthspan)

Plifespan <-fread("GWAS_outcome/parental_lifespan.tsv.gz",header=T)
Plifespan <- as.data.frame(Plifespan)

resilience <- fread("GWAS_outcome/resilience.txt",header=T,fill=T)
resilience <- as.data.frame(resilience)


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


# Start: Main analysis ----
# Results in Table S9
setwd("MVMR_IV_bw_adjobs_main")
filename <- dir()
all_results <- list()

for (i in 1:length(filename)) {
  exp <- read.csv(filename[i])
  rawdatafinal <- subset(exp, select = c(effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure, SNP, exposure, id.exposure))
  
  outcomes <- list(
    list(data = GIP1_1, name = "GIP1", snp_col = "rsid", beta_col = "beta1", se_col = "se", eaf_col = "freq1", effect_allele_col = "a1", other_allele_col = "a0", pval_col = "p", chr_col = "chr", pos_col = "pos", samplesize_col = "n"),
    list(data = mvAge, name = "mvAge", snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "MAF", effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "Pvalue", chr_col = "CHR", pos_col = "BP", samplesize_col = "n"),
    list(data = healthspan, name = "healthspan", snp_col = "SNPID", beta_col = "beta", se_col = "se", eaf_col = "EAF", effect_allele_col = "EA", other_allele_col = "RA", pval_col = "pval", chr_col = "chr", pos_col = "pos", samplesize_col = "n"),
    list(data = resilience, name = "resilience", snp_col = "SNP", beta_col = "BETA", se_col = "SE", eaf_col = "A1FREQ", effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", pval_col = "P_BOLT_LMM", chr_col = "CHR", pos_col = "BP", samplesize_col = "n"),
    list(data = Plifespan, name = "Plifespan", snp_col = "rsid", beta_col = "beta1", se_col = "se", eaf_col = "freq1", effect_allele_col = "a1", other_allele_col = "a0", pval_col = "p", chr_col = "chr", pos_col = "pos", samplesize_col = "n")
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

mvmr_results <- grep("_mvmr$", names(all_results), value = TRUE)
combined_mvmr <- do.call(rbind, all_results[mvmr_results])

mvegger_results <- grep("_mvegger$", names(all_results), value = TRUE)
combined_mvegger <- do.call(rbind, all_results[mvegger_results])

# End: Main analysis ----



# Start: Replication analysis ----
# Results in Table S10
setwd("MVMR_IV_bw_adjobs_replication")
filename <- dir()
all_results <- list()

for (i in 1:length(filename)) {
  exp <- read.csv(filename[i])
  rawdatafinal <- subset(exp, select = c(effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure, SNP, exposure, id.exposure))
  
  outcomes <- list(
    list(data = GIP1_1, name = "GIP1", snp_col = "rsid", beta_col = "beta1", se_col = "se", eaf_col = "freq1", effect_allele_col = "a1", other_allele_col = "a0", pval_col = "p", chr_col = "chr", pos_col = "pos", samplesize_col = "n"),
    list(data = mvAge, name = "mvAge", snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "MAF", effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "Pvalue", chr_col = "CHR", pos_col = "BP", samplesize_col = "n"),
    list(data = healthspan, name = "healthspan", snp_col = "SNPID", beta_col = "beta", se_col = "se", eaf_col = "EAF", effect_allele_col = "EA", other_allele_col = "RA", pval_col = "pval", chr_col = "chr", pos_col = "pos", samplesize_col = "n"),
    list(data = resilience, name = "resilience", snp_col = "SNP", beta_col = "BETA", se_col = "SE", eaf_col = "A1FREQ", effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", pval_col = "P_BOLT_LMM", chr_col = "CHR", pos_col = "BP", samplesize_col = "n"),
    list(data = Plifespan, name = "Plifespan", snp_col = "rsid", beta_col = "beta1", se_col = "se", eaf_col = "freq1", effect_allele_col = "a1", other_allele_col = "a0", pval_col = "p", chr_col = "chr", pos_col = "pos", samplesize_col = "n") 
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

mvmr_results <- grep("_mvmr$", names(all_results), value = TRUE)
combined_mvmr <- do.call(rbind, all_results[mvmr_results])

mvegger_results <- grep("_mvegger$", names(all_results), value = TRUE)
combined_mvegger <- do.call(rbind, all_results[mvegger_results])

# End: Main analysis ----
