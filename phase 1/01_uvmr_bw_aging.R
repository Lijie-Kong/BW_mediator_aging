# UVMR Analysis: Birthweight to healthy aging phenotypes

# --------------------------
# Section 1: Environment Setup
# --------------------------

# Load required libraries
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


# --------------------------
# Section 3: Outcome Data Preparation
# --------------------------
GIP <-fread("data/outcomes/GIP.txt",header=T) 
mvAge <-fread("data/outcomes/mvAge.txt",header=T) 
healthspan <-fread("data/outcomes/healthspan.txt",header=T) 
Plifespan <-fread("data/outcomes/Plifespan.txt",header=T) 
longevity <-fread("data/outcomes/longevity.txt",header=T) 
srh <-fread("data/outcomes/srh.txt",header=T) 
resilience <-fread("data/outcomes/resilience.txt",header=T) 
PhenoAge <-fread("data/outcomes/PhenoAge.txt",header=T) 

# format outcome GWAS ----
GIP_out <- format_data(
  GIP,
  type = "outcome",
  snps = exposure_dat$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "rsid",
  beta_col = "beta1",
  se_col = "se",
  eaf_col = "freq1",
  effect_allele_col = "a1",
  other_allele_col = "a0",
  pval_col = "p",
  chr_col = "chr",
  pos_col = "pos",
  samplesize_col = "n"
)
GIP_out$outcome <- "aging-GIP"

mvAge_out <- format_data(
  mvAge,
  type = "outcome",
  snps = exposure_dat$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  # eaf_col = "MAF",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "Pvalue",
  chr_col = "CHR",
  pos_col = "BP"
)
mvAge_out$outcome <- "mvAge"

healthspan_out <- format_data(
  healthspan,
  type = "outcome",
  snps = exposure_dat$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNPID",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "EAF",
  effect_allele_col = "EA",
  other_allele_col = "RA",
  pval_col = "pval",
  chr_col = "chr",
  pos_col = "pos",
  samplesize_col = "n"
)
healthspan_out$outcome <- "healthspan"

resilience_out <- format_data(
  resilience,
  type = "outcome",
  snps = exposure_dat$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "A1FREQ",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  pval_col = "P_BOLT_LMM",
  chr_col = "CHR",
  pos_col = "BP",
  samplesize_col = "n"
)
resilience_out$outcome <- "resilience"

longevity_out <- format_data(
  longevity,
  type = "outcome",
  snps = exposure_dat$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "rsID",
  beta_col = "Beta",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "EA",
  other_allele_col = "NEA",
  pval_col = "P-value",
  chr_col = "Chr",
  pos_col = "Position",
  samplesize_col = "Effective_N"
)
longevity_out$outcome <- "longevity"


srh_out <- format_data(
  srh,
  type = "outcome",
  snps = exposure_dat$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "SNP",
  beta_col = "Beta",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  chr_col = "CHR",
  pos_col = "POS",
  samplesize_col = "N"
)
srh_out$outcome <- "Self-rated health"

Plifespan_out <- format_data(
  Plifespan,
  type = "outcome",
  snps = exposure_dat$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "rsid",
  beta_col = "beta1",
  se_col = "se",
  eaf_col = "freq1",
  effect_allele_col = "a1",
  other_allele_col = "a0",
  pval_col = "p",
  chr_col = "chr",
  pos_col = "pos",
  samplesize_col = "n"
)
Plifespan_out$outcome <- "Parental lifespan"

PhenoAge_out <- format_data(
  PhenoAge,
  type = "outcome",
  snps = exposure_dat$SNP,
  header = TRUE,
  phenotype_col = "Phenotype",
  snp_col = "rsID",
  beta_col = "Effect",
  se_col = "SE",
  eaf_col = "Freq1",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  chr_col = "chr",
  pos_col = "bp",
  samplesize_col = "N"
)
PhenoAge_out$outcome <- "PhenoAge deceleration"


# format data
common_columns <- Reduce(intersect, lapply(list(GIP_out, mvAge_out,healthspan_out,resilience_out,longevity_out,srh_out,Plifespan_out,PhenoAge_out), colnames))
outdata <- rbind(GIP_out[, common_columns],
                 mvAge_out[, common_columns],
                 healthspan_out[, common_columns],
                 resilience_out[, common_columns],
                 longevity_out[, common_columns],
                 srh_out[, common_columns],
                 Plifespan_out[, common_columns],
                 PhenoAge_out[, common_columns]
)

# --------------------------
# Section 4: Data Harmonization
# --------------------------

harmonized_data <- harmonise_data(
  exposure_dat = exposure_dat,
  outcome_dat = outdata
)

# --------------------------
# Section 5: Main MR Analysis
# --------------------------

# Perform MR with multiple methods
mr_results<-mr(harmonized_data)
res <- generate_odds_ratios(mr_results)

# sensitivity analysis
pleiotropy <- mr_pleiotropy_test(harmonized_data)
heterogeneity <- mr_heterogeneity(harmonized_data)

# MR-PRESSO
presso_results <- harmonized_data %>%
  group_by(exposure) %>%
  filter(n() > 3) %>%
  run_mr_presso(NbDistribution = 1000)
