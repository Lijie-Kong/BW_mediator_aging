########################################################
# UVMR Analysis: Birthweight to healthy aging phenotypes
########################################################

# Load required libraries
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


# Start: Main analysis ----
# Results in Table S3 and Table S4
harmonized_data <- fread("Hdata_BW_aging_main.csv")

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

# End: Main analysis ----


# Start: Replication analysis ----
# Results in Table S7 and Table S8
harmonized_data <- fread("Hdata_BW_aging_replication.csv")

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

# End: Main analysis ----
