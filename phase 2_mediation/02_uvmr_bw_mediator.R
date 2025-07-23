
##############################################
# UVMR Analysis: Birthweight to Mediators
# Figure 3, Figure 4_left panel, Tables S11-12
##############################################

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
harmonized_data <- fread("Hdata_bw_pmed.csv")

# Perform MR with multiple methods
mr_results<-mr(harmonized_data)
res <- generate_odds_ratios(mr_results)

# sensitivity analysis
pleiotropy <- mr_pleiotropy_test(harmonized_data)
heterogeneity <- mr_heterogeneity(harmonized_data)

# End: Main analysis ----




