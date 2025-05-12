# UVMR Analysis: Mediators to healthy aging phenotypes

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

ANALYSIS_CONFIG <- list(
  iv_input = list(
    online_ids = c('ieu-a-89','ebi-a-GCST003156','ebi-a-GCST005195','ebi-a-GCST006414','ebi-a-GCST006906','ebi-a-GCST009541','ebi-a-GCST011365','ebi-a-GCST011494','ebi-a-GCST90000025','ebi-a-GCST90000047','ebi-a-GCST90000050','ebi-a-GCST90002232','ebi-a-GCST90002238','ebi-a-GCST90002244','ebi-a-GCST90012104','ebi-a-GCST90012111','ebi-a-GCST90012114','ebi-a-GCST90013405','ebi-a-GCST90013534','ebi-a-GCST90018890','ebi-a-GCST90027158','ebi-a-GCST90038604','finn-b-AB1_GASTROENTERITIS_NOS','finn-b-ASTHMA_ACUTE_RESPIRATORY_INFECTIONS','finn-b-ASTHMA_CHILD','finn-b-C3_COLORECTAL','finn-b-F5_DEMENTIA','finn-b-I9_ANGINA','finn-b-I9_CORATHER','finn-b-I9_ISCHHEART','finn-b-I9_UAP','finn-b-I9_VTE','finn-b-J10_ASTHMA','finn-b-J10_COPD','finn-b-M13_LOWBACKPAIN','finn-b-M13_OSTEOPOROSIS','finn-b-M13_POLYARTHROPATHIES','ieu-a-1126','ieu-a-1183','ieu-a-1239','ieu-a-31','ieu-a-7','ieu-a-966','ieu-b-102','ieu-b-142','ieu-b-29','ieu-b-30','ieu-b-31','ieu-b-32','ieu-b-33','ieu-b-34','ieu-b-35','ieu-b-38','ieu-b-39','ieu-b-4877','ieu-b-5102','ieu-b-7','ieu-b-73','ieu-b-85','met-c-842','met-c-843','ukb-b-1209','ukb-b-1489','ukb-b-1996','ukb-b-3881','ukb-b-4711','ukb-b-5237','ukb-b-6066','ukb-b-7408','ukb-b-7478','ukb-b-8476'),
    local_dir = "data/mediator/",
    output_path = "IV/mediator/"
  ),
  
  clump_params = list(r2 = 0.001),
  fstat_threshold = 10,
  
  outcomes = list(
    GIP = list(
      data = GIP,
      params = list(
        snp_col = "rsid",
        beta_col = "beta1",
        se_col = "se",
        eaf_col = "freq1",
        effect_allele_col = "a1",
        other_allele_col = "a0"
      )
    ),
    mvAge = list(
      data = mvAge,
      params = list(
        snp_col = "SNP",
        beta_col = "beta",
        se_col = "se",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele"
      )
    )
  )
)

# -------------------------------
# 2. Core Functions
# -------------------------------
prepare_iv_data <- function(config) {
  tryCatch({
    # online
    online_iv <- extract_instruments(
      outcomes = config$iv_input$online_ids,
      r2 = config$clump_params$r2
    )
    
    # local
    local_files <- list.files(
      config$iv_input$local_dir, 
      pattern = "\\.csv$", 
      full.names = TRUE
    )
    
    local_iv <- map_dfr(local_files, ~{
      df <- fread(.x) 
      df %>% mutate(exposure = tools::file_path_sans_ext(basename(.x)))
    })
    
    bind_rows(online_iv, local_iv) %>%
      write_csv(config$iv_input$output_path)
    
    message("IV data prepared successfully")
    return(read_exposure_data(config$iv_input$output_path))
    
  }, error = function(e) {
    message("IV preparation failed: ", e$message)
    return(NULL)
  })
}

analyze_outcome <- function(exposure_data, outcome_config) {
  tryCatch({
    formatted <- do.call(format_data, c(
      list(
        dat = outcome_config$data,
        type = "outcome",
        snps = exposure_data$SNP,
        header = TRUE,
        phenotype_col = "Phenotype"
      ),
      outcome_config$params
    ))
    
    harmonized <- harmonise_data(exposure_data, formatted)

    mr_res <- mr(harmonized) %>% generate_odds_ratios()
    pleio <- mr_pleiotropy_test(harmonized)
    heter <- mr_heterogeneity(harmonized)

    fstat <- harmonized %>%
      mutate(fval = (beta.exposure/se.exposure)^2) %>%
      group_by(exposure,outcome) %>%
      summarise(Fstat = mean(fval, na.rm = TRUE))
    
    list(
      harmonized = harmonized,
      mr_result = mr_res,
      pleiotropy = pleio,
      heterogeneity = heter,
      fstat = fstat
    )
    
  }, error = function(e) {
    message("Analysis failed for ", outcome_config$name, ": ", e$message)
    return(NULL)
  })
}

# MR-PRESSO analysis
run_mr_presso_analysis <- function(harmonized_data) {
  tryCatch({
    filtered_data <- harmonized_data %>%
      group_by(exposure) %>%
      filter(n() > 3) %>%
      ungroup()
    
    run_mr_presso(filtered_data, 
                  NbDistribution = 1000, 
                  SignifThreshold = 0.05)
    
  }, error = function(e) {
    message("MR-PRESSO failed: ", e$message)
    return(NULL)
  })
}

# -------------------------------
# 3. Main Analysis Pipeline
# -------------------------------
main_analysis <- function(config) {
  exposure_data <- prepare_iv_data(config)
  if(is.null(exposure_data)) return()
  
  results <- map(config$outcomes, ~analyze_outcome(exposure_data, .x))
  
  final_results <- map_dfr(results, ~{
    if(!is.null(.x)) {
      compile_results(.x$mr_result, 
                      .x$pleiotropy, 
                      .x$heterogeneity, 
                      .x$fstat)
    }
  })
  
  presso_results <- run_mr_presso_analysis(
    bind_rows(map(results, ~.x$harmonized))
  )
  
  write_csv(final_results, "final_mr_results.csv")
  saveRDS(presso_results, "presso_results.rds")
  
  message("Analysis completed successfully")
}

# -------------------------------
# 4. Execute Analysis
# -------------------------------
main_analysis(ANALYSIS_CONFIG)