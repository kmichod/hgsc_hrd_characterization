library(survminer)
library(broom)
library(tidyr)
library(survival)
library(dplyr)

directory <- "" #set working directory
source(file.path(directory, "code/1.utils.R"))
all_data <- readxl::read_excel(file.path(directory, "data/all_data.xlsx")) %>% filter(!SUID %in% neo_ids)
all_data_VUS <- readxl::read_excel(file.path(directory, "data/all_data_VUS.xlsx")) %>% filter(!SUID %in% neo_ids)

all_data <- all_data %>%
  mutate_at(vars("any_HRD",
                "any_HRD_mutation",
                "germline_variant",
                "germline_carrier_BRCA1", 
                "germline_carrier_BRCA2", 
                "germline_other", 
                "somatic_mutation",
                "somatic_BRCA1", 
                "somatic_BRCA2", 
                "somatic_other",
                "SBS3", 
                "any_HRD_all_genes",
                "any_HRD_mutation_all_genes",
                "germline_variant_all_genes",
                "germline_other_all_genes",
                "somatic_mutation_all_genes",
                "somatic_other_all_genes", 
                "vitalstatus"), as.numeric)

all_data_VUS <- all_data_VUS %>%
  mutate_at(vars("any_HRD",
                 "any_HRD_mutation",
                 "germline_variant",
                 "germline_carrier_BRCA1", 
                 "germline_carrier_BRCA2", 
                 "germline_other", 
                 "somatic_mutation",
                 "somatic_BRCA1", 
                 "somatic_BRCA2", 
                 "somatic_other",
                 "SBS3", 
                 "any_HRD_all_genes",
                 "any_HRD_mutation_all_genes",
                 "germline_variant_all_genes",
                 "germline_other_all_genes",
                 "somatic_mutation_all_genes",
                 "somatic_other_all_genes", 
                 "vitalstatus"), as.numeric)


all_data_B <- all_data %>% filter(Study == "Schildkraut-B" & !is.na(stage_b))
all_data_W <- all_data %>% filter(Study == "Schildkraut-W" & !is.na(stage_b))

all_data_B_VUS <- all_data_VUS %>% filter(Study == "Schildkraut-B" & !is.na(stage_b))
all_data_W_VUS <- all_data_VUS %>% filter(Study == "Schildkraut-W" & !is.na(stage_b))

#Define function that estimates the main effects of hrd features as input "x" on all cause mortality
survival_function <- function(x_1_values, df) {

    results_df <- data.frame()
  
  for (x in x_1_values) {
      survival_model_1 <- coxph(formula = as.formula(paste("Surv(timelastfu, vitalstatus) ~ ", x)), data = df)
      survival_model_2 <- coxph(formula = as.formula(paste("Surv(timelastfu, vitalstatus) ~ ", x, "+ refage + stage_b")), data = df)

      tidy_results_1 <- tidy(survival_model_1, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) %>% 
        select(term, estimate, conf.low, conf.high) %>% 
        dplyr::slice(1) %>%
        mutate("HR_1" = ifelse(is.na(conf.high), "-", 
                               ifelse(conf.high == "Inf", "-", paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")))) %>%
        select(-c("estimate", "conf.low", "conf.high"))
      
      tidy_results_2 <- tidy(survival_model_2, exponentiate = TRUE, conf.int = TRUE, conf.level = 0.95) %>% 
        select(term, estimate, conf.low, conf.high) %>% 
        dplyr::slice(1) %>%
        mutate("HR_2" = ifelse(is.na(conf.high), "-", 
                               ifelse(conf.high == "Inf", "-", paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")))) %>%
        select(-c("estimate", "conf.low", "conf.high"))
      
      merge_1 <- merge(tidy_results_1, tidy_results_2)
      results_df <- rbind(results_df, merge_1)
    }
  return(results_df)
  }

all_features <- c("germline_variant",
                  "germline_carrier_BRCA1", 
                  "germline_carrier_BRCA2", 
                  "germline_other", 
                  "somatic_mutation",
                  "somatic_BRCA1", 
                  "somatic_BRCA2", 
                  "somatic_other",
                  "SBS3",
                  "any_HRD",
                  "germline_variant_all_genes",
                  "somatic_mutation_all_genes")

#apply function 
main_effect_B <- survival_function(all_features, all_data_B)
main_effect_W <- survival_function(all_features, all_data_W)

main_effect_B_VUS <- survival_function(all_features, all_data_B_VUS)
main_effect_W_VUS <- survival_function(all_features, all_data_W_VUS)

#Define function that estimates the counts for cases and deaths for survival analysis
main_effect_counts <- function(feature_list, df) {
  
  results_df <- data.frame(term = character(),
                           Value = numeric(), 
                           No_cases = numeric(), 
                           No_deaths = numeric(),
                           stringsAsFactors = FALSE)
  
  for (x in feature_list) {
    
    surv_model <- survfit(formula = as.formula(paste("Surv(timelastfu, vitalstatus) ~ ", x)), data = df)
    
    results_df <- rbind(results_df, data.frame(term = c(x, x),
                                               Value = c(0, 1), 
                                               No_cases = c(sum(surv_model[1]$n), sum(surv_model[2]$n)), 
                                               No_deaths = c(sum(surv_model[1]$n.event), sum(surv_model[2]$n.event))))
  }
  
  return(results_df)
}

#Apply function 
main_effect_counts_B <- main_effect_counts(all_features, all_data_B)
main_effect_counts_W <- main_effect_counts(c("any_HRD",
                                             "germline_variant",
                                             "germline_carrier_BRCA1",
                                             "germline_carrier_BRCA2",
                                             "germline_other", 
                                             "somatic_mutation",
                                             "somatic_BRCA1",
                                             "somatic_BRCA2",
                                             "SBS3", 
                                             "germline_variant_all_genes",
                                             "somatic_mutation_all_genes"), all_data_W) 

main_effect_counts_ALL <- main_effect_counts(all_features, all_data)
survival_main_effect_B <- merge(main_effect_counts_B, main_effect_B, by = "term", all = TRUE) %>% mutate(HR_1 = ifelse(Value == 0, "Ref", HR_1),
                                                                                                         HR_2 = ifelse(Value == 0, "Ref", HR_2))
                                                                                                                              
colnames(survival_main_effect_B) <- paste0("B_", colnames(survival_main_effect_B))
survival_main_effect_B <- survival_main_effect_B %>% dplyr::rename("term" = "B_term", 
                                                            "Value" = "B_Value")

survival_main_effect_W <- merge(main_effect_counts_W, main_effect_W, by = "term", all = TRUE) %>% mutate(HR_1 = ifelse(Value == 0, "Ref", HR_1),
                                                                                                         HR_2 = ifelse(Value == 0, "Ref", HR_2))
colnames(survival_main_effect_W) <- paste0("W_", colnames(survival_main_effect_W))
survival_main_effect_W <- survival_main_effect_W %>% dplyr::rename("term" = "W_term", 
                                                            "Value" = "W_Value")

merge_main_effect <- merge(survival_main_effect_B, survival_main_effect_W, all = TRUE)
merge_main_effect <- merge_main_effect %>% arrange(factor(term, levels = c("germline_carrier_BRCA1", 
                                                                           "germline_carrier_BRCA2", 
                                                                           "germline_other", 
                                                                           "germline_variant",
                                                                           "somatic_BRCA1", 
                                                                           "somatic_BRCA2", 
                                                                           "somatic_other",
                                                                           "somatic_mutation",
                                                                           "SBS3",
                                                                           "any_HRD",
                                                                           "total_perMB",
                                                                           "any_HRD_all_genes",
                                                                           "any_HRD_mutation_all_genes",
                                                                           "germline_variant_all_genes",
                                                                           "germline_other_all_genes",
                                                                           "somatic_mutation_all_genes",
                                                                           "somatic_other_all_genes"))) %>%
  filter(!is.na(Value))

#Repeat among VUS's 
main_effect_counts_ALL_VUS <- main_effect_counts(all_features, all_data_VUS)
surv_main_effect_all_VUS <- merge(main_effect_counts_ALL_VUS, main_effect_ALL_VUS, by = "term", all = TRUE) %>% mutate(HR_1 = ifelse(Value == 0, "Ref", HR_1),
                                                                                                         HR_2 = ifelse(Value == 0, "Ref", HR_2))

main_effect_counts_B_VUS <- main_effect_counts(all_features, all_data_B_VUS)
surv_main_effect_B_VUS <- merge(main_effect_counts_B_VUS, main_effect_B_VUS, by = "term", all = TRUE) %>% mutate(HR_1 = ifelse(Value == 0, "Ref", HR_1),
                                                                                                                       HR_2 = ifelse(Value == 0, "Ref", HR_2))
colnames(surv_main_effect_B_VUS) <- paste0("B_", colnames(surv_main_effect_B_VUS))
surv_main_effect_B_VUS <- surv_main_effect_B_VUS %>% dplyr::rename("term" = "B_term", 
                                                                   "Value" = "B_Value")
somatic_other_W_VUS <- survfit(formula = as.formula(paste("Surv(timelastfu, vitalstatus) ~ ", "somatic_other")), data = all_data_W_VUS)
main_effect_counts_W_VUS <- main_effect_counts(c("any_HRD",
                                               "germline_variant",
                                               "germline_carrier_BRCA1",
                                               "germline_carrier_BRCA2",
                                               "germline_other", 
                                               "somatic_mutation",
                                               "somatic_BRCA1",
                                               "somatic_BRCA2",
                                               "SBS3", 
                                               "germline_variant_all_genes",
                                               "somatic_mutation_all_genes"), all_data_W_VUS) %>%
  rbind(data.frame(term = c("somatic_other", "somatic_other"),
                   Value = c(0, 1), 
                   No_cases = c(sum(somatic_other_W_VUS[1]$n), 0), 
                   No_deaths = c(sum(somatic_other_W_VUS[1]$n.event), 0)))

surv_main_effect_W_VUS <- merge(main_effect_counts_W_VUS, main_effect_W_VUS, by = "term", all = TRUE) %>% mutate(HR_1 = ifelse(Value == 0, "Ref", HR_1),
                                                                                                                 HR_2 = ifelse(Value == 0, "Ref", HR_2))
colnames(surv_main_effect_W_VUS) <- paste0("W_", colnames(surv_main_effect_W_VUS))
surv_main_effect_W_VUS <- surv_main_effect_W_VUS %>% dplyr::rename("term" = "W_term", 
                                                                   "Value" = "W_Value")

merge_main_effect_temp_VUS <- merge(surv_main_effect_all_VUS, surv_main_effect_B_VUS)
merge_main_effect_VUS <- merge(merge_main_effect_temp_VUS, surv_main_effect_W_VUS, all = TRUE)
merge_main_effect_VUS <- merge_main_effect_VUS %>% arrange(factor(term, levels = c("any_HRD",
                                                                                   "germline_variant",
                                                                                   "germline_carrier_BRCA1",
                                                                                   "germline_carrier_BRCA2",
                                                                                   "germline_other", 
                                                                                   "somatic_mutation",
                                                                                   "somatic_BRCA1",
                                                                                   "somatic_BRCA2",
                                                                                   "SBS3", 
                                                                                   "germline_variant_all_genes",
                                                                                   "somatic_mutation_all_genes"))) %>%
  filter(!is.na(Value))










