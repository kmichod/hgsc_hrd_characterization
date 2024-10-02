library(broom)
library(tidyr)
library(dplyr)

#change directory to hrd folder location 
directory <- "" #set working directory
source(file.path(directory, "code/1.utils.R"))

all_data <- readxl::read_excel(file.path(directory, "/data/all_data.xlsx")) %>% filter(!SUID %in% neo_ids)
all_data_VUS <- readxl::read_excel(file.path(directory, "data/all_data_VUS.xlsx")) %>% filter(!SUID %in% neo_ids)

#recode character variables to numeric for analysis 
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

#stratify data by self-reported race
all_data_B <- all_data %>% filter(Study == "Schildkraut-B")
all_data_W <- all_data %>% filter(Study == "Schildkraut-W")

all_data_VUS_B <- all_data_VUS %>% filter(Study == "Schildkraut-B")
all_data_VUS_W <- all_data_VUS %>% filter(Study == "Schildkraut-W")

#define function that will calculate the median change in age and tumor mutation burden for input x
mediator_function <- function(x, df) {
  
    df_filtered <- df %>% filter(!is.na({{x}}))
    age_model <- lm(formula = paste("refage ~", x), data = df)
    tmb_model <- lm(formula = paste("total_perMB ~", x), data = df)
  
    tidy_results_1 <- tidy(age_model, exponentiate = FALSE, conf.int = TRUE, conf.level = 0.95) %>% 
    select(term, estimate, conf.low, conf.high) %>% 
    dplyr::slice(2) %>%
    mutate("RR_age" = ifelse(is.na(conf.high), "-", 
                             ifelse(conf.high == "Inf", "-", paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")))) %>%
    select(-c("estimate", "conf.low", "conf.high"))
    
    tidy_results_2 <- tidy(tmb_model, exponentiate = FALSE, conf.int = TRUE, conf.level = 0.95) %>% 
      select(term, estimate, conf.low, conf.high) %>% 
      dplyr::slice(2) %>%
      mutate("RR_tmb" = ifelse(is.na(conf.high), "-", 
                               ifelse(conf.high == "Inf", "-", paste0(round(estimate, 2), " (", round(conf.low, 2), ", ", round(conf.high, 2), ")")))) %>%
      select(-c("estimate", "conf.low", "conf.high"))
    
    merge_1 <- merge(tidy_results_1, tidy_results_2) %>% mutate(value = 1)
    temp_df <- data.frame(term = {{x}}, 
                          value = 0, 
                          RR_age = "Ref", 
                          RR_tmb = "Ref")
    
    results <- merge_1 %>% mutate(value = 1) %>% rbind(temp_df)
    
  return(results)
}

germline_all_B <- mediator_function("germline_variant", all_data_B) %>% mutate("Study" = "SchildkrautB")
germline_brca1_all_B <- mediator_function("germline_carrier_BRCA1", all_data_B) %>% mutate("Study" = "SchildkrautB")
germline_brca2_all_B <- mediator_function("germline_carrier_BRCA2", all_data_B) %>% mutate("Study" = "SchildkrautB")
germline_other_all_B <- mediator_function("germline_other", all_data_B) %>% mutate("Study" = "SchildkrautB")

germline_all_W <- mediator_function("germline_variant", all_data_W) %>% mutate("Study" = "SchildkrautW")
germline_brca1_all_W <- mediator_function("germline_carrier_BRCA1", all_data_W) %>% mutate("Study" = "SchildkrautW")
germline_brca2_all_W <- mediator_function("germline_carrier_BRCA2", all_data_W) %>% mutate("Study" = "SchildkrautW")
germline_other_all_W <- mediator_function("germline_other", all_data_W) %>% mutate("Study" = "SchildkrautW")

somatic_B <- mediator_function("somatic_mutation", all_data_B) %>% mutate("Study" = "SchildkrautB")
somatic_brca1_B <- mediator_function("somatic_BRCA1", all_data_B) %>% mutate("Study" = "SchildkrautB")
somatic_brca2_B <- mediator_function("somatic_BRCA2", all_data_B) %>% mutate("Study" = "SchildkrautB")
somatic_other_B <- mediator_function("somatic_other", all_data_B) %>% mutate("Study" = "SchildkrautB")

somatic_W <- mediator_function("somatic_mutation", all_data_W) %>% mutate("Study" = "SchildkrautW")
somatic_brca1_W <- mediator_function("somatic_BRCA1", all_data_W) %>% mutate("Study" = "SchildkrautW")
somatic_brca2_W <- mediator_function("somatic_BRCA2", all_data_W) %>% mutate("Study" = "SchildkrautW")
somatic_other_W <- mediator_function("somatic_other", all_data_W) %>% mutate("Study" = "SchildkrautW")

hrdsig_B <- mediator_function("SBS3", all_data_B) %>% mutate("Study" = "SchildkrautB")
hrdsig_W <- mediator_function("SBS3", all_data_W) %>% mutate("Study" = "SchildkrautW")

anyhrd_B <- mediator_function("any_HRD", all_data_B) %>% mutate("Study" = "SchildkrautB")
anyhrd_W <- mediator_function("any_HRD", all_data_W) %>% mutate("Study" = "SchildkrautW")

bind_B <- rbind(germline_all_B, germline_brca1_all_B, germline_brca2_all_B, germline_other_all_B, somatic_B , somatic_brca1_B, somatic_brca2_B, somatic_other_B, anyhrd_B, hrdsig_B)
bind_W <- rbind(germline_all_W, germline_brca1_all_W, germline_brca2_all_W, germline_other_all_W, somatic_W, somatic_brca1_W, somatic_brca2_W, somatic_other_W, hrdsig_W, anyhrd_W)

bind_B <- bind_B %>% select(term, value, RR_age, RR_tmb)
colnames(bind_B) <- paste0(colnames(bind_B), "_B")
bind_B <- bind_B %>% dplyr::rename(value = value_B, 
                                   term = term_B)

bind_W <- bind_W %>% select(term, value, RR_age, RR_tmb)
colnames(bind_W) <- paste0(colnames(bind_W), "_W")
bind_W <- bind_W %>% dplyr::rename(value = value_W, 
                                   term = term_W)

merge_final <- merge(bind_B, bind_W) %>% select(term, value, RR_age_B, RR_age_W, RR_tmb_B, RR_tmb_W)
merge_final <- merge_final %>% 
  arrange(factor(term, levels = c("germline_carrier_BRCA1", 
                                  "germline_carrier_BRCA2", 
                                  "germline_other",
                                  "germline_variant",
                                  "somatic_BRCA1", 
                                  "somatic_BRCA2", 
                                  "somatic_other",
                                  "somatic_mutation",
                                  "SBS3", 
                                  "any_HRD")))


# Sensitivity including only VUS
germline_all_B_VUS <- mediator_function("germline_variant", all_data_VUS_B) %>% mutate("Study" = "SchildkrautB")
germline_brca1_all_B_VUS <- mediator_function("germline_carrier_BRCA1", all_data_VUS_B) %>% mutate("Study" = "SchildkrautB")
germline_brca2_all_B_VUS <- mediator_function("germline_carrier_BRCA2", all_data_VUS_B) %>% mutate("Study" = "SchildkrautB")
germline_other_all_B_VUS <- mediator_function("germline_other", all_data_VUS_B) %>% mutate("Study" = "SchildkrautB")

germline_all_W_VUS <- mediator_function("germline_variant", all_data_VUS_W) %>% mutate("Study" = "SchildkrautW")
germline_brca1_all_W_VUS <- mediator_function("germline_carrier_BRCA1", all_data_VUS_W) %>% mutate("Study" = "SchildkrautW")
germline_brca2_all_W_VUS <- mediator_function("germline_carrier_BRCA2", all_data_VUS_W) %>% mutate("Study" = "SchildkrautW")
germline_other_all_W_VUS <- mediator_function("germline_other", all_data_VUS_W) %>% mutate("Study" = "SchildkrautW")

somatic_B_VUS <- mediator_function("somatic_mutation", all_data_VUS_B) %>% mutate("Study" = "SchildkrautB")
somatic_brca1_B_VUS <- mediator_function("somatic_BRCA1", all_data_VUS_B) %>% mutate("Study" = "SchildkrautB")
somatic_brca2_B_VUS <- mediator_function("somatic_BRCA2", all_data_VUS_B) %>% mutate("Study" = "SchildkrautB")
somatic_other_B_VUS <- mediator_function("somatic_other", all_data_VUS_B) %>% mutate("Study" = "SchildkrautB")

somatic_W_VUS <- mediator_function("somatic_mutation", all_data_VUS_W) %>% mutate("Study" = "SchildkrautW")
somatic_brca1_W_VUS <- mediator_function("somatic_BRCA1", all_data_VUS_W) %>% mutate("Study" = "SchildkrautW")
somatic_brca2_W_VUS <- mediator_function("somatic_BRCA2", all_data_VUS_W) %>% mutate("Study" = "SchildkrautW")
somatic_other_W_VUS <- mediator_function("somatic_other", all_data_VUS_W) %>% mutate("Study" = "SchildkrautW")

anyhrd_B_VUS <- mediator_function("any_HRD", all_data_VUS_B) %>% mutate("Study" = "SchildkrautB")
anyhrd_W_VUS <- mediator_function("any_HRD", all_data_VUS_W) %>% mutate("Study" = "SchildkrautW")

bind_B_VUS <- rbind(germline_all_B_VUS, germline_brca1_all_B_VUS, germline_brca2_all_B_VUS, germline_other_all_B_VUS, somatic_B_VUS , somatic_brca1_B_VUS, somatic_brca2_B_VUS, somatic_other_B_VUS, anyhrd_B_VUS)
bind_W_VUS <- rbind(germline_all_W_VUS, germline_brca1_all_W_VUS, germline_brca2_all_W_VUS, germline_other_all_W_VUS, somatic_W_VUS, somatic_brca1_W_VUS, somatic_brca2_W_VUS, somatic_other_W_VUS, anyhrd_W_VUS)

bind_B_VUS <- bind_B_VUS %>% select(term, value, RR_age, RR_tmb)
colnames(bind_B_VUS) <- paste0(colnames(bind_B_VUS), "_B")
bind_B_VUS <- bind_B_VUS %>% dplyr::rename(value = value_B, 
                                   term = term_B)

bind_W_VUS <- bind_W_VUS %>% select(term, value, RR_age, RR_tmb)
colnames(bind_W_VUS) <- paste0(colnames(bind_W_VUS), "_W")
bind_W_VUS <- bind_W_VUS %>% dplyr::rename(value = value_W, 
                                   term = term_W)

merge_final_VUS <- merge(bind_B_VUS, bind_W_VUS) %>% select(term, value, RR_age_B, RR_age_W, RR_tmb_B, RR_tmb_W)
merge_final_VUS <- merge_final_VUS %>% 
  arrange(factor(term, levels = c("germline_carrier_BRCA1", 
                                  "germline_carrier_BRCA2", 
                                  "germline_other",
                                  "germline_variant",
                                  "somatic_BRCA1", 
                                  "somatic_BRCA2", 
                                  "somatic_other",
                                  "somatic_mutation",
                                  "any_HRD")))

# Sensitivity including expanded gene list
germline_all_B_expandedlist <- mediator_function("germline_variant_all_genes", all_data_B) %>% mutate("Study" = "SchildkrautB")
germline_all_W_expandedlist <- mediator_function("germline_variant_all_genes", all_data_W) %>% mutate("Study" = "SchildkrautW")
somatic_B_expandedlist <- mediator_function("somatic_mutation_all_genes", all_data_B) %>% mutate("Study" = "SchildkrautB")
somatic_W_expandedlist <- mediator_function("somatic_mutation_all_genes", all_data_W) %>% mutate("Study" = "SchildkrautW")

bind_B_expandedlist <- rbind(germline_all_B_expandedlist, somatic_B_expandedlist)
bind_W_expandedlist <- rbind(germline_all_W_expandedlist, somatic_W_expandedlist)

bind_B_expandedlist <- bind_B_expandedlist %>% select(term, value, RR_age, RR_tmb)
colnames(bind_B_expandedlist) <- paste0(colnames(bind_B_expandedlist), "_B")
bind_B_expandedlist <- bind_B_expandedlist %>% dplyr::rename(value = value_B, 
                                                             term = term_B)

bind_W_expandedlist <- bind_W_expandedlist %>% select(term, value, RR_age, RR_tmb)
colnames(bind_W_expandedlist) <- paste0(colnames(bind_W_expandedlist), "_W")
bind_W_expandedlist <- bind_W_expandedlist %>% dplyr::rename(value = value_W, 
                                                             term = term_W)

merge_final_expanded_list <- merge(bind_B_expandedlist, bind_W_expandedlist) %>% select(term, value, RR_age_B, RR_age_W, RR_tmb_B, RR_tmb_W)


