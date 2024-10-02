library(tidyr)
library(dplyr)
library(gtsummary)
library(gt)

directory <- "" #set working directory
source(file.path(directory, "code/1.utils.R"))

all_data <- readxl::read_excel(file.path(directory, "/data/all_data.xlsx")) %>% filter(!SUID %in% neo_ids)

#Make table 1
table_1 <- all_data %>%
  tbl_summary(
    by = Study,
    statistic = list(all_categorical() ~ "{n} :({p}%)", all_continuous() ~ "{median} :({p0}, {p100})"),
    digits = list(all_categorical() ~ c(0, 0), 
                  all_continuous() ~ 2),
    include=c(refage_5_year_recode,
              stage_b,
              ethnicity,
              dblkstat_recode,
              famhxov_recode,
              famhxbr_recode, 
              germline_carrier_BRCA1, 
              germline_carrier_BRCA2,
              germline_carrier_RAD51C, 
              germline_carrier_RAD51D, 
              germline_carrier_BRIP1, 
              germline_carrier_PALB2, 
              germline_variant,
              germline_other, 
              somatic_BRCA1, 
              somatic_BRCA2,
              somatic_RAD51C, 
              somatic_RAD51D, 
              somatic_BRIP1, 
              somatic_PALB2, 
              somatic_mutation,
              somatic_other,
              any_HRD_mutation,
              SBS3, 
              any_HRD,
              total_perMB)) %>%
  add_overall() %>%
  add_p(pvalue_fun = label_style_pvalue(digits = 2)) %>%
  bold_labels() %>%
  modify_caption("HRD categories are not mutually exclusive. Some cases had variants in >1 genes. Non-BRCA category includes alterations in RAD51C, RAD51D, BRIP1, and PALB2. P-value calculated using Fisher’s exact test; Pearson’s Chi-squared test; Wilcoxon rank sum test. FIGO: International Federation of Gynecology and Obstetrics; HRD: homologous recombination deficiency; TMB: tumor mutation burden")

#Make supplemental table 3
sensitivity_all_genes <- all_data %>%
  tbl_summary(
    by = Study,
    include=c(any_HRD_all_genes,
              any_HRD_mutation_all_genes,
              germline_variant_all_genes,
              germline_carrier_BRCA1, 
              germline_carrier_BRCA2,  
              germline_other_all_genes,
              hrd_all_germline_cols_other,
              somatic_mutation_all_genes,
              somatic_other_all_genes, 
              somatic_BRCA1, 
              somatic_BRCA2, 
              hrd_all_somatic_cols_other), 
   ) %>% 
  add_overall() %>% 
  add_p()