library(dplyr)
library(knitr)
library(kableExtra)

load(file="G:/.shortcut-targets-by-id/1SeBOMb4GZ2Gkldxp4QNEnFWHOiAqtRTz/Miami_IMPC/data/v10.1/AllControls_small.Rdata")
dim(allpheno)
head(allpheno)
phen_summary <- allpheno %>%
  summarize("Mice count" = length(unique(biological_sample_id)),
            "Phenotype domain" = length(unique(procedure_name)),
            "")

head(phen_summary)
unique(allpheno$procedure_name)

examples_summary <- allpheno %>%
  filter(procedure_name %in% c("Body Composition (DEXA lean/fat)", "Clinical Chemistry", "Open Field")) %>%
  group_by(procedure_name) %>%
  summarize("Observation Count" = n(),
            "Parameter Count" = length(unique(parameter_name)))
colnames(examples_summary)[1] <- "Procedure Name"
examples_summary

## Table for summary stats - how many genes are there for each domain
