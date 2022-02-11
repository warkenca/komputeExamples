library(dplyr)
library(ComplexHeatmap)

############################################
## Preprocess IMPC Control Phenotype data ##
############################################

### Read all control phenotype data
load(file="~/Google Drive Miami/Miami_IMPC/data/v10.1/AllControls_small.Rdata")
#load(file="G:/.shortcut-targets-by-id/1SeBOMb4GZ2Gkldxp4QNEnFWHOiAqtRTz/Miami_IMPC/data/v10.1/AllControls_small.Rdata")
dim(allpheno)
head(allpheno)
allpheno.org <- allpheno
#allpheno <- allpheno.org

### Correct procedure and phenotype names, filter out time series data
allpheno = allpheno %>%
  filter(procedure_name=="Body Composition (DEXA lean/fat)") %>%
  mutate(proc_short_name=recode(procedure_name, "Body Composition (DEXA lean/fat)"="BC")) %>%
  #mutate(parameter_name=recode(parameter_name, "Triglyceride"="Triglycerides")) %>%
  mutate(proc_param_name=paste0(proc_short_name,"_",parameter_name)) %>%
  mutate(proc_param_name_stable_id=paste0(proc_short_name,"_",parameter_name,"_",parameter_stable_id))

## Extract time series data and find out parameter names
ts <- allpheno %>% filter(observation_type=="time_series")
table(ts$proc_param_name)
# Filter out time series data
allpheno <- allpheno %>% filter(observation_type!="time_series")
table(allpheno$proc_param_name)
dim(allpheno)

BC.data <- allpheno[,c("biological_sample_id","procedure_name","parameter_name","data_point","sex","phenotyping_center","strain_name","proc_short_name",
                       "proc_param_name","proc_param_name_stable_id")]
dim(BC.data)
saveRDS(BC.data, file="data/BC.data.rds")




#######################################
## Preprocess IMPC summary stat data ##
#######################################

## Read KOMPv10.1
KOMPv10.1.file = "~/Google Drive Miami/Miami_IMPC/data/v10.1/IMPC_ALL_statistical_results.csv.gz"
#KOMPv10.1.file = "G:/.shortcut-targets-by-id/1SeBOMb4GZ2Gkldxp4QNEnFWHOiAqtRTz/Miami_IMPC/data/v10.1/IMPC_ALL_statistical_results.csv.gz"
KOMPv10.1 = fread(KOMPv10.1.file, header=TRUE, sep=",")
KOMPv10.1$parameter_name <- trimws(KOMPv10.1$parameter_name) #remove white spaces
KOMPv10.1$proc_param_name <- paste0(KOMPv10.1$procedure_name,"_",KOMPv10.1$parameter_name)

table(KOMPv10.1$procedure_name, KOMPv10.1$data_type)
#dat <- KOMPv10.1 %>% select(procedure_name=="Gross Pathology and Tissue Collection")

## extract unidimensional data only.
dim(KOMPv10.1)
KOMPv10.1 <- KOMPv10.1 %>% filter(data_type=="unidimensional")
dim(KOMPv10.1)

## count the number of tests in each phenotype
proc.list <- table(KOMPv10.1$procedure_name)
#proc.list <- proc.list[proc.list>1000]
proc.list
length(proc.list)

pheno.list <- table(KOMPv10.1$proc_param_name)
pheno.list <- pheno.list[pheno.list>1000] # find list of phenotypes with more than 1000 tests (i.e. 1000 mutants tested)
pheno.list <- names(pheno.list)
pheno.list
length(pheno.list) #122

# Use phenotypes with more than 1000 tests (i.e. 1000 mutants tested)
dim(KOMPv10.1)
all.stat <- KOMPv10.1 %>% filter(proc_param_name %in% pheno.list)
dim(all.stat)

mtest <- table(all.stat$proc_param_name, all.stat$marker_symbol)
mtest <-as.data.frame.matrix(mtest)
dim(mtest)

if(FALSE){
  nmax <-max(mtest)
  library(circlize)
  col_fun = colorRamp2(c(0, nmax), c("white", "red"))
  col_fun(seq(0, nmax))
  ht = Heatmap(as.matrix(mtest), cluster_rows = FALSE, cluster_columns = FALSE, show_column_names = F, col = col_fun,
               row_names_gp = gpar(fontsize = 5), name="Count")
  draw(ht)
}

table(all.stat$procedure_name)

BC.stat = all.stat %>%
  dplyr::select(phenotyping_center, procedure_name, parameter_name, zygosity, allele_symbol,
                genotype_effect_parameter_estimate, genotype_effect_stderr_estimate,
                genotype_effect_p_value, phenotyping_center, allele_name, marker_symbol) %>%
  filter(procedure_name == "Body Composition (DEXA lean/fat)") %>%
  mutate(procedure_name=recode(procedure_name, "Body Composition (DEXA lean/fat)"="BC")) %>%
  mutate(z_score = genotype_effect_parameter_estimate/genotype_effect_stderr_estimate,
         proc_param_name=paste0(procedure_name,"_",parameter_name),
         gene_pheno = paste0(parameter_name, "_", allele_symbol))

table(BC.stat$parameter_name, BC.stat$procedure_name)

BC.stat <- BC.stat[,c("procedure_name","parameter_name","zygosity","allele_symbol","marker_symbol","z_score",
                      "proc_param_name","gene_pheno")]

head(BC.stat)
dim(BC.stat)
saveRDS(BC.stat, file="data/BC.stat.rds")
