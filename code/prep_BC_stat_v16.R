library(dplyr)
library(ComplexHeatmap)
library(data.table)

#######################################
## Preprocess IMPC summary stat data ##
#######################################

## Read KOMPv16
KOMPv16.file = "~/Google Drive Miami/Miami_IMPC/data/v16/statistical-results-ALL.csv.gz"
#KOMPv16.file = "G:/.shortcut-targets-by-id/1SeBOMb4GZ2Gkldxp4QNEnFWHOiAqtRTz/Miami_IMPC/data/v16/IMPC_ALL_statistical_results.csv.gz"
KOMPv16 = fread(KOMPv16.file, header=TRUE, sep=",")
KOMPv16$parameter_name <- trimws(KOMPv16$parameter_name) #remove white spaces
KOMPv16$proc_param_name <- paste0(KOMPv16$procedure_name,"_",KOMPv16$parameter_name)

table(KOMPv16$procedure_name, KOMPv16$data_type)
#dat <- KOMPv16 %>% select(procedure_name=="Gross Pathology and Tissue Collection")

## extract unidimensional data only.
dim(KOMPv16)
KOMPv16 <- KOMPv16 %>% filter(data_type=="unidimensional")
dim(KOMPv16)

## count the number of tests in each phenotype
proc.list <- table(KOMPv16$procedure_name)
#proc.list <- proc.list[proc.list>1000]
proc.list
length(proc.list)

pheno.list <- table(KOMPv16$proc_param_name)
pheno.list <- pheno.list[pheno.list>1000] # find list of phenotypes with more than 1000 tests (i.e. 1000 mutants tested)
pheno.list <- names(pheno.list)
pheno.list
length(pheno.list) #152

# Use phenotypes with more than 1000 tests (i.e. 1000 mutants tested)
dim(KOMPv16)
all.stat <- KOMPv16 %>% filter(proc_param_name %in% pheno.list)
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
saveRDS(BC.stat, file="data/BC.stat.v16.rds")
