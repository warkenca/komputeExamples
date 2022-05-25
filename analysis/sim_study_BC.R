library(MASS)

# Source functions and phenotypic correlation matrix for BC domain
source("code/svd_impute.R")
source("code/sim_function.R")
pheno.cor <- readRDS("data/BC_pheno_cor.RDS")

r <- 6
n.samples <- 20000
set.seed(05152022)

# Matrix completion simulations
mc.res20 <- simulation(n.samples=n.samples, pheno.cor=pheno.cor, mask.prop=.2)
mc.res40 <- simulation(n.samples=n.samples, pheno.cor=pheno.cor, mask.prop=.4)
mc.res60 <- simulation(n.samples=n.samples, pheno.cor=pheno.cor, mask.prop=.6)

# KOMPUTE simulations without info cutoff
komp.res20 <- simulation(n.samples=n.samples, pheno.cor=pheno.cor, mask.prop=.2, method="kompute")
komp.res40 <- simulation(n.samples=n.samples, pheno.cor=pheno.cor, mask.prop=.4, method="kompute")
komp.res60 <- simulation(n.samples=n.samples, pheno.cor=pheno.cor, mask.prop=.6, method="kompute")

# KOMPUTE simulations with 0.8 info cutoff
komp.info.res20 <- simulation(n.samples=n.samples, pheno.cor=pheno.cor, mask.prop=.2, method="kompute", info.cutoff = 0.8)
komp.info.res40 <- simulation(n.samples=n.samples, pheno.cor=pheno.cor, mask.prop=.4, method="kompute", info.cutoff = 0.8)
komp.info.res60 <- simulation(n.samples=n.samples, pheno.cor=pheno.cor, mask.prop=.6, method="kompute", info.cutoff = 0.8)

# Create table of correlations
cor.table <- data.frame(c(mc.res20$cor, komp.res20$cor, komp.info.res20$cor),
                        c(mc.res40$cor, komp.res40$cor, komp.info.res40$cor),
                        c(mc.res60$cor, komp.res60$cor, komp.info.res60$cor))
colnames(cor.table) <- c("20% Removed", "40% Removed", "60% Removed")
rownames(cor.table) <- c("MC", "KPT", "KPT w/ info > 0.8")

cor.table

