library(dplyr)
library(circlize)
library(ComplexHeatmap)

BC.stat <- readRDS("data/BC.stat.rds")
CC.stat <- readRDS("data/CC.stat.rds")
OF.stat <- readRDS("data/OF.stat.rds")
head(OF.stat)
full.stat <- BC.stat %>%
  rbind(CC.stat) %>%
  rbind(OF.stat)
# marker_symbol, z_score, proc_param_name


### Figure 1: Heatmap of Z-scores

## heatmap of gene - phenotype pairs (red: tested, white: untested)
z.freq <- table(full.stat$proc_param_name, full.stat$marker_symbol)
dim(z.freq)
z.id.mat <- ifelse(z.freq==0, 0, 1)

ht = Heatmap(z.id.mat,
             cluster_rows = T, clustering_distance_rows ="binary",
             cluster_columns = T, clustering_distance_columns = "binary",
             show_row_dend = F, show_column_dend = F,  # do not show dendrogram
             show_column_names = F, show_row_names = F, col = c("gray90","red"),
             heatmap_legend_param = list(title_gp=gpar(fontsize=18), labels_gp=gpar(fontsize=15), labels=c("Yes", "No")), name="Missing")
png(file="analysis/Fig1.png", width=600, height=350)
draw(ht)
dev.off()
