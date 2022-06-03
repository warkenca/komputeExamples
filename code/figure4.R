library(patchwork)
library(ggpubr)

load("docs/figure/figures.Rmd/sim_BC.rdata")
load("docs/figure/figures.Rmd/sim_CC.rdata")
load("docs/figure/figures.Rmd/sim_OF.rdata")



supp_fig3 <- ggarrange(fig3.1, fig3.2, fig3.3, nrow=1)
supp_fig3

pdf(file="docs/figure/figures.Rmd/sim_all_supp.pdf", width=12, height=4)
supp_fig3
dev.off()
