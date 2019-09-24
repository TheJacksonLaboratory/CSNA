####################################################################################################################
#   This script will save the html files for plotting founder props across all generations in each chr
#
#   Author: Hao He
#   Date:   02/13/2019
#   E-mails: hao.he@jax.org
####################################################################################################################
library(htmlwidgets)
library(plotly)

setwd("/projects/heh/csna_workflow")
load("data/Jackson_Lab_Bubier_MURGIGV01/gm_DO2437_qc.RData")
load("output/prop_across_generation_chr_p.RData")

setwd("/projects/heh/")
for(c in unique(names(gm_DO2437_qc$geno))){
  print(c)
  htmlwidgets::saveWidget(as_widget(p[[c]]), paste0("prop_across_generation_chr",c,".html"))
}
