####################################################################################################################
#   This script performs first genoprobs, maxmarg, crossover, errorlod, snpg for the qc diagnositics
#
#   Author: Hao He
#   Date:   02/13/2019
#   E-mails: hao.he@jax.org
####################################################################################################################

library(qtl2)

# # Generate json file for all batches ------------------------------------
setwd("/projects/heh/csna_workflow/")
#load json file for the nine batches
gm <- read_cross2("data/Jackson_Lab_Bubier_MURGIGV01/gm.json")

gm
#Let’s omit markers without any genotype data and re-code the genotypes 
#so that the first allele is that which is most common among the founders.
gm <- drop_nullmarkers(gm)
#Dropping 140 markers with no data
for(chr in seq_along(gm$founder_geno)) {
  fg <- gm$founder_geno[[chr]]
  g <- gm$geno[[chr]]
  f1 <- colSums(fg==1)/colSums(fg != 0)
  
  fg[fg==0] <- NA
  g[g==0] <- NA
  
  fg[,f1 < 0.5] <- 4 - fg[,f1 < 0.5]
  g[,f1 < 0.5]  <- 4 - g[,f1 < 0.5]
  
  fg[is.na(fg)] <- 0
  g[is.na(g)] <- 0
  
  gm$founder_geno[[chr]] <- fg
  gm$geno[[chr]] <- g
}

gm

save(gm, file = "data/Jackson_Lab_Bubier_MURGIGV01/gm_2514.RData")

# pr,m,nxo calculation ----------------------------------------------------
pr <- calc_genoprob(gm, cores=20)
m <- maxmarg(pr, cores=20)
nxo <- count_xo(m, cores=20)
save(pr, file = "data/Jackson_Lab_Bubier_MURGIGV01/pr.RData")
str(pr)
save(m, file = "data/Jackson_Lab_Bubier_MURGIGV01/m.RData")
save(nxo, file = "data/Jackson_Lab_Bubier_MURGIGV01/nxo.RData")

e <- calc_errorlod(gm, pr, cores=20)
str(e)
e <- do.call("cbind", e)
save(e, file = "data/Jackson_Lab_Bubier_MURGIGV01/e.RData")

snpg <- predict_snpgeno(gm, m, cores=20)
snpg <- do.call("cbind", snpg)
save(snpg, file = "data/Jackson_Lab_Bubier_MURGIGV01/snpg.RData")
