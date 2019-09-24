####################################################################################################################
#   This script will run the genoprobs calculation after qc and do interpolatino to 69K grid
#
#   Author: Hao He
#   Date:   02/13/2019
#   E-mails: hao.he@jax.org
####################################################################################################################

### Load required library packages
options(stringsAsFactors = FALSE)
library(qtl2) #(0.18 version)
library(foreach)
library(doParallel)
library(parallel)
library(abind)
library(DOQTL)

setwd("/projects/heh/csna_workflow/")
load("data/Jackson_Lab_Bubier_MURGIGV01/gm_DO2437_qc.RData")
#pr <- calc_genoprob(gm_DO2437_qc, cores=20)
#save(pr, file = "data/Jackson_Lab_Bubier_MURGIGV01/pr_DO2437.RData")
#
#apr <- genoprob_to_alleleprob(pr, cores = 20)
#save(apr, file = "data/Jackson_Lab_Bubier_MURGIGV01/apr_DO2437.RData")

load("data/Jackson_Lab_Bubier_MURGIGV01/pr_DO2437.RData")
load("data/Jackson_Lab_Bubier_MURGIGV01/apr_DO2437.RData")

#x chr
if(dim(pr$X)[2] == 44){
  #X chromosome should have 36 states
  pr$X <- pr$X[,(-37:-44),]
}

#Create a large 3D array with samples/states/SNPs in dimensions 1,2,3.
#combine all chrs into one 3d array
pr.3d <- do.call("abind",list(pr,along = 3))

# impute x to 69k grid -----------------------------------------------------------------------
#gigamuga snps
load(url("ftp://ftp.jax.org/MUGA/GM_snps.Rdata"))
#define from
from = GM_snps[intersect(dimnames(pr.3d)[[3]],GM_snps$marker),1:4]
#define to
grid = read.delim("data/marker_grid_0.02cM_plus.txt", stringsAsFactors = F)
rownames(grid) <- grid$marker
to = grid[,1:3]
#subset pr.3d to the marker
sub.pr.3d <- pr.3d[,,from$marker]
#interplote to 69K
#interplote function
interplote_fun <- function(i){
  z <- t(interpolate.markers(data = as.matrix(t(sub.pr.3d[i,,])), from = from, to = to))
  return(z)
}

print(paste0("Number of samples, ",dim(sub.pr.3d)[1]))
system.time({
  pr.69k <- mclapply(X = 1:dim(sub.pr.3d)[1], FUN = interplote_fun, mc.preschedule = F, mc.cores = 10)
})
names(pr.69k) <- dimnames(sub.pr.3d)[[1]]
pr.69k <- do.call("abind",list(pr.69k,along = 3))
pr.69k <- aperm(pr.69k, perm = c(3,1,2))
save(pr.69k, file = "data/Jackson_Lab_Bubier_MURGIGV01/pr_69k_DO2437.RData")

print("pr_69k_DO2437 done")
# condense pr_69k_DO2437 to allele probs -----------------------------------------------------------------------
#split into chr
pr.69kchr <- list()
chr.names <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19")
for (chr in chr.names){
  pr.69kchr[[chr]] <- pr.69k[,,grid[grid$chr == chr,"marker"]]
}
#xchr
xgeno <- array(0, dim = c(dim(pr.69k)[1], 8, length(grid[grid$chr == "X","marker"])))
dimnames(xgeno)[[2]] <- c("AY", "BY", "CY", "DY", "EY", "FY", "GY", "HY")
pr.69kchr$`X` <- abind(pr.69k[,,grid[grid$chr == "X","marker"]],xgeno,along = 2)
#add attributes
attr(pr.69kchr, "crosstype") <- "do"
attr(pr.69kchr, "is_x_chr") <- structure(c(rep(FALSE,19),TRUE), names=c(1:19,"X"))
attr(pr.69kchr, "alleles") <- LETTERS[1:8]
attr(pr.69kchr, "alleleprobs") <- FALSE
attr(pr.69kchr, "class") <- c("calc_genoprob", "list")
save(pr.69kchr, file = "data/Jackson_Lab_Bubier_MURGIGV01/pr_69kchr_DO2437.RData")
#condense
apr.69kchr <- genoprob_to_alleleprob(pr.69kchr,cores = 20)
# Save apr.69kchr to a file
save(apr.69kchr, file = "data/Jackson_Lab_Bubier_MURGIGV01/apr_69kchr_DO2437.RData")
