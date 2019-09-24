####################################################################################################################
#   This script performs scan1blup for DO mice
#   It takes  args in this script
#   -phenocsv= location of phenotype .csv file
#   -pheno.idx= column number of phenotype in the csv
#   -model= choose a model name: m1, m2, m3
#   -outdir= output directory ending with "/"
#
#   output is one RData file. It contains
#
#   It will be used for visualization.
#   Author: Hao He
#   Date:   02/13/2019
#   E-mails: hao.he@jax.org
####################################################################################################################
# Library -----------------------------------------------------------------
# Load packages
library(qtl2)
library(tidyr)
library(dplyr)
library(data.table)
library(foreach)
library(doParallel)
library(parallel)
library(abind)
library(qtl2convert)
library(DOQTL)
library(gap)
library(regress)
library(ggplot2)
library(corrplot)
library(lubridate)
library(broman)
library(qtlcharts)
library(pheatmap)
library(lme4)
options(stringsAsFactors = FALSE)

# args --------------------------------------------------------------------
options(warn=1)
args <- commandArgs()
args

# -phenocsv=
str.look <- "-phenocsv="
i.look <- c(1:length(args))[substr(args,1,nchar(str.look)) == str.look]
phenocsv <- as.character(substr(args[i.look], nchar(str.look)+1, nchar(args[i.look])))
phenocsv
name.phenocsv <- gsub(".*[/]([^.]+)[.].*", "\\1", phenocsv)
name.phenocsv

# -pheno.idx=
str.look <- "-pheno.idx="
i.look <- c(1:length(args))[substr(args,1,nchar(str.look)) == str.look]
pheno.idx <- as.numeric(substr(args[i.look], nchar(str.look)+1, nchar(args[i.look])))
pheno.idx

# -model=
str.look <- "-model="
i.look <- c(1:length(args))[substr(args,1,nchar(str.look)) == str.look]
model <- as.character(substr(args[i.look], nchar(str.look)+1, nchar(args[i.look])))
model

# -outdir=
str.look <- "-outdir="
i.look <- c(1:length(args))[substr(args,1,nchar(str.look)) == str.look]
outdir <- as.character(substr(args[i.look], nchar(str.look)+1, nchar(args[i.look])))
outdir

setwd("/projects/heh/csna_workflow/")
# pheno data --------------------------------------------------------------
pheno <- read.csv(phenocsv, header = TRUE)
pheno$Sex[pheno$Sex == "Male"] <- 1
pheno$Sex[pheno$Sex == "Female"] <- 0

# geno data ---------------------------------------------------------------
load("data/Jackson_Lab_Bubier_MURGIGV01/gm_DO2437_qc.RData")
load("data/Jackson_Lab_Bubier_MURGIGV01/apr_DO2437.RData")
#new geno id
geno.id <- as.character(do.call(rbind.data.frame, strsplit(ind_ids(gm_DO2437_qc), "_"))[,6])

#replace new ids in apr
nodup.geno.id <- data.frame(geno.id = geno.id,
                            nodup.geno.id = NA)
nodup.geno.id <- nodup.geno.id %>% mutate(nodup.geno.id = make.unique(as.character(geno.id)))
nodup.geno.id <- nodup.geno.id$nodup.geno.id
names(nodup.geno.id) <- ind_ids(gm_DO2437_qc)
apr <- replace_ids(apr,nodup.geno.id)

#kinship
k = calc_kinship(probs = apr, type = "loco", use_allele_probs = TRUE, cores = 20)

# permutation for qtlmapping  -------------------------------------------------------------
#pheno list
pheno.list <- colnames(pheno)[pheno.idx]
pheno.list

# one of three models
if(model == "m1"){ # ~Sex+Generation
  #blup
  m1.blup <- list()
  load(paste0("/projects/heh/csna_workflow/output/m1.",name.phenocsv,".RData"))
  for(i in pheno.list){
    print(i)
    overlap.id <- intersect(nodup.geno.id, pheno[!is.na(pheno[,i]), "Mouse.ID"])
    pheno.i <- pheno[pheno$Mouse.ID %in% overlap.id,]
    rownames(pheno.i) <- pheno.i$Mouse.ID
    phe.i = pheno.i[,i]
    names(phe.i) <- pheno.i[,"Mouse.ID"]
    addcovar.i = model.matrix(~Sex+Generation, data = pheno.i)[,-1]
    rownames(addcovar.i) <- pheno.i$Mouse.ID
    #subset apr
    apr.i <- apr[as.character(overlap.id),]
    #subset k
    k.i <- lapply(k, function(y) y[dimnames(y)[[1]] %in% overlap.id, dimnames(y)[[2]] %in% overlap.id])
    #subset chr
    chr <- as.character(m1.peak[[i]][m1.peak[[i]]$lodcolumn == i,"chr"])
    print(chr)
    m1.blup[[i]] <- scan1blup(genoprobs    =  apr.i[,chr],
                              pheno        =  phe.i,
                              kinship      =  k.i[[chr]],
                              addcovar     =  addcovar.i,
                              cores        =  20)
  }
  save(m1.blup,
       file = paste0(outdir,"m1_blup_", name.phenocsv, "_", pheno.list, ".RData"))
  
}else if(model == "m2"){ # ~Sex
  #blup
  m2.blup <- list()
  load(paste0("/projects/heh/csna_workflow/output/m2.",name.phenocsv,".RData"))
  for(i in pheno.list){
    print(i)
    overlap.id <- intersect(nodup.geno.id, pheno[!is.na(pheno[,i]), "Mouse.ID"])
    pheno.i <- pheno[pheno$Mouse.ID %in% overlap.id,]
    rownames(pheno.i) <- pheno.i$Mouse.ID
    phe.i = pheno.i[,i]
    names(phe.i) <- pheno.i[,"Mouse.ID"]
    addcovar.i = model.matrix(~Sex, data = pheno.i)[,-1,drop=FALSE]
    rownames(addcovar.i) <- pheno.i$Mouse.ID
    #subset apr
    apr.i <- apr[as.character(overlap.id),]
    #subset k
    k.i <- lapply(k, function(y) y[dimnames(y)[[1]] %in% overlap.id, dimnames(y)[[2]] %in% overlap.id])
    #subset chr
    chr <- as.character(m2.peak[[i]][m2.peak[[i]]$lodcolumn == i,"chr"])
    print(chr)
    m2.blup[[i]] <- scan1blup(genoprobs    =  apr.i[,chr],
                              pheno        =  phe.i,
                              kinship      =  k.i[[chr]],
                              addcovar     =  addcovar.i,
                              cores        =  20)
  }
  save(m2.blup,
       file = paste0(outdir,"m2_blup_", name.phenocsv, "_", pheno.list, ".RData"))
  
}else{ #m3 model ~Sex and intcovar
  #blup
  m3.blup <- list()
  load(paste0("/projects/heh/csna_workflow/output/m3.",name.phenocsv,".RData"))
  for(i in pheno.list){
    print(i)
    overlap.id <- intersect(nodup.geno.id, pheno[!is.na(pheno[,i]), "Mouse.ID"])
    pheno.i <- pheno[pheno$Mouse.ID %in% overlap.id,]
    rownames(pheno.i) <- pheno.i$Mouse.ID
    phe.i = pheno.i[,i]
    names(phe.i) <- pheno.i[,"Mouse.ID"]
    addcovar.i = model.matrix(~Sex, data = pheno.i)[,-1,drop=FALSE]
    rownames(addcovar.i) <- pheno.i$Mouse.ID
    #subset apr
    apr.i <- apr[as.character(overlap.id),]
    #subset k
    k.i <- lapply(k, function(y) y[dimnames(y)[[1]] %in% overlap.id, dimnames(y)[[2]] %in% overlap.id])
    #subset chr
    chr <- as.character(m3.peak[[i]][m3.peak[[i]]$lodcolumn == i,"chr"])
    print(chr)
    m3.blup[[i]] <- scan1blup(genoprobs    =  apr.i[,chr],
                              pheno        =  phe.i,
                              kinship      =  k.i[[chr]],
                              addcovar     =  addcovar.i,
                              cores        =  20)
  }
  save(m3.blup,
       file = paste0(outdir,"m3_blup_", name.phenocsv, "_", pheno.list,".RData"))
}
