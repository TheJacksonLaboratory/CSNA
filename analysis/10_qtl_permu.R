####################################################################################################################
#   This script performs permutation qtl mapping for DO mice
#   It takes  args in this script
#   -phenocsv= location of phenotype .csv file
#   -pheno.idx= column number of phenotype for permutation
#   -model= choose a model name: m1, m2, m3
#   -permu.id= number id for permutation
#   -outdir= output directory ending with "/"
#
#   output is one RData file. It contains
#   *.permu.out

#   It will be used for cutoff and visualization.
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

# -permu.id=
str.look <- "-permu.id="
i.look <- c(1:length(args))[substr(args,1,nchar(str.look)) == str.look]
permu.id <- as.character(substr(args[i.look], nchar(str.look)+1, nchar(args[i.look])))
permu.id

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

# one of three models
if(model == "m1"){ # ~Sex+Generation
  #permu
  m1.permu.out <- list()
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
    m1.permu.out[[i]] <- scan1perm(genoprobs  =  apr.i,
                                 pheno        =  phe.i,
                                 kinship      =  k.i,
                                 addcovar     =  addcovar.i,
                                 n_perm       =  100,
                                 cores        =  20)
  }
  save(m1.permu.out,
       file = paste0(outdir,"m1_", name.phenocsv, "_", pheno.list, "_", permu.id, ".RData"))
  
}else if(model == "m2"){ # ~Sex
  #permu
  m2.permu.out <- list()
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
    m2.permu.out[[i]] <- scan1perm(genoprobs    =  apr.i,
                                   pheno        =  phe.i,
                                   kinship      =  k.i,
                                   addcovar     =  addcovar.i,
                                   n_perm       =  100,
                                   cores        =  20)
  }
  save(m2.permu.out,
       file = paste0(outdir,"m2_", name.phenocsv, "_", pheno.list, "_", permu.id, ".RData"))
  
}else{ #m3 model ~Sex and intcovar
  #permu
  m3.permu.out <- list()
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
    m3.permu.out[[i]] <- scan1perm(genoprobs    =  apr.i,
                                   pheno        =  phe.i,
                                   kinship      =  k.i,
                                   addcovar     =  addcovar.i,
                                   intcovar     =  addcovar.i,
                                   n_perm       =  100,
                                   cores        =  20)
  }
  save(m3.permu.out,
       file = paste0(outdir,"m3_", name.phenocsv, "_", pheno.list, "_", permu.id, ".RData"))
}
