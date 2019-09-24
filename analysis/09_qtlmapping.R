####################################################################################################################
#   This script performs qtl mapping for DO mice
#   It takes  args in this script
#   -phenocsv= location of phenotype .csv file
#   -firstpheno.idx= column number of first phenotype
#   -lastpheno.idx= column number of last phenotype
#   -model= choose a model name: m1, m2, m3
#   -outdir= output directory ending with "/"
#
#   output is a RData file for corresponding model and phenotype file. It contains
#   *.qtl.out, 
#   *.peak, 
#   *.coef,
#   *.snps,
#   *.genes
#   They will be used for visualization.
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

# -firstpheno.idx=
str.look <- "-firstpheno.idx="
i.look <- c(1:length(args))[substr(args,1,nchar(str.look)) == str.look]
pheno.first <- as.character(substr(args[i.look], nchar(str.look)+1, nchar(args[i.look])))
pheno.first

# -lastpheno.idx=
str.look <- "-lastpheno.idx="
i.look <- c(1:length(args))[substr(args,1,nchar(str.look)) == str.look]
pheno.last <- as.character(substr(args[i.look], nchar(str.look)+1, nchar(args[i.look])))
pheno.last

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

# intersection of geno and pheno ------------------------------------------
#pheno list
pheno.list <- colnames(pheno)[pheno.first:pheno.last]
num.geno.pheno <- data.frame(pheno = pheno.list,
                             num.pheno = 0,
                             num.geno = 0,
                             num.overlap = 0)
for(i in 1:length(pheno.list)){
  num.geno.pheno[i,2] <- length(na.omit(pheno[,pheno.list[i]]))
  num.geno.pheno[i,3] <- length(geno.id)
  num.geno.pheno[i,4] <- length(intersect(geno.id, pheno[!is.na(pheno[,pheno.list[i]]), "Mouse.ID"]))
}
write.csv(num.geno.pheno, file = paste0(outdir, "num.geno.pheno.in.", name.phenocsv,".csv"), row.names = F, quote = F)

# qtlmapping  -------------------------------------------------------------
query_variants <- create_variant_query_func("data/cc_variants.sqlite")
query_genes <- create_gene_query_func("data/mouse_genes_mgi.sqlite")

# one of three models
if(model == "m1"){ # ~Sex+Generation
  #scan
  m1.qtl.out <- list()
  m1.peak    <- list()
  m1.coef    <- list()
  m1.snps    <- list()
  m1.genes   <- list()
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
    m1.qtl.out[[i]] <- scan1(genoprobs    =  apr.i,
                             pheno        =  phe.i,
                             kinship      =  k.i,
                             addcovar     =  addcovar.i,
                             cores        =  20)
    #peak
    m1.peak[[i]] <- find_peaks(m1.qtl.out[[i]], gm_DO2437_qc$gmap,
                            threshold=maxlod(m1.qtl.out[[i]])-0.001,
                            drop = 1.5)
    m1.peak[[i]]$lodcolumn <- i
    print(m1.peak[[i]])
    #effect plot
    #subset chr
    chr <- as.character(m1.peak[[i]][m1.peak[[i]]$lodcolumn == i,"chr"])
    #coeff
    m1.coef[[i]] <- scan1coef(genoprobs    =  apr.i[,chr],
                              pheno        =  phe.i,
                              kinship      =  k[[chr]],
                              addcovar     =  addcovar.i)
    #gene
    variants <- query_variants(chr, 
                               m1.peak[[i]][m1.peak[[i]]$lodcolumn == i,"ci_lo"], 
                               m1.peak[[i]][m1.peak[[i]]$lodcolumn == i,"ci_hi"])
    m1.snps[[i]] <-scan1snps(genoprobs = apr.i, 
                             map = gm_DO2437_qc$pmap, 
                             pheno = phe.i, 
                             kinship = k[[chr]],
                             addcovar = addcovar.i,
                             query_func=query_variants,
                             chr=chr, 
                             start=m1.peak[[i]][m1.peak[[i]]$lodcolumn == i,"ci_lo"], 
                             end=m1.peak[[i]][m1.peak[[i]]$lodcolumn == i,"ci_hi"], 
                             keep_all_snps=TRUE)
    
    m1.genes[[i]] <- query_genes(chr, 
                                 m1.peak[[i]][m1.peak[[i]]$lodcolumn == i,"ci_lo"], 
                                 m1.peak[[i]][m1.peak[[i]]$lodcolumn == i,"ci_hi"])
  }
  save(m1.qtl.out, 
       m1.peak, 
       m1.coef,
       m1.snps,
       m1.genes,
       file = paste0(outdir,"m1.", name.phenocsv, ".RData"))
  
}else if(model == "m2"){ # ~Sex
  #scan
  m2.qtl.out <- list()
  m2.peak    <- list()
  m2.coef    <- list()
  m2.snps    <- list()
  m2.genes   <- list()
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
    m2.qtl.out[[i]] <- scan1(genoprobs    =  apr.i,
                             pheno        =  phe.i,
                             kinship      =  k.i,
                             addcovar     =  addcovar.i,
                             cores        =  20)
    #peak
    m2.peak[[i]] <- find_peaks(m2.qtl.out[[i]], gm_DO2437_qc$gmap,
                               threshold=maxlod(m2.qtl.out[[i]])-0.001,
                               drop = 1.5)
    m2.peak[[i]]$lodcolumn <- i
    print(m2.peak[[i]])
    #effect plot
    #subset chr
    chr <- as.character(m2.peak[[i]][m2.peak[[i]]$lodcolumn == i,"chr"])
    #coeff
    m2.coef[[i]] <- scan1coef(genoprobs    =  apr.i[,chr],
                              pheno        =  phe.i,
                              kinship      =  k[[chr]],
                              addcovar     =  addcovar.i)
    #gene
    variants <- query_variants(chr, 
                               m2.peak[[i]][m2.peak[[i]]$lodcolumn == i,"ci_lo"], 
                               m2.peak[[i]][m2.peak[[i]]$lodcolumn == i,"ci_hi"])
    m2.snps[[i]] <-scan1snps(genoprobs = apr.i, 
                             map = gm_DO2437_qc$pmap, 
                             pheno = phe.i, 
                             kinship = k[[chr]],
                             addcovar = addcovar.i,
                             query_func=query_variants,
                             chr=chr, 
                             start=m2.peak[[i]][m2.peak[[i]]$lodcolumn == i,"ci_lo"], 
                             end=m2.peak[[i]][m2.peak[[i]]$lodcolumn == i,"ci_hi"], 
                             keep_all_snps=TRUE)
    
    m2.genes[[i]] <- query_genes(chr, 
                                 m2.peak[[i]][m2.peak[[i]]$lodcolumn == i,"ci_lo"], 
                                 m2.peak[[i]][m2.peak[[i]]$lodcolumn == i,"ci_hi"])
  }
  save(m2.qtl.out, 
       m2.peak, 
       m2.coef,
       m2.snps,
       m2.genes,
       file = paste0(outdir,"m2.", name.phenocsv, ".RData"))
  
}else{ #m3 model ~Sex and intcovar
  #scan
  m3.qtl.out <- list()
  m3.peak    <- list()
  m3.coef    <- list()
  m3.snps    <- list()
  m3.genes   <- list()
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
    m3.qtl.out[[i]] <- scan1(genoprobs    =  apr.i,
                             pheno        =  phe.i,
                             kinship      =  k.i,
                             addcovar     =  addcovar.i,
                             intcovar     =  addcovar.i,
                             cores        =  20)
    #peak
    m3.peak[[i]] <- find_peaks(m3.qtl.out[[i]], gm_DO2437_qc$gmap,
                               threshold=maxlod(m3.qtl.out[[i]])-0.001,
                               drop = 1.5)
    m3.peak[[i]]$lodcolumn <- i
    print(m3.peak[[i]])
    #effect plot
    #subset chr
    chr <- as.character(m3.peak[[i]][m3.peak[[i]]$lodcolumn == i,"chr"])
    #coeff
    m3.coef[[i]] <- scan1coef(genoprobs    =  apr.i[,chr],
                              pheno        =  phe.i,
                              kinship      =  k[[chr]],
                              addcovar     =  addcovar.i,
                              intcovar     =  addcovar.i)
    #gene
    variants <- query_variants(chr, 
                               m3.peak[[i]][m3.peak[[i]]$lodcolumn == i,"ci_lo"], 
                               m3.peak[[i]][m3.peak[[i]]$lodcolumn == i,"ci_hi"])
    m3.snps[[i]] <-scan1snps(genoprobs = apr.i, 
                             map = gm_DO2437_qc$pmap, 
                             pheno = phe.i, 
                             kinship = k[[chr]],
                             addcovar = addcovar.i,
                             intcovar = addcovar.i,
                             query_func=query_variants,
                             chr=chr, 
                             start=m3.peak[[i]][m3.peak[[i]]$lodcolumn == i,"ci_lo"], 
                             end=m3.peak[[i]][m3.peak[[i]]$lodcolumn == i,"ci_hi"], 
                             keep_all_snps=TRUE)
    
    m3.genes[[i]] <- query_genes(chr, 
                                 m3.peak[[i]][m3.peak[[i]]$lodcolumn == i,"ci_lo"], 
                                 m3.peak[[i]][m3.peak[[i]]$lodcolumn == i,"ci_hi"])
  }
  save(m3.qtl.out, 
       m3.peak, 
       m3.coef,
       m3.snps,
       m3.genes,
       file = paste0(outdir,"m3.", name.phenocsv, ".RData"))
}
