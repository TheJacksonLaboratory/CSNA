# Data Elements for the QTL Viewer
# 
# The QTL Viewer utilizes R and several different libraries in order to calculate the data for various types of QTL projects.  The following sections will explain each element in detail.
# 
# Please note that some data element must be pre-computed.
# 
# 
# RData Environment Overview
# 
# The following elements should be contained within the RData file.
# 
# ensembl.version – the numerical version of Ensembl
# genoprobs – the genotype probabilities
# K – the kinship matric
# map – list of one element per chromosome, with the genomic position of each marker
# markers – marker names and positions
# 
# The following element is a special element and there must be at least one per RData file.
# 
# dataset.*  - where * should be a very short, unique and informative name.  This element will contain most of the data and will be detailed in the section below.
# The dataset.* element is a list that should contain the following named elements:
#   
# annot.datatype – annotations, where datatype is one of mrna, protein, or phenotype
# annot.samples – annotation data for the samples
# covar.matrix  – a matrix of covariate data, samples (rows) x covariates (columns)
# covar.info  – information describing the covariates
# data  – either a matrix containing data or a list containing several kinds of data
# datatype  – one of mrna, protein, or phenotype
# display.name  – name of the dataset, for QTL Viewer display purposes
# lod.peaks  – a list of LOD peaks over a certain threshold
library(tidyverse)
library(data.table)
library(parallel)
library(foreach)
library(doParallel)
library(qtl2)

load("DO_str_2016_eQTL.rdata")

ensembl.version <- 90
genoprobs <- DO_str_2016_gprobs

K <- DO_str_2016_kinship

map <- DO_str_2016_cross$pmap

# You need two objects annotation and expr --------------------------------
annotation <- data.frame(gene_id = eQTL_maRt$ensembl_gene_id,
                         symbol  = eQTL_maRt$mgi_symbol,
                         chr     = eQTL_maRt$chromosome_name,
                         start   = eQTL_maRt$start_position,
                         end     = eQTL_maRt$end_position,
                         strand  = eQTL_maRt$strand, stringsAsFactors = F)
annotation$middle <- round((annotation$start + annotation$end)/2)
annotation <- annotation[annotation$chr %in% c(1:19,"X"),]
rownames(annotation) <- annotation$gene_id
annotation$nearest_marker <- NA


#nearest marker
cl <- makeCluster(6)
registerDoParallel(cl)
snps <- GM_snps[GM_snps$chr %in% c(1:19,"X"),]
idx <- foreach(i=1:nrow(annotation), .combine='c') %dopar% {
  dist.to <- abs(10^6*snps$pos - annotation$middle[i])
  min.dist <- min(dist.to[snps$chr == annotation$chr[i]])
  which(snps$chr == annotation$chr[i] & dist.to==min.dist)[1]
}
annotation$nearest_marker <- snps[idx,"marker"]
stopCluster(cl)

#dataset object
expr.0 <- DO_str_2016_cross$pheno
expr <- expr.0 [,colnames(expr.0) %in% annotation$gene_id]

#annot.mrna
annot.mrna <- as_tibble(annotation[colnames(expr),])

#annot.samples
covar <- DO_str_2016_cross$covar
rownames(covar) <- covar$subject
covar$sex[covar$sex == "M"] <- 1
covar$sex[covar$sex == "F"] <- 0
annot.samples <- data.frame(mouse.id = rownames(covar), covar[,2:4], stringsAsFactors = F)
rownames(annot.samples) <- annot.samples$mouse.id
annot.samples <- as_tibble(annot.samples)

#covar.matrix
covar.matrix <- model.matrix(~ngen + sex, covar)[, -1]  #covar[,2:4]

#covar.info
covar.info <- as_tibble(data.frame(sample.column   = colnames(covar[,2:3]),
                                   display.name    = c("Sex", "Generation"),
                                   interactive     = rep(FALSE, 2),
                                   primary         = c(TRUE, TRUE),
                                   lod.peaks       = c(NA,NA), stringsAsFactors = FALSE))

#data
data <- expr # you need the expression profile matrix

#datatype
datatype <- "mRNA"

#display.name
display.name <- "DrusgNaiveStriatum_DO"

#lod.peaks
#load cis eqtl
b <- cis <- list()
for(i in c(1:19,"X")){
  b[[i]] <- readRDS(paste0("DO_cis_eQTL_chr",i,".rds")) # this is an example. do it in a loop
  b[[i]] <- b[[i]][,apply(b[[i]], 2, function(x) sum(x >= 7)) > 0]
  b[[i]] <- b[[i]][apply(b[[i]], 1, function(x) sum(x >= 7)) > 0,]
  cis[[i]] <- data.frame(marker=rownames(b[[i]])[row(b[[i]])], 
                    gene=colnames(b[[i]])[col(b[[i]])],
                    lod=c(b[[i]]))
  cis[[i]] <- setDT(cis[[i]])[, .SD[which.max(lod)], by=gene]
  cis[[i]] <- cis[[i]][cis[[i]]$lod >= 7,] # define a cutoff peak
}
cis <- do.call(rbind, cis)

#lod.peaks
additive <- data.frame(annot.id = cis$gene, 
                       marker.id = cis$marker,
                       lod =cis$lod,stringsAsFactors = F)
rownames(additive) <- additive$annot.id
additive <- as_tibble(additive[annot.mrna$gene_id,])
lod.peaks <- list(additive = additive)

#dataset list
assign(paste0("dataset.","DO_Striatum_416"),list(annot.mrna    = annot.mrna,
                                            annot.samples      = annot.samples,
                                            covar.matrix       = covar.matrix,
                                            covar.info         = covar.info, 
                                            data               = data,      
                                            datatype           = datatype, 
                                            display.name       = display.name, 
                                            lod.peaks          = lod.peaks))

#markers
markers <- NULL
for(i in 1:20){
  info <- as.data.frame(DO_str_2016_cross$pmap[i])
  chrom <- rep(i, dim(info)[1])
  marker <- row.names(info)
  pos <- info[,1]
  chr <- cbind.data.frame(marker,chrom,pos)
  markers <- rbind.data.frame(markers, chr)
}
markers[markers$chrom=="20",2] <- "X"


dataset.DO_Striatum_416$annot.mrna <-
  dataset.DO_Striatum_416$annot.mrna %>%
  dplyr::rename(gene.id = gene_id,
                nearest.marker.id = nearest_marker)



dataset.DO_Striatum_416$lod.peaks$additive <-
  dataset.DO_Striatum_416$lod.peaks$additive %>%
  dplyr::rename(gene.id = annot.id)



markers <-
  markers %>%
  as_tibble() %>%
  dplyr::rename(marker.id = marker,
                chr = chrom)

for(i in 1:20){
  temp.genoprobs <- genoprobs[[i]]
  ids <- as.character(unlist(strsplit(row.names(temp.genoprobs), split = "_")))
  ids.1 <- ids[-grep(ids, pattern = "file")]
  ids.2 <- ids.1[-grep(ids.1, pattern = "200735410007|R09C01|200840140078|R05C02")]
  row.names(genoprobs[[i]]) <- ids.2
}


for(i in 1:20){
  temp.k <- K[[i]]
  ids <- as.character(unlist(strsplit(row.names(temp.k), split = "_")))
  ids.1 <- ids[-grep(ids, pattern = "file")]
  ids.2 <- ids.1[-grep(ids.1, pattern = "200735410007|R09C01|200840140078|R05C02")]
  row.names(temp.k) <- ids.2
  colnames(temp.k) <- ids.2
  K[[i]] <- temp.k
}


ids <- as.character(unlist(strsplit(row.names(dataset.DO_Striatum_416$data), split = "_")))
ids.1 <- ids[-grep(ids, pattern = "file")]
ids.2 <- ids.1[-grep(ids.1, pattern = "200735410007|R09C01|200840140078|R05C02")]
row.names(dataset.DO_Striatum_416$data) <- ids.2


genoprobs.36 <- genoprobs
genoprobs <- genoprob_to_alleleprob(genoprobs.36)

for (i in 1:20){
  temp.k <- K[[i]]
  row.names(temp.k) <- row.names(dataset.DO_Striatum_416$covar.matrix)
  colnames(temp.k) <- row.names(dataset.DO_Striatum_416$covar.matrix)
  K[[i]] <- temp.k
}

for (i in 1:20){
  temp.genoprobs <- genoprobs[[i]]
  row.names(temp.genoprobs) <- row.names(dataset.DO_Striatum_416$covar.matrix)
  genoprobs[[i]] <- temp.genoprobs
}



save(genoprobs, 
     K, 
     map, 
     markers, 
     dataset.DO_Striatum_416, 
     file = paste0("qtlviewer_DO_Striatum_416_02072020",".rdata")
)