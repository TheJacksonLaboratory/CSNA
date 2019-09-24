####################################################################################################################
#   This script performs to grab intensities from GeneSeek FinalReport.txt files
#   convert them to a single big data frame, and save for further use.
#
#   Author: Hao He
#   Date:   02/13/2019
#   E-mails: hao.he@jax.org
####################################################################################################################

setwd("/projects/heh/csna_workflow/data/FinalReport/")
# simple version of data.table::fread()
myfread <- function(filename) data.table::fread(filename, data.table=FALSE, skip=9)

# read the data
batch.name <- c(
  "Jackson_Lab_Bubier_MURGIGV01_20160908",   "Jackson_Lab_Bubier_MURGIGV01_20161227",
  "Jackson_Lab_Bubier_MURGIGIV01_20171001",  "Jackson_Lab_Bubier_MURGIGV01_20170904",
  "Jackson_Lab_Bubier_MURGIGV01_20180518",   "Jackson_Lab_Bubier_MURGIGV01_20181206",
  "Jackson_Lab_Bubier_MURGIGV01_20181207",   "Jackson_Lab_Bubier_MURGIGV01_20190108",
  "Jackson_Lab_Bubier_MURGIGV01_20190425")

ifiles <- paste0(batch.name,"/",batch.name,"_FinalReport.txt")
dat <- vector("list", length(ifiles))
for(i in ifiles) {
  print(i)
  dat[[i]] <- myfread(i)
  #sample name
  sample.map <- read.table(file = paste0(gsub( "/.*$", "", i),"/","Sample_Map.txt"), header = TRUE, sep = "\t")
  #check whether ID order match
  cat("---Sample ID order in Finalreport and sample map matched: ",
      all.equal(as.character(do.call(rbind.data.frame, strsplit(as.character(unique(dat[[i]]$`Sample ID`)), "_"))[,1]), 
                as.character(sample.map$ID))
  )
  
  #new sample ID
  dat[[i]][,"Sample ID"] <- rep(paste0(gsub( "/.*$", "", i),"_", sample.map[,"ID"],"_", sample.map[,"Well"]),
                                each = 143259)
  
  dat[[i]] <- dat[[i]][,c(1,2,10,11)]
}

# rbind the results together, saving selected columns
dat <- do.call("rbind", dat)

# create matrices that are snps x samples
snps <- unique(dat[,"SNP Name"])
samples <- unique(dat[,"Sample ID"])
X <- Y <- matrix(ncol=length(samples), nrow=length(snps))
dimnames(X) <- dimnames(Y) <- list(snps, samples)
for(i in seq(along=samples)) {
  message(i, " of ", length(samples))
  tmp <- dat[dat[,"Sample ID"]==samples[i],]
  X[,samples[i]] <- tmp[,"X"]
  Y[,samples[i]] <- tmp[,"Y"]
}

# bring together in one matrix
result <- cbind(snp=rep(snps, 2),
                channel=rep(c("X", "Y"), each=length(snps)),
                as.data.frame(rbind(X, Y)))
rownames(result) <- 1:nrow(result)

# bring SNP rows together
result <- result[as.numeric(t(cbind(seq_along(snps), seq_along(snps)+length(snps)))),]
rownames(result) <- 1:nrow(result)

save(result, file = "/projects/heh/csna_workflow/data/Jackson_Lab_Bubier_MURGIGV01/intensities.fst.RData")
