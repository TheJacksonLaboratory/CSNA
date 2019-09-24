####################################################################################################################
#   This script performs GCTA heritability for DO mice using 08_gcta_herit.R
#   It takes  args in this script
#   -phenocsv= location of phenotype .csv file
#   -firstpheno.idx= column number of first phenotype
#   -lastpheno.idx= column number of last phenotype
#
#   Author: Hao He
#   Date:   02/13/2019
#   E-mails: hao.he@jax.org
####################################################################################################################
library(ggplot2)

# args --------------------------------------------------------------------
options(warn=1)
args <- commandArgs()
args
options(stringsAsFactors = FALSE)

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

#mkdir folder for this file
out.dir <- paste0("/projects/heh/csna_workflow/data/GCTA/", name.phenocsv)
if(!dir.exists(out.dir)){
  system(paste0("mkdir ", out.dir))
}
setwd(out.dir)

#pheno id
pheno <- read.csv(phenocsv, header = TRUE)

#gcta id
gcta.id <- read.table("/projects/heh/csna_workflow/data/GCTA/nine_batches.fam", header = F, sep = " ")

#overlap
pheno.overlap <- pheno[pheno$Mouse.ID %in% intersect(gcta.id$V1, pheno$Mouse.ID),]
pheno.overlap$Sex[pheno.overlap$Sex == "Male"] <- 1
pheno.overlap$Sex[pheno.overlap$Sex == "Female"] <- 2

#subset id
idlist <- cbind(pheno.overlap$Mouse.ID, pheno.overlap$Mouse.ID)
write.table(idlist, file = "update_id.txt", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

#update sex
sex <- data.frame(id1 = pheno.overlap$Mouse.ID,
                  id2 = pheno.overlap$Mouse.ID,
                  sex = as.character(pheno.overlap$Sex))
write.table(sex, file = "update_sex.txt",sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

#update pheno
update.pheno <- cbind(data.frame(FID = pheno.overlap$Mouse.ID,
                                 IID = pheno.overlap$Mouse.ID),
                      pheno.overlap[,14:41]
)
write.table(update.pheno, file = "update.pheno.txt",sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


#subset id and update sex
system(paste0("/projects/heh/csna_workflow/data/GCTA/plink -bfile ", 
       "/projects/heh/csna_workflow/data/GCTA/","nine_batches ",
       "--keep update_id.txt ",
       "--update-sex update_sex.txt ",
       "--make-bed --out ",
       name.phenocsv))

# Estimate the GRM
system(paste0("/projects/heh/csna_workflow/data/GCTA/gcta64  --bfile ",
              name.phenocsv,  " --make-grm  --out ",
              name.phenocsv))

for(i in 1:length(colnames(pheno.overlap)[pheno.first:pheno.last])){
  system(
    paste0("/projects/heh/csna_workflow/data/GCTA/gcta64 --reml  --grm ",
           name.phenocsv,  " --pheno update.pheno.txt --covar update_sex.txt --mpheno ",
           i," --out ",colnames(pheno.overlap)[14:41][i],".out")
  )
}

#plot heritability by GCTA
hsq.gcta <- list()
phe.name <- colnames(pheno.overlap)[pheno.first:pheno.last]
for (i in phe.name){
  hsq.gcta[[i]] <- read.table(file = paste0(i,".out.hsq"),sep = "\t", header = FALSE,fill = TRUE, stringsAsFactors = FALSE)
}
h <- data.frame(
  Phenotype = phe.name,
  Heritability = as.numeric(as.vector(unlist(lapply(hsq.gcta, FUN = function(x){x[5,2]})))),
  Sample_size = as.numeric(as.vector(unlist(lapply(hsq.gcta, FUN = function(x){x[11,2]})))),
  Domain = sub("\\..*", "", phe.name),
  stringsAsFactors = FALSE)

#histgram
h$Heritability <- round(h$Heritability,2)
pdf(file = paste0(name.phenocsv,"_heritability_by_GCTA.pdf"), height = 10, width = 10)
p<-ggplot(data=h, aes(x=Phenotype, y=Heritability, fill=Domain, color = Domain)) +
  geom_bar(stat="identity") +
  scale_y_continuous(breaks=seq(0.0, 1.0, 0.1)) +
  geom_text(aes(label = Heritability, y = Heritability + 0.005), position = position_dodge(0.9),vjust = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
dev.off()
