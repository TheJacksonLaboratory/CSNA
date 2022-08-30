library(qtl2)
library(tidyverse)
library(SimDesign)
library(MASS)
library(EnvStats)

#rankz transform
rz.transform <- function(y) {
  rankY=rank(y, ties.method="average", na.last="keep")
  rzT=qnorm(rankY/(length(na.exclude(rankY))+1))
  return(rzT)
}

###############
# heritability
#####
# function to simulate data to estimate heritability intervals
h2_sim <- function(trait, K, covar, nsim=100){
  # trait, K and covar are matrices with row.names in qtl2-style
  # covar adjustment applied to data but not in simulated data

  # reduce mice to common set with kinship and trait values
  # mice.names <- intersect(row.names(K), row.names(trait))
  # K <- K[mice.names, mice.names]
  # trait <- as.matrix(trait[mice.names,],ncol=1)
  N <- dim(trait)[1]

  # estimate heritability
  h2 <- as.numeric(est_herit(trait, K, covar))

  # simulate from kinship matrix and re-estimate
  h2_sim <- rep(0,nsim)
  for(i in 1:nsim){
    trait_sim <- t(rmvnorm(n = 1, mean = rep(0,N), sigma = h2*2*K + (1-h2)*diag(1,N)))
    rownames(trait_sim) <- rownames(trait)
    h2_sim[i] <- as.numeric(est_herit(trait_sim, K))
  }

  # return estimated h2 and simulations
  list(N=N, h2=h2, sim=h2_sim)
}

setwd("/projects/csna/csna_workflow/")
#Prj01_RL_pheno_qttl2_DO_11092020
pheno <- readr::read_csv("data/pheno/Prj01_RL_pheno_qttl2_DO_11092020.csv", na = c("N/A", "NA"))
pheno <- pheno %>%
  filter(apply(pheno, 1, function(x) sum(is.na(x))) < 5) %>% # remove rows with 5 NAs
  mutate(
    Acq.Anticipatory.Correct.Responses.nooutlier = case_when(
      (Acq.Anticipatory.Correct.Responses > 5 | is.na(Acq.Anticipatory.Correct.Responses)) ~ NA_real_,
      TRUE ~ Acq.Anticipatory.Correct.Responses),
    Rev.Anticipatory.Incorrect.Responses.nooutlier = case_when(
      (Rev.Anticipatory.Incorrect.Responses > 6) ~ NA_real_,
      TRUE ~ Rev.Anticipatory.Incorrect.Responses
    )) %>% #boxcoxtransformation
  mutate(across(c(Rev.Anticipatory.Incorrect.Responses), ~ EnvStats::boxcoxTransform(.x, lambda = EnvStats::boxcox(.x, optimize = TRUE)$lambda),
                .names ="{.col}.boxcox.trans")) %>% #omit values >6 for Acq.Anticipatory.Correct.Responses and Rev.Anticipatory.Incorrect.Responses variables.
  mutate(across(Acq.Total.Trials:Rev.Anticipatory.Incorrect.Responses.nooutlier, ~ rz.transform(.x),
                .names ="{.col}.rz"))

#load gm data
load("data/Jackson_Lab_12_batches/gm_DO3173_qc_newid.RData")
#subset
overlap.id <- intersect(pheno$mouseID, ind_ids(gm_after_qc))
gm_after_qc <- gm_after_qc[overlap.id, ]
saveRDS(gm_after_qc, file = "output/RL_prj/RL_prj_gm_after_qc.RDS")

load("data/Jackson_Lab_12_batches/qc_info.RData")
#merge pheno with covar
pheno <- pheno %>%
  left_join(qc_info) %>%
  filter(mouseID %in% overlap.id) %>%
  filter(bad.sample == FALSE) %>%
  mutate(sex = case_when(
    sex == "M" ~ 1,
    sex == "F" ~ 0
  )) %>%
  mutate_at(vars(sex,
                 ngen),
            list(factor))

#pr
load("data/Jackson_Lab_12_batches/apr_69kchr_DO3173.RData") #apr.69kchr
apr = apr.69kchr[overlap.id,]
str(apr)
saveRDS(apr, file = "output/RL_prj/apr.RDS")

#kinship
k = calc_kinship(probs = apr, type = "loco", use_allele_probs = TRUE, cores = 20)
k_overall = calc_kinship(probs = apr, type = "overall", use_allele_probs = TRUE, cores = 20)

#pgmap
load("/projects/csna/csna_workflow/data/69k_grid_pgmap.RData")

#pheno list
pheno.list <- colnames(pheno)[4:18]
pheno <- as.data.frame(pheno)

# qtlmapping  -------------------------------------------------------------
query_variants <- create_variant_query_func("/projects/csna/csna_workflow/data/cc_variants.sqlite")
query_genes <- create_gene_query_func("/projects/csna/csna_workflow/data/mouse_genes_mgi.sqlite")

# #heritability -----------------------------------------------------------
# loop through traits
pheno.names <- pheno.list
pheno.herit <- data.frame(Name = pheno.names, h2=0, N=0, ll=0, ul=0)
covar = model.matrix(~sex, data = pheno)[,-1, drop=F]
rownames(covar) <- pheno$mouseID
colnames(covar)[1] = "sex"
#loop
for(i in 1:length(pheno.names)){
  j = pheno.names[[i]]
  print(j)
  phe = pheno[,j,drop=F]
  rownames(phe) = pheno$mouseID
  tmp <- h2_sim(phe, k_overall, covar, 100)
  pheno.herit$N[i] <- tmp$N
  pheno.herit$h2[i] <- tmp$h2
  pheno.herit$ll[i] <- quantile(tmp$sim, 0.1)
  pheno.herit$ul[i] <- quantile(tmp$sim, 0.9)
}
write.csv(pheno.herit, file = "output/Prj01_RL_11092020_cutoff6.pheno.herit.csv",
          row.names = FALSE, quote = F)


# m2 model only Generation
model = "m2"
if(model == "m2"){ # ~Generation
  #scan
  m2.qtl.out <- list()
  m2.permu   <- list()
  m2.sigqtl.peak <- list()
  m2.sigqtl.snps <- list()
  m2.sigqtl.genes <- list()
  m2.sigqtl.blup <- list()
  m2.sigqtl.coef <- list()
  for(i in pheno.list){
    print(i)
    phe.i = pheno[,i]
    names(phe.i) <- pheno$mouseID
    #covarite
    addcovar.i = model.matrix(~sex+ngen, data = pheno)[,-1]
    rownames(addcovar.i) <- pheno$mouseID
    colnames(addcovar.i)[1] = "sex"

    #check the order
    print(all.equal(rownames(k$`1`), rownames(apr$`1`)))
    print( all.equal(names(phe.i), rownames(apr$`1`)))
    print( all.equal(names(phe.i), rownames(addcovar.i)))

    #mapping
    m2.qtl.out[[i]] <- scan1(genoprobs    =  apr,
                             pheno        =  phe.i,
                             kinship      =  k,
                             addcovar     =  addcovar.i,
                             cores        =  20) # changed to cluster list
    #permu
    print("starting permutation")
    print(Sys.time())
    m2.permu[[i]] <- scan1perm(genoprobs        =  apr,
                               pheno            =  phe.i,
                               kinship          =  k,
                               addcovar         =  addcovar.i,
                               n_perm           =  1000,
                               cores            =  20) # changed to cluster list
    # print("done permutation")
    print(Sys.time())
    cutoff <- 6#cutoff_line[[i]]
    print(cutoff)
    #peaks significant at p 0.1
    m2.sigqtl.peak[[i]] <- find_peaks(m2.qtl.out[[i]], pmap,
                                      threshold=cutoff,
                                      drop = 1.5)
    print(m2.sigqtl.peak[[i]])
    m2.sigqtl.peak[[i]]$diff_ci <- m2.sigqtl.peak[[i]]$ci_hi - m2.sigqtl.peak[[i]]$ci_lo
    #for inveral greater than 50
    indx <- which(m2.sigqtl.peak[[i]]$diff_ci >= 50)
    if(length(indx)>0){
      for(j in 1:length(indx)){
        lod_intval <- lod_int(m2.qtl.out[[i]],
                              pmap,
                              chr = m2.sigqtl.peak[[i]][indx[j],"chr"],
                              peakdrop = 5)
        lod_intval <- lod_intval[lod_intval[,2]==m2.sigqtl.peak[[i]][indx[j],"pos"],]
        print(lod_intval)
        m2.sigqtl.peak[[i]][indx[j],"ci_lo"] <- lod_intval[1]
        m2.sigqtl.peak[[i]][indx[j],"ci_hi"] <- lod_intval[3]

      }
      m2.sigqtl.peak[[i]]$diff_ci_fixed <-  m2.sigqtl.peak[[i]]$ci_hi - m2.sigqtl.peak[[i]]$ci_lo
      print(m2.sigqtl.peak[[i]])
    }

    #blup, gene and snp
    if(dim(m2.sigqtl.peak[[i]])[1] != 0){
      m2.sigqtl.snps[[i]] <- list()
      m2.sigqtl.genes[[i]] <- list()
      m2.sigqtl.blup[[i]] <- list()
      m2.sigqtl.coef[[i]] <- list()
      for(j in 1:dim(m2.sigqtl.peak[[i]])[1]){
        sigqtl.chr <- as.character(m2.sigqtl.peak[[i]][j,"chr"])
        print(sigqtl.chr)
        m2.sigqtl.blup[[i]][[j]] <- scan1blup(genoprobs    =  apr[,sigqtl.chr],
                                              pheno        =  phe.i,
                                              kinship      =  k[[sigqtl.chr]],
                                              addcovar     =  addcovar.i,
                                              cores        =  20) # changed to cluster list

        #coeff
        m2.sigqtl.coef[[i]][[j]] <- scan1coef(genoprobs    =  apr[,sigqtl.chr],
                                              pheno        =  phe.i,
                                              kinship      =  k[[sigqtl.chr]],
                                              addcovar     =  addcovar.i)

        m2.sigqtl.snps[[i]][[j]] <- scan1snps(genoprobs = apr,
                                              map = pmap,
                                              pheno = phe.i,
                                              kinship = k[[sigqtl.chr]],
                                              addcovar = addcovar.i,
                                              query_func=query_variants,
                                              chr=sigqtl.chr,
                                              start=m2.sigqtl.peak[[i]][j,"ci_lo"],
                                              end=m2.sigqtl.peak[[i]][j,"ci_hi"],
                                              keep_all_snps=TRUE)
        m2.sigqtl.genes[[i]][[j]] <- query_genes(sigqtl.chr,
                                                 m2.sigqtl.peak[[i]][j,"ci_lo"],
                                                 m2.sigqtl.peak[[i]][j,"ci_hi"])
      }
    }else{
      m2.sigqtl.snps[[i]] <- list()
      m2.sigqtl.genes[[i]] <- list()
      m2.sigqtl.blup[[i]] <- list()
      m2.sigqtl.coef[[i]] <- list()
    }
  }
  save(m2.qtl.out,
       m2.permu,
       m2.sigqtl.peak,
       m2.sigqtl.snps,
       m2.sigqtl.genes,
       m2.sigqtl.blup,
       m2.sigqtl.coef,
       file = paste0("output/","Prj01_RL_11092020_cutoff6.qtlout.69k.RData"))
}
