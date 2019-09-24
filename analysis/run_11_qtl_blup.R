####################################################################################################################
#This script submits jobs to perform bluptation for qtl mapping in DO mice using 11_qtl_blup.R

phenocsv <- c("Novelty_resids_datarelease_07302918.csv",
              "Novelty_residuals_RankNormal_datarelease_07302918.csv",
              "Novelty_no_outliers_datarelease_07302918.csv",
              "Novelty_raw_datarelease_07302918.csv")

model <- c("m1",
  	   "m2")
           #"m3")

#mkdir folder for bluptation
out.dir <- "/projects/heh/csna_workflow/output/blup/"
if(!dir.exists(out.dir)){
  system(paste0("mkdir ", out.dir))
}

for(p in phenocsv[1:2]){
  for(m in model){
    for(i in 14:41){ #index col number for the phenotype in the csv
      pwd <- "/projects/heh/csna_workflow/analysis/"
      job.name <- paste0("blup_",i,"_",m,"_", gsub(".csv","",p))
      R.file      <- paste0(pwd, "11_qtl_blup.R")
      script.file <- paste0(pwd, job.name, ".script")
      Rout.file   <- paste0(pwd, job.name, ".Rout")
      stderr      <- paste0(pwd, job.name, ".stderr")
      stdout      <- paste0(pwd, job.name, ".stdout")
      
      ##Write script
      write("#PBS -S /bin/bash", script.file, append = F)
      write(paste0("#PBS -e ",stderr), script.file, append = T)
      write(paste0("#PBS -o ",stdout), script.file, append = T)
      write("#PBS -l nodes=1:ppn=20,walltime=02:00:00:00", script.file, append =T)
      #write("#PBS -l mem=256gb", script.file, append =T)
      write(paste0("#PBS -N ",job.name), script.file, append = T)
      
      write("module load R/3.3.2", script.file, append = T)
      write(paste0("R CMD BATCH --no-restore --no-save --no-readline",
                   " -phenocsv=", paste0("/projects/heh/csna_workflow/data/pheno/",p),
                   " -pheno.idx=", i,
                   " -model=", m,
                   " -outdir=", out.dir,
                   " ", R.file,
                   " ", Rout.file),# Rout file
            script.file, append = T)
      
      write("date", script.file, append = T)
      write("sleep 3",script.file,append=T)
      
      njobs.running <- as.numeric(system("qstat -u heh | wc -l", intern=T, wait=T))-5
      while(njobs.running>=200){
        system("sleep 100", wait=T)
        njobs.running <- as.numeric(system("qstat -u heh | wc -l", intern=T, wait=T))-5
      }
      
      system(paste0("qsub -V ",script.file),wait=T)
      system("sleep 3", wait=T)
    }
  }
}
