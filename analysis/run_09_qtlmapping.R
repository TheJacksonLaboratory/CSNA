####################################################################################################################
#This script submits jobs to perform qtl mapping for DO mice using 08_qtlmapping.R

phenocsv <- c("Novelty_resids_datarelease_07302918.csv",
              "Novelty_residuals_RankNormal_datarelease_07302918.csv",
              "Novelty_no_outliers_datarelease_07302918.csv",
              "Novelty_raw_datarelease_07302918.csv")

model <- c(#"m1",
	   #"m2",
           "m3")

for(p in phenocsv[1:2]){
  for(m in model){
    pwd <- "/projects/heh/csna_workflow/analysis/"
    job.name    <- paste0(gsub(".csv","",p),"_", m)
    R.file      <- paste0(pwd, "09_qtlmapping.R")
    script.file <- paste0(pwd, job.name, ".script")
    Rout.file   <- paste0(pwd, job.name, ".Rout")
    stderr      <- paste0(pwd, job.name, ".stderr")
    stdout      <- paste0(pwd, job.name, ".stdout")
    
    ##Write script
    write("#PBS -S /bin/bash", script.file, append = F)
    write(paste0("#PBS -e ",stderr), script.file, append = T)
    write(paste0("#PBS -o ",stdout), script.file, append = T)
    write("#PBS -l nodes=1:ppn=20,walltime=02:00:00:00", script.file, append =T)
    write("#PBS -l mem=256gb", script.file, append =T)
    write(paste0("#PBS -N ",job.name), script.file, append = T)
    
    write("module load R/3.3.2", script.file, append = T)
    write(paste0("R CMD BATCH --no-restore --no-save --no-readline",
                 " -phenocsv=", paste0("/projects/heh/csna_workflow/data/pheno/",p),
                 " -firstpheno.idx=", 14,
                 " -lastpheno.idx=", 41,
                 " -model=", m,
                 " -outdir=", "/projects/heh/csna_workflow/output/",
                 " ", R.file,
                 " ", Rout.file),# Rout file
          script.file, append = T)
    
    write("date", script.file, append = T)
    write("sleep 3",script.file,append=T)
    
    njobs.running <- as.numeric(system("qstat -u heh | wc -l", intern=T, wait=T))-5
    while(njobs.running>=100){
      system("sleep 100", wait=T)
      njobs.running <- as.numeric(system("qstat -u heh | wc -l", intern=T, wait=T))-5
    }
    
    system(paste0("qsub -V ",script.file),wait=T)
    system("sleep 3", wait=T)
  }
}
