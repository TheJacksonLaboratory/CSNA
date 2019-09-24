####################################################################################################################
#   This script will run the 07_recomb_size_founder_prop.Rmd
#
#   Author: Hao He
#   Date:   02/13/2019
#   E-mails: hao.he@jax.org
####################################################################################################################

library(workflowr)
library(rmarkdown)
setwd(dir="/projects/heh/csna_workflow/")
Sys.setenv(RSTUDIO_PANDOC="/usr/lib/rstudio-server/bin/pandoc")
workflowr::wflow_status()
# publish a single file
file = c("analysis/index.Rmd",
         "analysis/about.Rmd",
         "analysis/license.Rmd",
	 "analysis/07_recomb_size_founder_prop.Rmd")
msg =  sub(".*/", "", file[4])
wflow_publish(files = file,
              message =msg)
