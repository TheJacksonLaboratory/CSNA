## Diversity Outcross data from Recla et al (2014) and Logan et al. (2013)

### Source

DO Phenotypes from this study <https://phenome.jax.org/projects/Chesler4>

Founder genotypes from <ftp://ftp.jax.org/MUGA/>


### Files

- [`DO2014.json`](DO2014.json) - Control file in JSON format
- [`DO2014_covar.csv`](DO2014_covar.csv) - covariate data (individuals x
  covariates)
- [`DO2014_pheno.csv`](DO2014_pheno.csv) - phenotype data (individuals x
  phenotypes)
- [`DO2014_geno.csv`](DO2014_geno.csv) - genotype data (markers x individuals)
- [`do_foundergeno.csv`](DO2014_foundergeno.csv) - founder genotype data
  (markers x founders)
- [`DO2014_gmap.csv`](DO2014_gmap.csv) - Genetic map of markers (positions in
  cM)
- [`DO2014_pmap.csv`](DO2014_pmap.csv) - Physical map of markers (positions in
  NCBI38/mm10 Mbp)

The data are also available as a zip file, [`DO2014.zip`](DO2014.zip).

### File format

See the [R/qtl2 input file format](https://kbroman.org/qtl2/assets/vignettes/input_files.html).


### Citations

Recla JM, Robledo RF, Gatti DM, Bult CJ, Churchill GA, Chesler EJ (2014)
[Precise genetic mapping and integrative bioinformatics in Diversity Outbred mice reveals Hydin as a novel pain gene](https://www.ncbi.nlm.nih.gov/pubmed/24700285).
Mamm Genome 25:211-222

Logan RW, Robledo RF, Recla JM, Philip VM, Bubier JA, Jay JJ, Harwood
C, Wilcox T, Gatti DM, Bult CJ, Churchill GA, Chesler EJ (2013)
[High-precision genetic mapping of behavioral traits in the diversity outbred mouse population](https://www.ncbi.nlm.nih.gov/pubmed/23433259).
Genes Brain Behav 12:424-437


### Use with [R/qtl2](https://kbroman.org/qtl2)

Load these data into R directly from the web as follows:

```r
library(qtl2)
file <- paste0("https://raw.githubusercontent.com/rqtl/",
               "qtl2data/master/DO_Recla/DO2014.zip")
DO2014 <- read_cross2(file)
```
