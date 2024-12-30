---
output:
  pdf_document: default
  html_document: default
---

# Installation of required R packages

The following R packages are required for the successful completion of the following workshop modules:
1. Gene Expression 
2. QTL Mapping 

## Gene Expression
In the gene expression module we will be undertaking a differential gene expression analysis. For this, we will make use of the [DESeq2](https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html) R package. To install this package, copy and paste the commands below in your R console:

`source("https://bioconductor.org/biocLite.R")`

`biocLite("DESeq2")`

In addition to DESeq2 we will need the following R packages as well:

```
source("https://bioconductor.org/biocLite.R")
#libraries for gene expression analysis
biocLite("DESeq2")
biocLite("vsn")

#libraries for table manipulations
install.packages("DT")
install.packages("plyr")

#libraries for visualization
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("RColorBrewer")

#libraries for gene annotation and enrichement analysis
biocLite("org.Mm.eg.db")
biocLite("topGO")
```

## QTL Mapping

QTL mapping workshop will require the installation of the following  R libraries. Copy and paste the commands below in you R console:

```
install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
install.packages("GGally")
```


# External datasets that need to be downloaded 

The QTL mapping workshop, particularly the one on Diversity Outbred mice, has a section on **SNP Association Mapping** that requires the following two files:

- cc_variants.sqlite [Download here](https://doi.org/10.6084/m9.figshare.5280229.v2) : These are the variants in the Collaborative Cross founders (3 GB)
- mouse_genes.sqlite [Download here](https://doi.org/10.6084/m9.figshare.5280238.v4) : full set of mouse gene annotations (677 MB)

