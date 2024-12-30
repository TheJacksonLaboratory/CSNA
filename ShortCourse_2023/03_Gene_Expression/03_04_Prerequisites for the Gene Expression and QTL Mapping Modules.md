---
output:
  pdf_document: default
  html_document: default
---

# Installation of required R packages
The following R packages are required for the successful completion of the following workshop modules: 

1. Gene Expression 

2. QTL Mapping

Most of the packages will be installed using BiocManager. Our first step therefore, is to install BiocManager (<https://bioconductor.org/install>),

## Install BiocManager

```         
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install()
```

## Gene Expression

In the gene expression module we will be undertaking a differential gene expression analysis. For this, we will make use of the [DESeq2](https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html) R package. To install this package, copy and paste the commands below in your R console:

```         
BiocManager::install("DESeq2")
```

In addition to DESeq2 we will need the following R packages as well:

### libraries for gene expression analysis

```         
BiocManager::install("vsn", force = TRUE)
```

### libraries for table manipulations

```         
BiocManager::install("DT", force = TRUE)
BiocManager::install("plyr", force = TRUE)
```

### libraries for visualization

```         
BiocManager::install("ggplot2", force = TRUE)
BiocManager::install("pheatmap", force = TRUE)
BiocManager::install("RColorBrewer", force = TRUE)
```

###libraries for gene annotation and enrichement analysis

```         
BiocManager::install("org.Mm.eg.db", force = TRUE)
BiocManager::install("topGO", force = TRUE)
```

## QTL Mapping

QTL mapping workshop will require the installation of the following R libraries. Copy and paste the commands below in you R console:

```         
BiocManager::install("qtl2", force = TRUE)
BiocManager::install("GGally", force = TRUE)
```

# External datasets that need to be downloaded

The QTL mapping workshop, particularly the one on Diversity Outbred mice, has a section on **SNP Association Mapping** that requires the following two files:

-   cc_variants.sqlite [Download here](https://doi.org/10.6084/m9.figshare.5280229.v2) : These are the variants in the Collaborative Cross founders (3 GB)
-   mouse_genes.sqlite [Download here](https://doi.org/10.6084/m9.figshare.5280238.v4) : full set of mouse gene annotations (677 MB)
