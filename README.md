PRISMA - Automated Polygenic Risk Score modeling and evaluation
================================================================================

prisma is an R package for the semi-automated construction and evaluation of 
Polygenic Risk Scores (PRS) with data from multiple summary statistics and PRS 
extraction algorithms. prisma supports several Genome Wide Association (GWA) 
tests for the calculation of summary statistics and integrates the results of 
multiple PRS algorithms. In addition, it offers multiple PRS evaluation criteria 
and a semi-automated optimal PRS selection process through the local 
maximization of evaluation criteria. Finally, it integrates GWA effects with 
linear programming and helps the user select the best PRS through a rich report 
of the results and online search in the GWAS Catalog.

# Installation

prisma requires the presence of several R packages but also additional 3rd party
software tools for running GWA tests to derive summary statistics as well as PRS
extraction. See the more detailed 
[installation](#) 
guide for more information and step-by-step process.

## R packages

### Basic installation of Bioconductor

```
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.15")
```

### Additional packages

The following packages are either on CRAN or on Bioconductor:

```
pkgs <- c("akmedoids","Biobase","biomaRt","BSgenome","data.table",
    "GenomicRanges","ggplot2","gwasrapidd","harmonicmeanp","islasso",
    "jsonlite","kableExtra","knitr","liftOver","Matrix","magrittr","methods",
    "openxlsx","parallel","PhenotypeSimulator","quincunx","R.utils",
    "rmarkdown","rmdformats","RcppArmadillo","RMTstat","rrBLUP","rrcov","rsnps",
    "rtracklayer","S4Vectors","scales","scrime","snpStats","SNPRelate",
    "statgenGWAS","SummarizedExperiment","survcomp","testthat","tseries",
    "utils")    

BiocManager::install(pkgs)
```

Also, `lassosum` is not on CRAN:

```
library(devtools)
install_github("tshmak/lassosum")
```

## Other external tools and resources

The following tools should be findable by R, for example they should be included
in the system's `PATH` environmental variable.

### Summary statistics

* SNPTEST
* PLINK

### PRS extraction

* PRSice-2

### Other tools

* GTOOL
* QCTOOL
* IMPUTE2
