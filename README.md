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

## System packages

Some or all of these may be already installed on the system. There might be
additional ones depending on the Linux distribution prisma is being installed
to. This list will be maintained and updated.

```
sudo apt install -y apt-transport-https software-properties-common \
    build-essential zlib1g-dev libdb-dev libcurl4-openssl-dev libssl-dev \
    libxml2-dev apache2 libfontconfig1-dev libjpeg-dev
```

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
pkgs <- c("BSgenome.Hsapiens.UCSC.hg19","BSgenome.Hsapiens.UCSC.hg38",
    "Biobase","biomaRt","BSgenome","data.table","GenomicRanges","ggplot2",
    "gwasrapidd","harmonicmeanp","islasso","jsonlite","kableExtra","knitr",
    "liftOver","Matrix","magrittr","methods","openxlsx","pander","parallel",
    "PhenotypeSimulator","quincunx","R.utils","rmarkdown","rmdformats",
    "RcppArmadillo","RMTstat","rrBLUP","rrcov","rsnps","rtracklayer",
    "S4Vectors","scales","scrime","snpStats","SNPRelate",
    "SNPlocs.Hsapiens.dbSNP151.GRCh38","statgenGWAS","SummarizedExperiment",
    "survcomp","testthat","tseries","utils")

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

* [SNPTEST](https://www.well.ox.ac.uk/~gav/snptest/)
* [PLINK](https://www.cog-genomics.org/plink/)

### PRS extraction

* [PRSice-2](https://www.prsice.info/)

### Other tools

* [GTOOL](https://www.well.ox.ac.uk/~cfreeman/software/gwas/gtool.html)
* [QCTOOL](https://www.well.ox.ac.uk/~gav/qctool_v2/)
* [IMPUTE2](http://mathgen.stats.ox.ac.uk/impute/impute_v2.html)

## Installation of PRISMA

PRISMA can be installed from GitHub for the time being:

```
library(devtools)
install_github("pmoulos/prisma")
```

# Quickstart

The following depict the very basic PRISMA PRS extraction pipeline. For full
documentation see the package vignettes and also 
[online](#) documentation.

## Step by step pipeline

```
library(prisma)

data(toy,package="prisma")

# Lose filters for the toy dataset
filts <- getDefaults("filters")
filts$IBD <- NA
filts$hwe <- 1e-3
filts$pcaOut <- FALSE
filts$inbreed <- NA

# Workspace
d1 <- tempfile()
wspace1 <- file.path(tempdir(),d1)

# Traits
response <- "BMI"
covariates <- c("Case_Control","Gender","Age")

# Complete 4-iteration PRS pipeline - ~1 min to run with 4 cores - no PCA
prismaOut <- prisma(
    gwe=toy,
    phenotype=response,
    covariates=covariates,
    pcs=FALSE,
    trainSize=0.8,
    niter=4,
    resolution="frequency",
    minSnps=2,
    filters=filts,
    pcaMethod="snprel",
    imputeMissing=FALSE,
    gwaMethods=c("glm","statgen"),
    prsiceOpts=list(clump_p=0.75,perm=100),
    prsWorkspace=wspace1,
    cleanup="intermediate",
    logging="file",
    output="normal",
    rc=0.25
)

# Below, individual wrapper steps implemented in prisma as whole
# Output only of the discovery PRS pipeline with statgenGWAS
d2 <- tempfile()
wspace2 <- file.path(tempdir(),d2)
dnList <- prsPipeline(
    gwe=toy,
    phenotype=response,
    covariates=covariates,
    pcs=FALSE,
    trainSize=0.8,
    niter=4,
    filters=filts,
    pcaMethod="snprel",
    imputeMissing=FALSE,
    gwaMethods="statgen",
    prsiceOpts=list(clump_p=0.75,perm=100),
    prsWorkspace=wspace2,
    cleanup="intermediate",
    logging="file",
    rc=0.25
)

# After, calculate the evaluation metrics
d3 <- tempfile()
wspace3 <- file.path(tempdir(),d3)
evalMetrics <- prsSelection(
    dnList=dnList,
    gwe=toy,
    phenotype=response,
    covariates=covariates,
    pcs=FALSE,
    trainSize=0.8,
    niter=4,
    resolution="frequency",
    minSnps=2,
    filters=filts,
    pcaMethod="snprel",
    imputeMissing=FALSE,
    gwaMethods="statgen",
    prsiceOpts=list(clump_p=0.75,perm=100),
    prsWorkspace=wspace3,
    cleanup="intermediate",
    logging="file",
    output="summary",
    useDenovoWorkspace=dnList,
    rc=0.25
)

candidates <- selectPrs(
    metrics=evalMetrics$metrics,
    snpSelection=evalMetrics$pgs,
    gwe=toy,
    base=evalMetrics$baseline
)
```

## Complete pipeline

```
if (Sys.which("PRSice_linux") != "") {
    outPath <- tempfile()
    wspace1 <- tempfile()
    
    prismaPipeline(
        gwe=toy,
        phenotype=response,
        covariates=covariates,
        pcs=FALSE,
        trainSize=0.8,
        niter=4,
        resolution="frequency",
        minSnps=2,
        filters=filts,
        pcaMethod="snprel",
        imputeMissing=FALSE,
        gwaMethods=c("glm","statgen"),
        prsiceOpts=list(clump_p=0.75,perm=100),
        prsWorkspace=wspace4,
        cleanup="intermediate",
        logging="file",
        output="normal",
        rc=0.25,
        outPath=outPath
    )
}
```
