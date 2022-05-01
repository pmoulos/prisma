\name{importGWAS}
\alias{importGWAS}
\title{Create GWASExperiments from PLINK files}
\usage{
    importGWAS(input, phenos = NULL, backend = c("snpStats",
        "bigsnpr"), selection = NULL, genome = NA_character_,
        alleleOrder = c("plink", "reverse"), writeGds = TRUE,
        gdsfile = ifelse(backend == "snpStats", tempfile(), NA),
        gdsOverwrite = TRUE)
}
\arguments{
    \item{input}{A list with the full path to a PLINK file 
    triplet (BED/BIM/FAM) or a string with the triplet 
    basename which will be suffixed with bed/bim/fam. See
    als Details.}
    
    \item{phenos}{A data frame or a file with the dataset
    phenotypes. It must have the same number of rows as the
    samples present in the PLINK triplet.}
    
    \item{backend}{the R backend to use for object management.
    See Details.}
    
    \item{selection}{A list with members named \code{snps}
    and \code{samples} to be used as the related 
    \code{select.snps} and \code{select.subjects} arguments
    of the \code{\link{read.plink}} function from snpStats.}
    
    \item{genome}{The reference genome which the SNPs are
    mapped to. Typically one of \code{"hg19"} or 
    \code{"hg38"} if this information is available, otherwise
    \code{NA}.}
    
    \item{alleleOrder}{Should the alleles be reversed? If yes
    choose \code{"reverse"}.}
    
    \item{writeGds}{Write a GDS supporting file while 
    importing PLINK files? In some cases such as reading
    files for external imputation is not needed. Defaults 
    to \code{TRUE}.}
    
    \item{gdsfile}{A GDS file supported by the package
    SNPRelate is necessary for some operations like
    SNPRelate PCA and LD. Where should this be written?
    It is adivsed to explicitly provide this path.}
    
    \item{gdsOverwrite}{overwrite the provided GDS file if
    exists? Defaults to \code{TRUE}.}
}
\value{
    A GWASExperiment object with the contents of the PLINK
    files triplet and optionally, additional associated
    phenotypes.
}
\description{
    This function creates a \code{GWASExperiment} object
    from PLINK BED/BIM/FAM files using functions from the
    package snpStats. In the future, importing with the
    package bigSNP will be supported.
}
\details{
    Regarding the \code{input} argument, when this is a list,
    the latter must have named members with names \code{bed},
    \code{bim}, \code{fam}. When this is a single string, it
    must point to a full path where the triplet is found and
    must not end with any of \code{bim}, \code{bed}, 
    \code{fam} as these are appended to the base.
    
    Regarding the \code{backend} argument, current PRISMA
    design goals are to support both importing and processing
    data with the R package snpStats and the R package 
    bigSNP. Currently, only snpStats is supported.
}
\examples{
    fam <- system.file("extdata/sample.fam",package="snpStats")
    bim <- system.file("extdata/sample.bim",package="snpStats")
    bed <- system.file("extdata/sample.bed",package="snpStats")
    input <- list(fam=fam,bim=bim,bed=bed)
    sample <- snpStats::read.plink(input$bed,input$bim,input$fam)
    
    # Without phenotypes
    gwe <- importGWAS(input,backend="snpStats")
    
    # With phenotypes
    set.seed(42)
    pseudopheno <- data.frame(
        case_control=sample(c(0,1),nrow(sample$fam),replace=TRUE),
        other_pheno=sample(c("drug","nodrug"),nrow(sample$fam),replace=TRUE),
        cont=round(runif(nrow(sample$fam)),3),
        row.names=rownames(sample$fam)
    )
    gwe <- importGWAS(input,phenos=pseudopheno,backend="snpStats")
}
\author{
    Panagiotis Moulos
}

