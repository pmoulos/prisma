\name{writeFuns}
\alias{writeGdsFile}
\alias{writePlink}
\title{PLINK and GDS files from GWASExperiment}
\usage{
    writeGdsFile(obj, gdsfile = NULL)
    
    writePlink(obj, pheno = NULL, outBase = NULL, salvage = FALSE,
        reverse = FALSE, perChr = FALSE, overwrite = TRUE)
}
\arguments{
    \item{obj}{An object of GWASExperiment class with at
    least genotypic data.}
    
    \item{gdsfile}{The path to the GDS file to be written.
    If \code{NULL} (unadvised!), will be \code{tempfile}.}
    
    \item{pheno}{The name of the phenotype to be written to
    the respective FAM column.}
    
    \item{outBase}{The full path to the base of the PLINK 
    files triplet. I.e. not with any extension \code{.bim},
    \code{.bed} or \code{.fam}.}
    
    \item{salvage}{try to salvage SNPs with missing locations?
    Defaults to \code{FALSE}. See also Details.}
    
    \item{reverse}{Reverse allele genotypes prior to write?
    Defaults to \code{FALSE}.}
    
    \item{perChr}{Write one PLINK triplet per chromosome?
    Defaults to \code{FALSE}.}
    
    \item{overwrite}{Overwrite previous files with the same
    name? Defaults to \code{TRUE}.}
}
\value{
    This function does not return anything, just writes the
    PLINK files to the desired location.
}
\description{
    This function writes a \code{GWASExperiment} object to
    a triplet of PLINK BID/BED/BAM file(s). A split per
    chromosome is also possible. The \code{pheno} argument
    controls which phenotype will be written (if desired)
    to the respective column of the FAM file.
}
\examples{
    data(toy,package="prisma")
    out <- tempfile()
    writePlink(toy,pheno="BMI",outBase=out)
    
    gdsf <- tempfile()
    writeGdsFile(toy,gdsfile=gdsf)
}
\author{
    Panagiotis Moulos
}
