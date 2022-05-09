\name{gweConvertFuns}
\alias{GWASExperiment2gData}
\alias{GWASExperiment2GDS}
\alias{GWASExperimentLiftOver}
\alias{guessHumanGenomeVersion}
\title{Functions to convert GWASExperiment to other classes}
\usage{
    GWASExperiment2gData(obj, covariates = NULL, pcs = FALSE,
        reverse = FALSE)
    
    GWASExperiment2GDS(obj)
    
    GWASExperimentLiftOver(obj, from, to)
    
    guessHumanGenomeVersion(obj)
}
\arguments{
    \item{obj}{a \code{GWASExperiment} object.}
    
    \item{covariates}{the covariates to include in the
    \code{gData} object (package statgenGWAS) including the
    response phenotype.}
    
    \item{pcs}{include the PCs present in \code{obj} in the
    \code{gData} object? Defaults to \code{FALSE}.}
    
    \item{reverse}{reverse the genotypes in the construction
    of a \code{gData} object? Defaults to \code{FALSE}.}
    
    \item{from}{the source genome version, one of 
    \code{"hg19"} or \code{"hg38"}.}
    
    \item{to}{the target genome version, one of 
    \code{"hg19"} or \code{"hg38"}.}
}
\value{
    For \code{GWASExperiment2gData}, a \code{gData} object
    that can be used with functions from the statgenGWAS 
    package.
    For \code{GWASExperiment2GDS}, a GDS file is written to
    the location specified by \code{gdsfile(obj).}
    For \code{GWASExperimentLiftOver}, the same 
    GWASExperiment object is returned, with changed SNP 
    positions from one genome version to the other. Some
    SNPs may be dropped because of inability of the liftOver
    algorithm to map between versions (rare).
    For \code{guessHumanGenomeVersion}, a string denoting
    the guessed human genome version (\code{"hg19"} or 
    \code{"hg38"}).
}
\description{
    These functions either convert a GWASExperiment to 
    another object or file required for functionalities 
    in the supporting packages, or convert between human
    genome coordinate systems. Finally, to facilitate this
    conversion process when required, a functionality to
    "guess" the human genome version is provided.
}
\details{
    Regarding the \code{guessHumanGenomeVersion} function,
    it queries dbSNP for 5 random SNPs present in \code{obj}
    using the package rsnps. Depending on the result, the
    genome version is guessed.
}
\examples{
    data(toy,package="prisma")
    
    # There are some non-numeric phenotypes - ignore
    gd <- suppressWarnings(GWASExperiment2gData(toy))
    
    gdsfile(toy) <- tempfile()
    GWASExperiment2GDS(toy)
    file.exists(gdsfile(toy)) # TRUE
    
    genome(toy) # hg19
    toy2 <- GWASExperimentLiftOver(toy,from="hg19",
        to="hg38")
    
    guessHumanGenomeVersion(toy) # hg19
}
\author{
    Panagiotis Moulos
}
