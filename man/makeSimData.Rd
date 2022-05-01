\name{makeSimData}
\alias{makeSimData}
\title{Function to create basic simulated GWAS data}
\usage{
    makeSimData(nsnp = 10000, nsam = 100, nphe = 3,
        csnp = 100, psnp = c(0.01, 0.05, 0.1, 0.15, 0.20, 0.25,
        0.3, 0.4))
}
\arguments{
    \item{nsnp}{The number of SNPs to simulate genotypes 
    for. Defaults to \code{10000}.}
    
    \item{nsam}{The number of samples to simulate. Defaults
    to \code{100}.}
    
    \item{nphe}{The number of continuous phenotypes to
    simulate. Defaults to \code{3}.}
    
    \item{csnp}{The number of causative SNPs to simulate.
    Defaults to 100.}
    
    \item{psnp}{The MAFs to sample from.}
}
\value{
    A list with members:
    \itemize{
        \item \code{snp}: A SnpMatrix object with genotypes.
        \item \code{sample}: The sample map (\code{gsamples} or
        the equivalent of a FAM PLINK file).
        \item \code{feature}: The feature map (\code{gfeatures} 
        or the equivalent of a BIM PLINK file).
        \item \code{pheno}: The simulated data frame of 
        phenotypes.
    }
    This list can be used to create a \code{GWASExperiment}
    object.
}
\description{
    This function is losely based on functionalities from
    the R package PhenotypeSimulator and can be used to 
    create a basic simulated GWAS dataset with continuous
    phenotypes.
}
\examples{
    dat <- makeSimData(nsnp=100,nsam=100,nphe=3,csnp=10)
    gwe <- GWASExperiment(
        genotypes=dat$snp,
        features=dat$feature,
        samples=dat$sample,
        phenotypes=dat$pheno
    )
}
\author{
    Panagiotis Moulos
}
