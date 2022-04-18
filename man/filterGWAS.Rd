\name{filterGWAS}
\title{Aggregate PRS candidate SNPs from iterative extraction}
\usage{
    filterGWAS(obj, filters = getDefaults("filters"),
        imputeMissing = TRUE, imputeMode = c("single","split"),
        rc = NULL)
}
\arguments{
    \item{obj}{An object of class \code{GWASExperiment} holding
    the genotypic data to be filtered.}
    
    \item{filters}{A set of genotype filters, see also Details
    and \code{\link{getDefaults}}.}
    
    \item{imputeMissing}{Impute the remaining missing values
    after filtering? An iterative imputation method based
    on the package snpStats is used, see also 
    \code{\link{imputeGWAS}}.}
    
    \item{imputeMode}{Perform the imputation within the same
    dataset (\code{"single"}, default) or split the dataset
    into a part containing missing values and a part with
    no missing values (\code{"split"}), may happen in large
    dataset and could be more accurate.}
    
    \item{rc}{The fraction of available cores to use for 
    parallel calculations. Default is \code{NULL} for no
    parallelization.}
}
\value{
    A filtered \code{GWASExperiment} object with or without
    missing values according to \code{imputeMissing}.
}
\description{
    This function implements various genotype and sample
    filters based on properties of the genotypic data
    present in the input GWAS object. Most filters are
    losely based in the protocol by Reed et al., 2015
    (PubMed ID: 26343929).
}
\details{
    Regarding the \code{filters} argument, this is a list
    whose members denote SNP or sample-based filters and
    their values are filter cutoffs. A more detailed
    explanation of these filters can be found in Reed et 
    al., 2015 (PubMed ID: 26343929).
    
    The following filters are available:
    \itemize{
        \item \code{snpCallRate}: The SNPs call rate within
        a sample. Defaults to \code{0.98}.
        \item \code{sampleCallRate}: The sample call rate of
        SNPs accross samples. Defaults to \code{0.95}.
        \item \code{maf}: The Minimum Allele Frequency across
        samples, below which SNPs are filtered. Defaults to 
        \code{0.05}.
        \item \code{hwe}: The Hardy-Weinberg equilibrium 
        p-value. Defaults to \code{1e-6}.
        \item \code{heteroStat}: The statistic with which
        heterozygosity distribution center is defined. 
        Defaults to \code{"median"}.
        \item \code{heteroFac}: The factor with which the 
        scale (standard deviation, MAD) of the heterozygosity
        distribution is multiplied and added before and after
        the \code{heteroStat} to create a range of acceptable
        heterozygosity values, e.g. mean +/- 3*std. Defaults
        to \code{3}. If \code{heteroStat="mean"}, the scale
        is the standard deviation. If 
        \code{heteroStat="median"}, the scale is the MAD.
        \item \code{heteroHard}: A hard filter for 
        heterozygosity values. If set, \code{heteroStat} and
        \code{heteroFac} are ignored.
        \item \code{pcaOut}: Automatically remove samples 
        based on PCA filtering? Defaults to \code{TRUE} but
        should be avoided for datasets with a few samples.
        Also, it is applied if a robust PCA method is chosen.
        \item \code{pcaRobust}: The robust PCA method to be
        applied for sample exclustion. One of \code{"grid"}
        or \code{"hubert"} (default). See the package rrcov
        for more details.
        \item \code{nPC}: Fixed number of PCs to use for PCA
        filtering. If \code{NA} (default), they are 
        determined using the \code{\link{twTest}.}
        \item \code{LD}: Linkage Disequilibrium filter 
        (calculated using package SNPRelate fuctions).
        Defaults to \code{0.2}.
        \item \code{IBD}: Identity By Descent filter 
        (calculated using package SNPRelate fuctions).
        Defaults to \code{0.2}.
        \item \code{inbreed}: Inbreed coefficient filter.
        Defaults to \code{0.1}.
    }
    
    See also the output of \code{getDefaults("filters")}.
}
\examples{
    data(toy,package="prisma")

    # Lose filters for the toy dataset
    filts <- getDefaults("filters")
    filts$IBD <- NA
    filts$hwe <- 1e-3
    filts$pcaOut <- FALSE
    filts$inbreed <- NA
    
    toyF <- filterGWAS(toy,filters=filts,
        imputeMissing=FALSE)
}
\author{
    Panagiotis Moulos
}
