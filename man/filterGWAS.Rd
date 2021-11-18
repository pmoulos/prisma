\name{filterGWAS}
\title{Aggregate PRS candidate SNPs from iterative extraction}
\usage{
    filterGWAS(obj, filters = getDefaults("filters"),
        imputeMissing = TRUE, imputeMode = c("single","split"),
        rc = NULL)
}
\arguments{
    \item{obj}{an object of class \code{GWASExperiment} holding
    the genotypic data to be filtered.}
    
    \item{filters}{a set of genotype filters, see also Details
    and \code{\link{getDefaults}}.}
    
    \item{imputeMissing}{}
    
    \item{imputeMode}{}
    
    \item{rc}{The fraction of available cores to use for 
    parallel calculations. Default is \code{NULL} for no
    parallelization.}
}
\value{
    
}
\description{
    
}
\details{
    
}
\examples{
    
}
\author{
    Panagiotis Moulos
}

