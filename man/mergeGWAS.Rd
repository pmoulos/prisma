\name{mergeGWAS}
\title{Merge GWASExperiments from different populations}
\usage{
    mergeGWAS(gwe1, gwe2, gdsfile = NULL, writegds = FALSE)
}
\arguments{
    \item{gwe1}{A \code{GWASExperiment} object from the first
    dataset.}
    
    \item{gwe1}{A \code{GWASExperiment} object from the second
    dataset.}
    
    \item{gdsfile}{path to the new GDS file. If \code{NULL}, a 
    temporary file will be written (unadvised!).}
    
    \item{writegds}{write the new GDS file? Defaults to 
    \code{FALSE}. It can be written anytime with the
    \code{\link{GWASExperiment2GDS}} function.}
}
\value{
    A new GWASExperiment object with merged genotypes,
    samples and phenotypes.
}
\description{
    This function attempts to merge two GWASExperiments with
    various goals, such as merging similar populations with
    common phenotypic and/or anthropometric data to increase
    statistical power. See Details for a more rigorous 
    explanation of the merging process.
}
\details{
    
}
\examples{
    
}
\author{
    Panagiotis Moulos
}
