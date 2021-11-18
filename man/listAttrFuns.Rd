\name{listAttrFuns}
\alias{listGwa}
\alias{listPrs}
\title{Functions to list analyses in a GWASExperiment}
\usage{
    listGwa(obj)
    listPrs(obj)
}
\arguments{
    \item{obj}{a \code{GWASExperiment} object.}
}
\value{
    These functions only print information on screen.
}
\description{
    These functions list the associated analyses with a
    GWASExperiment, i.e. GWA tests for association summary
    statistics, or PRS weights from PRS extraction analysis.
    A GWASExperiment may have multiple GWA and PRS analysis
    performed on it, such as GWA with multiple or no PCs as
    covariates. Each analysis is attached to the main
    GWASExperiment to avoid repetition of large genotypic
    data.
}
\examples{
    
}
\author{
    Panagiotis Moulos
}

