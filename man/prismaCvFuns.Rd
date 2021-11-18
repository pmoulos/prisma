\name{prismaCvFuns}
\alias{prismaCrossValidate}
\alias{prsCrossValidate}
\title{Functions to perform PRS cross-validation evaluations}
\usage{
    prismaCrossValidate(prismaOut, gwe, response, 
        covariates = NULL, pcs = FALSE, 
        leaveOut = seq(0.05, 0.5, 0.05), times = 10,
        prsCalc = c("avg", "sum", "std"), family = NULL, 
        rc = NULL, ...)
    
    prsCrossValidate(snpSelection, gwe, response, 
        covariates = NULL, pcs = FALSE, leaveOut = 0.05,
        times = 10, prsCalc = c("avg", "sum", "std"),
        family = NULL, rc = NULL, ...)
}
\arguments{
    \item{prismaOut}{a \code{GWASExperiment} object.}
    
    \item{gwe}{a \code{GWASExperiment} object (same as the one
    used to generate \code{prismaOut}).}
    
    \item{response}{}
    
    \item{covariates}{}
    
    \item{leaveOut}{}
    
    \item{times}{}
    
    \item{prsCalc}{}
    
    \item{family}{}
    
    \item{pcs}{}
    
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

