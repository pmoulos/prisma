\name{splitFuns}
\alias{createSplit}
\alias{partitionGWAS}
\title{Functions to split datasets in training and test samples}
\usage{
    createSplit(y, n = 1, frac = 0.5, ngrp = min(5, length(y)),
        replace = FALSE, out = c("index", "binary"), rc = NULL)
    
    partitionGWAS(obj, by, n = 1, frac = 0.5, 
        ngrp = min(5, ncol(obj)), replace = FALSE,
        out = c("index", "train", "ttboth"), rc = NULL)
}
\arguments{
    \item{y}{}
    
    \item{obj}{}
    
    \item{by}{}
    
    \item{n}{}
    
    \item{frac}{}
    
    \item{ngrp}{}
    
    \item{replace}{}
    
    \item{out}{}
    
    \item{rc}{The fraction of available cores to use for 
    parallel calculations. Default is \code{NULL} for no
    parallelization.}
}
\value{
    Depending on \code{out}
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

