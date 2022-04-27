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
    \item{y}{A vector of numeric values to be split to a
    "training" and a "test" set and return the respective
    indices.}
    
    \item{obj}{Instead of \code{y} above, a GWASExperiment
    object.}
    
    \item{by}{The phenotype in \code{obj} to be used for
    splitting. It should be a scalar character.}
    
    \item{n}{How many splits to create? Defaults to \code{1}.}
    
    \item{frac}{The fraction of values to remain in the
    "training" set. Defaults to \code{0.5}.}
    
    \item{ngrp}{In how many "virtual" groups should \code{y}
    be split? Controls the uniformity of values in the two
    resulting sets. See also Details.}
    
    \item{replace}{Sample values/samples with replacement?
    Defaults to \code{FALSE}.}
    
    \item{out}{For \code{createSplit}, either \code{"index"}
    (default) to return a list of length \code{n} with 
    "training" index splits, or \code{"binary"} to return
    the same list but instead of indices, it contains a
    mask of \code{0}s and \code{1}s indicating which samples
    have been selected and which discarded in the "training"
    set. For \code{partitionGWAS}, either \code{"index"} 
    (default and same as \code{createSplit}m but now the
    indices refer to the columns of \code{obj}), or 
    \code{"train"} to return a \code{GWASExperiment} with
    the training samples, or \code{"ttboth"} to return
    a list of \code{GWASExperiment}s. In the latter case,
    the list has named members with names \code{train} and 
    \code{test} respectively and the value of \code{n} is
    ignored, returning only one split.}
    
    \item{rc}{The fraction of available cores to use for 
    parallel calculations. Default is \code{NULL} for no
    parallelization.}
}
\value{
    For \code{createSplit}, a list of length \code{n}
    with the splits in the form of vector indices or
    a binary mask with \code{1}s the training samples
    and \code{0}s the test samples.
    
    For \code{partitionGWAS}, a list of indices as above
    or a list of \code{GWASExperiment} objects to be used
    as training samples, or a list with the training and
    testing \code{GWASExperiment} objects.
}
\description{
    These functions are used to create training and testing
    splits to be used in the iterative PRS extraction process
    of PRISMA. They are losely based on the dataset 
    partitoning functions of the R package caret and can also 
    be used independently of the general pipeline.
}
\examples{
    data(toy,package="prisma")
    
    trainTest <- partitionGWAS(toy,by="BMI",out="ttboth")
    
    y <- phenotypes(toy)$BMI
    z <- createSplit(y,frac=0.8,n=5)
}
\author{
    Panagiotis Moulos
}
