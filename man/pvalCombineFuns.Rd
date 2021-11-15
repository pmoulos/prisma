\name{pvalCombineFuns}
\alias{combineHarmonic}
\alias{combineMaxp}
\alias{combineMinp}
\alias{combineSimes}
\alias{combineWeight}
\alias{fisherMethod}
\title{GWA test p-value combination functions}
\usage{
    combineHarmonic(p, w = NULL, zerofix = NULL)
    combineMaxp(p)
    combineMinp(p)
    combineSimes(p)
    combineWeight(p, w = NULL, zerofix = NULL)
    
    fisherMethod(p, zerofix = NULL)
}
\arguments{
    \item{p}{a numeric vector with p-values or a matrix
    for \code{fisherMethod}.}
    
    \item{w}{a weights vector, must sum to 1.}
    
    \item{zerofix}{\code{NULL} (default) or a fixed 
    numeric value between 0 and 1.}
}
\value{
    A new p-value, derived from the combination of 
    \code{p}.
}
\description{
    This group of fuctions can be used to calculate
    combined p-values, where each p-value is typically
    derived from a genome-wide association test for a
    SNP or a group of SNPs against a phenotype.
}
\details{
    The \code{combineHarmonic} function combines a vector
    of p-values using the method of Wilson, 2019 
    \link{https://www.pnas.org/content/116/4/1195} and
    user- or prisma- provided weights.
    
    The \code{combineMaxp}, \code{combineMinp} functions 
    combine a vector of p-values using the maximum p-value
    (equivalent to the interesection of lists of significant
    SNPs returned by several GWA tests) or using the
    minimum p-value (equivalent to the union of lists of 
    significant SNPs returned by several GWA tests).
    
    The \code{combineSimes} function combines a vector
    of p-values using the method of Simes, 1986 method.
    
    The \code{combineWeight} function combines a vector
    of p-values using the PANDORA algorithm described by
    Moulos and Hatzis, 2015 (PMID: 25452340) and Fanidis
    and Moulos, 2020 (PMID: 32778872), which is based on
    weighted p-value combination based on simulations from
    real data.
    
    The \code{fisherMethod} function combines a vector
    of p-values using the respective method of Fisher. 
    
    The argument \code{zerofix} is used to correct for
    the case of a p-value which is equal to 0 as a result
    of internal numerical and approximation procedures.
    When \code{NULL}, random numbers greater than 0 and
    less than or equal to 0.5 are used to multiply the
    offending p-values with the lowest provided non-zero
    p-value, maintaining thus a virtual order of 
    significance, avoiding having the same p-values for 
    two tests and assuming that all zero p-values represent
    extreme statistical significance. When a numeric
    between 0 and 1, this number is used for the above
    multiplication instead.
}
\examples{
    w <- c(0.2,0.5,0.3)
    p <- runif(3)
    
    pha <- combineHarmonic(p,w)
    pma <- combineMax(p)
    pmi <- combineMin(p)
    psi <- combineSimes(p)
    pwe <- combineWeight(p,w)
    
    pm <- matrix(runif(300),100,3)
    
    pha <- apply(pm,1,combineHarmonic,w)
    pma <- apply(pm,1,combineMax)
    pmi <- apply(pm,1,combineMin)
    psi <- apply(pm,1,combineSimes)
    pwe <- apply(pm,1,combineWeight,w)
    
    pwf <- fisherMethod(pm)
}
\author{
    Panagiotis Moulos
}

