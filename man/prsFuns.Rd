\name{prsFuns}
\alias{lassosumPRS}
\alias{prsicePRS}
\alias{runPRS}
\title{Functions to extract PRS candidates with various algorithms}
\usage{
    lassosumPRS(base, target = base, response, covariates = NULL,
        pcs = FALSE, wspace = NULL, anc = c("eur", "asn", "afr"),
        valid = c("auto", "stack", "split"))
    
    prsicePRS(base, target = base, response, covariates = NULL,
        pcs = FALSE, mode = c("calculate", "apply"), wspace = NULL,
        clump_kb = 250, clump_r2 = 0.1, clump_p = 1,
        score = c("avg", "sum", "std", "con-std"), perm = 10000,
        seed = 42, rc = NULL)
    
    runPRS(base, target = base, response, covariates = NULL,
        pcs = FALSE, methods = c("lassosum", "prsice"),
        lassosumOpts = getDefaults("lassosum"),
        prsiceOpts = getDefaults("prsice"), wspace = NULL, 
        rc = NULL)
}
\arguments{
    \item{base}{}
    
    \item{target}{}
    
    \item{response}{}
    
    \item{covariates}{}
    
    \item{pcs}{}
    
    \item{wspace}{}
    
    \item{anc}{}
    
    \item{valid}{}
    
    \item{clump_kb}{}
    
    \item{clump_r2}{}
    
    \item{clump_p}{}
    
    \item{score}{}
    
    \item{perm}{}
    
    \item{seed}{}
    
    \item{methods}{}
    
    \item{lassosumOpts}{}
    
    \item{prsiceOpts}{}
    
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

