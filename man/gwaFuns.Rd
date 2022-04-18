\name{gwaFuns}
\alias{gwa}
\alias{gwaBlup}
\alias{gwaGlm}
\alias{gwaPlink}
\alias{gwaSnptest}
\alias{gwaStatgen}
\title{Functions to perform GWA tests with various algorithms}
\usage{
    gwa(obj, response, covariates = NULL, pcs = FALSE,
        psig = 0.05, methods = c("glm", "rrblup", "statgen",
            "snptest", "plink", "lasso"),
        combine = c("fisher", "simes", "max", "min", "harmonic",
            "whitlock","pandora"), glmOpts = getDefaults("glm"),
        rrblupOpts = getDefaults("rrblup"), 
        statgenOpts = getDefaults("statgen"),
        snptestOpts = getDefaults("snptest"),
        plinkOpts = getDefaults("plink"), rc = NULL, ...)
    
    gwaBlup(obj, response, covariates = NULL, 
        usepc = c("auto", "estim", "rrint", "fixed", "none"), 
        npcs = NULL, psig = 0.05, rc = NULL)
    
    gwaGlm(obj, response, covariates = NULL, pcs = FALSE,
        family = NULL, psig = 0.05, penalized = FALSE,
        size = 1, rc = NULL, ...)
        
    gwaPlink(obj, response, covariates = NULL, pcs = TRUE,
        psig = 0.05, seed = 42, wspace = NULL, rc = NULL)
    
    gwaStatgen(obj, response, covariates = NULL, pcs = TRUE, 
        psig = 0.05, rc = NULL)
    
    gwaSnptest(obj, response, covariates = NULL, pcs = TRUE,
        psig = 0.05, test = c("frequentist", "bayesian"), 
        model = c("additive", "dominant", "recessive", 
        "general", "heterozygote"), wspace = NULL)
}
\arguments{
    \item{obj}{an object of class \code{GWASExperiment} holding
    the genotypic and phenotypic data to be used for association
    tests.}
    
    \item{response}{The trait for which associations are 
    searched for. It should be a character of length 1 and 
    exist in the phenotypes of \code{obj}.}
    
    \item{covariates}{The covariates to include in association
    methods/regression models. It should be a character vector 
    of traits all existing in the phenotypes of \code{obj}.}
    
    \item{pcs}{Inlcude PCs in the association models? Defaults
    to \code{FALSE}.}
    
    \item{psig}{The significance level (p-value cutoff) for which
    an association between a SNP/marker and a phenotype 
    (\code{response}) is considered significant. Defaults to
    \code{0.05}.}
    
    \item{methods}{The GWA methods to apply on \code{obj} with
    respect to \code{response} and \code{covariates}. It can be
    one or more of \code{"glm"}, \code{"rrblup"}, \code{"statgen"},
    \code{"snptest"}, \code{"plink"}, \code{"lasso"}. See also
    Details for explanation of each case.}
    
    \item{combine}{How (if) p-values should be combined to
    produce a summarized p-value from all GWA methods for
    each association test. Not all methods are yet available,
    see also Details.}
    
    \item{glmOpts}{}
    
    \item{rrblupOpts}{}
    
    \item{statgenOpts}{}
    
    \item{snptestOpts}{}
    
    \item{plinkOpts}{}
    
    \item{usepc}{}
    
    \item{npcs}{}
    
    \item{family}{}
    
    \item{wspace}{}
    
    \item{test}{}
    
    \item{model}{}
    
    \item{rc}{The fraction of available cores to use for 
    parallel calculations. Default is \code{NULL} for no
    parallelization.}
}
\value{
    Depending on \code{out}
}
\description{
    These functions perform the association tests to create 
    summary statistics for a GWA dataset prior to the 
    extraction of PRS candidates. They are part of the 
    overall PRISMA workflow (\code{\link{prisma}}) but they
    are designed to be also used independently just to look
    for univariate association cases, i.e. look for 
    associations for a single trait using various tests.
}
\details{
    
}
\examples{
    
}
\author{
    Panagiotis Moulos
}

