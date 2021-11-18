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
    
    \item{response}{}
    
    \item{covariates}{}
    
    \item{pcs}{}
    
    \item{psig}{}
    
    \item{methods}{}
    
    \item{combine}{}
    
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
    
}
\details{
    
}
\examples{
    
}
\author{
    Panagiotis Moulos
}

