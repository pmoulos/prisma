\name{downloadExternalFuns}
\alias{downloadExternalTools}
\alias{downloadGtool}
\alias{downloadImpute}
\alias{downloadPlink}
\alias{downloadPrsice}
\alias{downloadSnptest}
\alias{downloadQctool}
\alias{download1000GP3}
\title{Functions to download external GWA and PRS tools}
\usage{
    downloadExternalTools(tools = c("snptest", "plink", "prsice"),
        ver = getDefaults("externalTools"), path = NULL)
    
    downloadPlink(ver = c("v1.90"), path = NULL)
    
    downloadPrsice(ver = c("v2.3.3"), path = NULL) 
    
    downloadSnptest( ver = c("v2.5.6", "v2.5.4", "v2.5.2"),
        path = NULL)
    
    
    download1000GP3(path = NULL)
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
    These functions do not return a value, only update the
    environment holding paths to these tools. 
}
\description{
    After the download and installation process is complete 
    (if any), these funcitons update the environmental variable 
    (an R environment essentially)
}
\details{
    
}
\examples{
    
}
\author{
    Panagiotis Moulos
}

