\name{adjustPrsWeights}
\alias{adjustPrsWeights}
\title{PRS weight adjustment through final regression}
\usage{
    adjustPrsWeights(snpSelection, gwe, response, 
        covariates = NULL, pcs = FALSE, family = NULL,
        rc = NULL)
}
\arguments{
    \item{snpSelection}{A data frame with SNPs comprising
    a PRS, properly formatted (see 
    \code{link{aggregatePrsMarkers})}.}
    
    \item{gwe}{The \code{GWASExperiment} object used to
    derive the PRS.}
    
    \item{response}{A character with the phenotypic response
    used to derive the PRS.}
    
    \item{covariates}{A character vector with covariares
    used to derive the PRS.}
    
    \item{pcs}{Use the PCs present in \code{gwe}? Defaults
    to \code{FALSE}. If the PRS has been extracted with PCs
    as covariates, it should be \code{TRUE}.}
    
    \item{family}{Regression family, e.g. \code{"gaussian"},
    see also \code{\link{gwa}}.}
    
    \item{rc}{The fraction of available cores to use for 
    parallel calculations. Default is \code{NULL} for no
    parallelization.}
}
\value{
    A new data frame with SNPs comprising the PRS with
    adjusted effects.
}
\description{
    This function may be used to adjust the SNP effects
    to derive the final PRS (e.g. with the total dataset)
    after the whole extraction and assessment process.
    The process is not obligatory.
}
\examples{
    data(toy,package="prisma")
    data(prisma_out,package="prisma")
    
    prs <- getPrsCandidates(prismaOut,"statgen",1)
    #prs <- prismaOut$results$statgen$candidates[[1]]
    
    prsa <- adjustPrsWeights(prs,toy,"BMI",
        c("Age","Gender","Case_Control"),pcs=FALSE)
    
    identical(prs$effect_weight,prsa$effect_weight)
    # FALSE
}
\author{
    Panagiotis Moulos
}
