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
    data(toy,package="prisma")
    
    # Only 1 PC will be returned because of dataset size
    toy <- suppressWarnings(calcPcaCovar(toy,method="snprel"))
    
    toy <- gwa(toy,response="BMI",covariates=c("Age","Gender"),
        pcs=TRUE,psig=0.05,methods=c("glm","statgen"))
    toy <- gwa(toy,response="BMI",covariates=c("Age","Gender"),
        pcs=FALSE,psig=0.05,methods=c("glm","statgen"))
    
    listGwa(toy)
    
    prsiceOpts <- getDefaults("prsice")
    prsiceOpts$perm <- 100
    toy <- runPRS(base=toy,response="BMI",covariates=c("Age","Gender"),
        pcs=TRUE)
    toy <- runPRS(base=toy,response="BMI",covariates=c("Age","Gender"),
        pcs=FALSE)
        
    listPrs(toy)
}
\author{
    Panagiotis Moulos
}

