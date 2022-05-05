\name{prsRegressionMetrics}
\alias{prsRegressionMetrics}
\title{Calculate PRS evaluation metrics}
\usage{
    prsRegressionMetrics(snpSelection, gwe, response, 
        covariates = NULL, pcs = FALSE, step = 10, family = NULL,
        rc = NULL, ...)
}
\arguments{
    \item{snpSelection}{A data frame with candidate SNPs, the
    output, or part of the output of 
    \code{\link{aggregatePrsMarkers}}, or one of the list 
    members of the output of \code{\link{selectPrs}}.}
    
    \item{gwe}{A \code{GWASExperiment} object (same as the one
    used to generate \code{prismaOut}).}
    
    \item{response}{The trait for which a PRS is sought to be
    extracted. It should be a character of length 1 and exist
    in the phenotypes of \code{gwe}.}
    
    \item{covariates}{The covariates to include in regression
    models. It should be a character vector of traits all 
    existing in the phenotypes of \code{gwe}.}
    
    \item{pcs}{Include PCs from \code{gwe} in the regression
    models? Defaults to \code{FALSE}.}
    
    \item{step}{Increment step for the number of SNPs from
    \code{snpSelection} to include in iterative evaluation.
    Defaults to \code{10}.}
    
    \item{family}{Regression family of functions, e.g. 
    \code{"gaussian"} for classical linear regression or
    \code{"binomial"} for logistic regression. See also
    \code{\link{gwa}}. Defaults to \code{NULL} for trying to
    set automatically.}
    
    \item{rc}{The fraction of available cores to use for 
    parallel calculations. Default is \code{NULL} for no
    parallelization.}
    
    \item{...}{Further arguments passed to \code{glm}
    regression functions.}
}
\value{
    A data frame with regression fit metrics, containing the
    following columns:
    \itemize{
        \item \code{n_snp}: The number of SNPs in the PRS 
        being evaluated
        \item \code{freq}: The frequency of appearance for
        these SNPs
        \item \code{full_r2}: The R^2 of the full regression
        model, i.e. the regression model including the 
        covariates and the PRS
        \item \code{full_pvalue}: The p-value of the full
        regression model
        \item \code{prs_pvalue}: The p-value of the PRS
        coefficient in the model
        \item \code{full_aic}: The AIC of the full 
        regression model
        \item \code{reduced_r2}: The R^2 of the reduced
        model, i.e. the regression model including only 
        the covariates and excluding the PRS
        \item \code{prs_r2}: The R^2 attributed to the PRS
        (\code{full_r2 - reduced_r2})
        \item \code{reduced_pvalue}: The p-value of the 
        reduced regression model
        \item \code{reduced_aic}: The AIC of the reduced
        regression model
        \item \code{prs_aic}: The AIC of the PRS 
        (\code{full_aic - reduced_aic})
    }
}
\description{
    This function fits a regression model given a number of
    SNPs from which a PRS is calculated, in an iterative 
    manner. Then, several metrics for goodness of fit are
    calculated and returned.
}
\examples{
    data(toy,package="prisma")
    data(prisma_out,package="prisma")
    
    response <- "BMI"
    covariates <- c("Case_Control","Gender","Age")
    
    can <- getPrsCandidates(prismaOut,"glm",1)
    m <- prsRegressionMetrics(can,toy,response=response,
        covariates=covariates,pcs=FALSE)
}
\author{
    Panagiotis Moulos
}

