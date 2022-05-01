\name{prismaCvFuns}
\alias{prismaCrossValidate}
\alias{prsCrossValidate}
\title{Functions to perform PRS cross-validation evaluations}
\usage{
    prismaCrossValidate(prismaOut, gwe, response, 
        covariates = NULL, pcs = FALSE, 
        leaveOut = seq(0.05, 0.5, 0.05), times = 10,
        prsCalc = c("avg", "sum", "std"), family = NULL, 
        rc = NULL, ...)
    
    prsCrossValidate(snpSelection, gwe, response, 
        covariates = NULL, pcs = FALSE, leaveOut = 0.05,
        times = 10, prsCalc = c("avg", "sum", "std"),
        family = NULL, rc = NULL, ...)
}
\arguments{
    \item{prismaOut}{an output from the \code{link{prisma}}
    function containing the results of an iterative PRS
    extraction.}
    
    \item{gwe}{a \code{GWASExperiment} object (same as the one
    used to generate \code{prismaOut}).}
    
    \item{snpSelection}{A data frame with candidate SNPs, the
    output, or part of the output of 
    \code{\link{aggregatePrsMarkers}}, or one of the list 
    members of the output of \code{\link{selectPrs}}.}
    
    \item{response}{The trait for which a PRS is sought to be
    extracted. It should be a character of length 1 and exist
    in the phenotypes of \code{gwe}.}
    
    \item{covariates}{The covariates to include in regression
    models. It should be a character vector of traits all 
    existing in the phenotypes of \code{gwe}.}
    
    \item{leaveOut}{The fraction (percentage) of samples to
    leave-out in a cross-validation run. This can be a vector
    of multiple fractions for \code{prismaCrossValidate}.}
    
    \item{times}{How many cross-validations to perform?
    Defaults to \code{10}.}
    
    \item{prsCalc}{How the final PRS is calculated from the
    data frame of the PRS candidate SNPs? One of \code{"avg"} 
    (default), \code{"sum"} or \code{"std"}. See also 
    Details.}
    
    \item{family}{Regression family of functions, e.g. 
    \code{"gaussian"} for classical linear regression or
    \code{"binomial"} for logistic regression. See also
    \code{\link{gwa}}. Defaults to \code{NULL} for trying to
    set automatically.}
    
    \item{pcs}{Inlcude PCs in the regression model? Defaults
    to \code{FALSE}. See also Details.}
    
    \item{rc}{The fraction of available cores to use for 
    parallel calculations. Default is \code{NULL} for no
    parallelization.}
    
    \item{...}{Further arguments passed to regression 
    functions, like \code{glm}.}
}
\value{
    For \code{prismaCrossValidate}, a list of 
    cross-validation metrics, for each GWA method used, for
    each PRS candidate and for each leave-out percentage.
    The list is nested, with the first level being the GWA 
    method, the second level being the SNP candidates (named
    after the number of SNPs in the PRS candidate) and the
    third level being the leave-out percentage of samples
    in cross-validation. The actual value of the third
    level is a data frame with \code{times} rows, containing
    the metrics collected after each cross-validation,
    leaving out \code{leaveOut[i]} fraction of samples,
    with \code{i=1:ntimes}. This data frame is the output
    of \code{prsCrossValidate}. For the collected metrics,
    see Details.
}
\description{
    These functions perform cross-validation regressions
    given either the output of \code{\link{prisma}} for
    full evaluation (\code{prismaCrossValidate}) or one
    PRS candidate data frame with SNPs in proper format
    (\code{prsCrossValidate}). Several metrics are 
    returned for inspection and evaluation and are all
    used in the PRISMA HTML report.
}
\details{
    Regarding the \code{prsCalc} argument, this tells the CV
    functions how to calculate the PRS in the model given the
    genotypic data. The three available method reflect the 
    three first ways of calculating the PRS in the PRSice2
    software (https://www.prsice.info/step_by_step/) and
    are named accordingly (see the \code{--score} option in
    PRSice2).
    
    Regarding the \code{pcs} argument, if \code{TRUE}, the PCs 
    in \code{gwe} will be used. If \code{gwe} does not contain 
    any PCs, a warning is issued and \code{pcs} is set to 
    \code{FALSE}.
    
    Regarding the collected metrics, these are:
    \itemize{
        \item R^2 of the reduced regression model 
        (\code{reduced_r2}): The R^2 of the reduced
        regression model, including all the covariates
        and PCs if requested, but not the PRS.
        \item R^2 of the full regression model 
        (\code{full_r2}): The R^2 of the full regression
        model including all the covariates and the PRS.
        \item R^2 of the PRS (\code{prs_r2}): The adjusted
        R^2 of the PRS (essentially R^2 of full model minus
        the R^2 of the reduced model).
        \item Statistical significance (p-value) of the PRS 
        effect coefficient in the regression model
        (\code{prs_pvalue}).
        \item R^2 of the reduced regression model 
        (\code{reduced_r2}): The R^2 of the reduced
        regression model, including all the covariates
        and PCs if requested, but not the PRS.
        \item RMSE (Root Mean Squared Error) of the reduced 
        regression model (\code{reduced_rmse}): The RMSE of 
        the reduced regression model, including all the 
        covariates and PCs if requested, but not the PRS.
        \item RMSE of the full regression model 
        (\code{full_rmse}): The RMSE of the full regression
        model including all the covariates and the PRS.
        \item MAE (Mean Absolute Error) of the reduced 
        regression model (\code{reduced_mae}): The MAE of 
        the reduced regression model, including all the 
        covariates and PCs if requested, but not the PRS.
        \item MAE of the full regression model 
        (\code{full_mae}): The MAE of the full regression
        model including all the covariates and the PRS.
        \item Correlation between observed and predicted
        phenotype values for the reduced regression model
        (\code{reduced_pred_cor}).
        \item Correlation between observed and predicted
        phenotype values for the full regression model
        (\code{full_pred_cor}).
        \item R^2 for the observed and predicted phenotype
        values for the reduced regression model
        (\code{reduced_pred_r2}).
        \item R^2 for the observed and predicted phenotype
        values for the full regression model
        (\code{full_pred_r2}).
        \item R^2 of the PRS for the observed and predicted 
        phenotype values for the reduced regression model
        (\code{prs_pred_r2}).
    }
}
\examples{
    data(toy,package="prisma")
    data(prisma_out,package="prisma")
    
    response <- "BMI"
    covariates <- c("Case_Control","Gender","Age")
    
    cvm <- prismaCrossValidate(prismaOut,toy,response=response,
        covariates=covariates,leaveOut=c(0.2,0.4),times=5)
    
    can <- prismaOut$results$glm$candidates[[1]]
    # can <- getPrsCandidates(prismaOut,"glm",1)
    m <- prsCrossValidate(can,toy,response=response,
        covariates=covariates,leaveOut=0.2,times=5)
}
\author{
    Panagiotis Moulos
}
