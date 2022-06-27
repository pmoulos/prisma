\name{prismaCvFuns}
\alias{prismaCrossValidate}
\alias{prsCrossValidate}
\alias{applyPRS}
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
    
    applyPRS(snpSelection, gwe, response, 
        covariates = NULL, pcs = FALSE, minSnps = 2,
        prsCalc = c("avg", "sum", "std"), times = 10,
        family = NULL, summary = TRUE, rc = NULL,...)
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
    
    \item{times}{For \code{prismaCrossValidate}, how many 
    cross-validations to perform. For \code{applyPRS}, the
    number of SNP randomizations to perform if summary 
    statistics are present to \code{gwe}. Defaults to 
    \code{10}. See also Details.}
    
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
    
    \item{minSnps}{The minimum number of SNPs present in
    \code{gwe} in order to continue with PRS application
    using the genotypes of \code{gwe}. Defaults to 
    \code{2}.}
    
    \item{rc}{The fraction of available cores to use for 
    parallel calculations. Default is \code{NULL} for no
    parallelization.}
    
    \item{summary}{Display the evaluation summarized metrics
    in more human readable format and make plots.}
    
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
    
    For \code{applyPRS} a list with four or five members:
    \itemize{
        \item \code{prs}: The input PRS (\code{snpSelection}
        data frame).
        \item \code{metrics}: The same as the cross
        validation metrics - see Details.
        \item \code{values}: A data frame with the following
        columns:
        \itemize{
            \item \code{raw_pheno}: The phenotype which
            the PRS contributes to its explanation
            \item \code{red_pheno}: The predicted phenotype
            in a regression model without the PRS as a
            covariate
            \item \code{full_pheno}: The predicted phenotype
            in a regression model including the PRS as a
            covariate
            \item \code{prs}: The PRS calculated according
            to \code{prsCalc}
        }
        \item \code{plots}: A list of two ggplot objects:
        \itemize{
            \item \code{scatter}: A scatterplot of
            \code{raw_pheno} against \code{prs}
            \item \code{hist}: A histogram of \code{prs}
        }
        \item \code{permutations}: This is available only
        when \code{times > 0} and only when summary 
        statistics (\code{effects} slot) are present in the
        \code{gwe} GWASExperiment object.
        \itemize{
            \item \code{metrics}: A matrics of evaluation
            metrics (same as the \code{metrics} upper level
            member) over \code{times} iterations (the number
            of rows)
            \item \code{bootp}: A bootstrap p-value denoting
            the statistical significance of the tested PRS
            being non-random (essentially the number of times
            the PRS R^2 of the random PRS is larger than
            the original PRS R^2 divided by \code{times})
        }
    }
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
    
    Regarding the \code{times} argument in the 
    \code{applyPRS} function, this instructs the function 
    to perform \code{times} random SNP selections from the
    input \code{gwe} and repeat the validation procedure 
    through the calculation of evaluation metrics. This 
    process is available only when summary statistics are
    present in \code{gwe} and is useful to test how a random
    PRS performs in explaining the examined phenotype. 
    Normally, the metrics should be worse. A boostrap 
    p-value is also calculated for some of the evaluation 
    metrics.
    
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
    
    can <- getPrsCandidates(prismaOut,"glm",1)
    m <- prsCrossValidate(can,toy,response=response,
        covariates=covariates,leaveOut=0.2,times=5)
}
\author{
    Panagiotis Moulos
}
