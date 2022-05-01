\name{cvSumFuns}
\alias{cvSummary}
\alias{plotCvMetrics}
\alias{summarizeCvMetrics}
\title{Functions to summarize PRS cross-validations}
\usage{
    summarizeCvMetrics(cvm)
    
    cvSummary(M, nsnp, scr = TRUE, out = FALSE)
    
    plotCvMetrics(cvo, what = c("r2", "rmse", "mae", "pr2","crl"),
        pl = c("pct", "prs"), silent = FALSE)
}
\arguments{
    \item{cvm}{the output of \code{link{prismaCrossValidate}}
    for one GWA method. See also Details.}
    
    \item{M}{a list member of \code{cvm} above for one
    PRS candidate, for one leave-out fraction. Contains
    metrics at this particular leave-out percentage of 
    total samples. See examples.}
    
    \item{nsnp}{number of SNPs in PRS candidate. If not
    provided, guessed from \code{M}.}
    
    \item{scr}{print the output of \code{cvSummary} to
    screen? See also Details.}
    
    \item{out}{store the output of \code{cvSummary} to
    a variable? See also Details.}
    
    \item{cvo}{same as \code{cvm}.}
    
    \item{what}{which metric to plot out of the available
    ones? See also Details.}
    
    \item{pl}{\code{"pct"} for plots according to % of 
    leave-out samples for all PRS candidates, \code{"prs"} 
    for plots according to numbr of SNPs in PRS for a given
    percentage. Essentially, the two types of CV plots in
    a PRISMA report.}
    
    \item{silent}{Do the plotting without messages? Defaults
    to \code{FALSE}.}
}
\value{
    For \code{summarizeCvMetrics}, a list of length the same
    as the PRS candidates for a specific GWA method. Each
    list member contains the metrics collected by
    \code{\link{prismaCrossValidate}} across the requested
    number of iterations. See also the man page of
    \code{\link{prismaCrossValidate}}.
    
    For \code{cvSummary}, the summary of CV metrics for a 
    PRS candidate across the requested iterations in
    \code{\link{prismaCrossValidate}}. See Details for
    a detailed explanation of what is calculated.
    
    For \code{plotCvMetrics}, a list of the requested
    ggplot2 objects, named according to the type of
    plots (leave-out percentage based, or PRS size
    based). See also examples.
}
\description{
    These functions facilitate the summarization of lenghty
    PRS candidate cross-validation runs using multiple
    linear regression samples agains the main phenotype under
    investigation. They are not mean to be used directly
    (although not prohibited if you know what to expect and
    how to interpret) but they are rather used to generate
    plots and numbers for the PRISMA final HTML report.
}
\details{
    Reagrding the \code{cvm} argument, this is a list of
    lists, containing several cross-validation metrics. The
    main members are named the number of SNPs in a PRS
    candidate and the sub-members are named after the 
    leave-out splits of the cross validation. These 
    structures are not meant to be used directly, although 
    provided. They are to feed the overall PRISMA HTML 
    report. Typically, this is a member of the list output
    produced by the function \code{\link{prismaCrossValidate}}
    which has as many \code{cvm}'s as the number of GWA
    methods that were used in the main \code{\link{prisma}}
    run.
    
    Regarding the \code{scr} and \code{out} variables, these
    control how should the summarized cross-validation 
    metrics should be displayed and/or stored. If \code{out}
    is \code{TRUE}, the output is stored to variable for 
    later use e.g. from \code{\link{prismaReport}} while if
    \code{scr} is \code{TRUE} the output is nicely printed
    on screen. If both are \code{TRUE}, then also a varable
    is returned. In this case be careful of the syntax as
    the screen may be flooded with information! Both cannot
    be set to \code{FALSE}. You will need some output!
    
    Regarding the \code{what} argument, the following 
    explain the different options:
    \itemize{
        \item \code{r2}: model R^2 values (train samples)
        \item \code{rmse}: Root Mean Squared Error
        \item \code{mae}: Mean Absolute Error
        \item \code{pr2}: predicted R^2 values (LOO samples)
        \item \code{crl}: correlations with phenotype
    }
    
    Regarding the output of \code{cvSummary}, the following
    are calculated for a PRS candidate, for a specific
    percentage of leave-out samples and for a specific
    number of cross-validation iterations:
    \itemize{
        \item General information
        \itemize{
            \item Number of SNPs in PRS (\code{n_snp}): The
            number of SNPs in the PRS candidate being 
            evaluated.
            \item Number of cross-validation (\code{cv_times}):
            The number of leave-out cross-validations. In each
            iteration, a random number of samples comprising
            the requested percentage of leave-out samples of
            the dataset is left out and the PRS validated
            based on regressing the phenotypes of the rest
            samples including the PRS as a covariate and
            evaluating the model based on the left-out 
            samples.
            \item Statistical significance of the PRS in
            model (\code{max_prs_pvalue}): The best p-value
            of the PRS as a covariate in the CV regression
            model.
        }
        \item R^2 summary statistics
        \itemize{
            \item R^2 of the full regression model including 
            the PRS  (\code{mean_full_r2}, \code{sd_full_r2}): 
            The average R^2 of the full model including the 
            PRS +/- standard deviation.
            \item R^2 of the reduced regression model 
            excluding the PRS (\code{mean_reduced_r2}, 
            \code{sd_reduced_r2}): The average R^2 of the 
            reduced regression model, including all other
            covariates but excluding the PRS, +/- standard 
            deviation.
            \item Adjusted R^2 of the PRS full model 
            contribution (\code{mean_prs_r2}, 
            \code{sd_prs_r2}): The average adjusted R^2
            of the full model, denoting the contribution
            of the PRS (essentially R^2 of the full model
            minus the R^2 of the reduced model), +/-
            standard deviation.
        }
        \item Statistical significance of full R^2 against 
        reduced R^2
        \itemize{
            \item p-value of a Student's t-test 
            (\code{ttest_full_reduced_r2}) for the hypothesis
            of the mean R^2 of the full model being larger
            than the mean R^2 of the reduced model.
            \item p-value of a Wiloxon rank sum test
            (\code{wilcox_full_reduced_r2}) for the 
            hypothesis of the median R^2 of the full model 
            being larger than the median R^2 of the reduced 
            model.
            \item Empirical p-value 
            (\code{empirical_full_reduced_r2}) denoting the
            number of times the R^2 of the full model was 
            smaller than the R^2 of the reduced model across
            the number of cross-validation iterations.
        }
        \item RMSE (Root Mean Squared Error) summary 
        statistics
        \itemize{
            \item RMSE of the full regression model 
            including the PRS (\code{mean_full_rmse}, 
            \code{sd_full_rmse}): Average RMSE of the 
            full model with the PRS +/- standard deviation.
            \item RMSE of the reduced regression model 
            excluding the PRS (\code{mean_reduced_rmse}, 
            \code{sd_reduced_rmse}): Average RMSE of the 
            reduced model including all other covariates 
            but excluding the PRS, +/- standard deviation.
        }
        \item Statistical significance of full RMSE against 
        reduced RMSE
        \itemize{
            \item p-value of a Student's t-test 
            (\code{ttest_full_reduced_rmse}) for the 
            hypothesis of the mean RMSE of the full model 
            being larger than the mean RMSE of the reduced
            model.
            \item p-value of a Wiloxon rank sum test
            (\code{wilcox_full_reduced_rmse}) for the 
            hypothesis of the median RMSE of the full model 
            being larger than the median RMSE of the reduced 
            model.
            \item Empirical p-value 
            (\code{empirical_full_reduced_rmse}) denoting the
            number of times the RMSE of the full model was 
            smaller than the RMSE of the reduced model across
            the number of cross-validation iterations.
        }
        \item MAE (Mean Absolute Error) summary statistics
        \itemize{
            \item MAE of the full regression model 
            including the PRS (\code{mean_full_mae}, 
            \code{sd_full_mae}): Average MAE of the 
            full model with the PRS +/- standard deviation.
            \item MAE of the reduced regression model 
            excluding the PRS (\code{mean_reduced_mae}, 
            \code{sd_reduced_mae}): Average MAE of the 
            reduced model including all other covariates 
            but excluding the PRS, +/- standard deviation.
        }
        \item Statistical significance of full MAE against 
        reduced MAE
        \itemize{
            \item p-value of a Student's t-test 
            (\code{ttest_full_reduced_mae}) for the 
            hypothesis of the mean MAE of the full model 
            being larger than the mean MAE of the reduced
            model.
            \item p-value of a Wiloxon rank sum test
            (\code{wilcox_full_reduced_mae}) for the 
            hypothesis of the median MAE of the full model 
            being larger than the median MAE of the reduced 
            model.
            \item Empirical p-value 
            (\code{empirical_full_reduced_mae}) denoting the
            number of times the MAE of the full model was 
            smaller than the MAE of the reduced model across
            the number of cross-validation iterations.
        }
        \item Correlation summary statistics
        \itemize{
            \item R between observed & predicted values with 
            the full model (\code{mean_full_pred_cor}, 
            \code{sd_full_pred_cor}): The average correlation
            between the actual observed phenotype values and
            the predicted by the full regression model,
            including the PRS, +/- standard deviation.
            \item R between observed & predicted values with 
            the reduced model (\code{mean_reduced_pred_cor}, 
            \code{sd_reduced_pred_cor}): The average 
            correlation between the actual observed phenotype 
            values and the predicted by the reduced regression 
            model, including all the other covariates but 
            excluding the PRS, +/- standard deviation.
        }
        \item Statistical significance of full R against 
        reduced R
        \itemize{
            \item p-value of a Student's t-test 
            (\code{ttest_full_reduced_pred_cor}) for the 
            hypothesis of the mean R of the full model 
            being larger than the mean R of the reduced
            model.
            \item p-value of a Wiloxon rank sum test
            (\code{wilcox_full_reduced_pred_cor}) for the 
            hypothesis of the median R of the full model 
            being larger than the median R of the reduced 
            model.
            \item Empirical p-value 
            (\code{empirical_full_reduced_pred_cor}) denoting 
            the number of times the R of the full model was 
            smaller than the R of the reduced model across
            the number of cross-validation iterations.
        }
        \item R^2 summary statistics for the predicted 
        test values
        \itemize{
            \item R^2 of the full regression model including 
            the PRS  (\code{mean_full_pred_r2}, 
            \code{sd_full_pred_r2}): The average R^2 of the 
            full model including the PRS +/- standard 
            deviation, for the predicted leave-out test 
            values.
            \item R^2 from the reduced regression model 
            excluding the PRS (\code{mean_reduced_pred_r2}, 
            \code{sd_reduced_pred_r2}): The average R^2 of 
            the reduced regression model, including all other
            covariates but excluding the PRS, +/- standard 
            deviation, for the predicted leave-out test 
            values.
            \item Adjusted R^2 from the PRS full model 
            contribution (\code{mean_prs_pred_r2}, 
            \code{sd_prs_pred_r2}): The average adjusted R^2
            of the full model, denoting the contribution
            of the PRS (essentially R^2 of the full model
            minus the R^2 of the reduced model), +/-
            standard deviation, for the predicted leave-out 
            test values.
        }
        \item Significance of full R^2 against reduced R^2 
        (predicted)
        \itemize{
            \item p-value of a Student's t-test 
            (\code{ttest_full_reduced_pred_r2}) for the 
            hypothesis of the mean R^2 of the full model being 
            larger than the mean R^2 of the reduced model, for 
            the predicted test values.
            \item p-value of a Wiloxon rank sum test
            (\code{wilcox_full_reduced_pred_r2}) for the 
            hypothesis of the median R^2 of the full model 
            being larger than the median R^2 of the reduced 
            model, for the predicted test values.
            \item Empirical p-value 
            (\code{empirical_full_reduced_pred_r2}) denoting 
            the number of times the R^2 of the full model was 
            smaller than the R^2 of the reduced model across
            the number of cross-validation iterations, for
            the predicted test values.
        }
    }
}
\examples{
    data(toy,package="prisma")
    data(prisma_out,package="prisma")
    
    response <- "BMI"
    covariates <- c("Case_Control","Gender","Age")
    
    cvm <- prismaCrossValidate(prismaOut,toy,response=response,
        covariates=covariates,leaveOut=c(0.2,0.4),times=5)
    s <- summarizeCvMetrics(cvm[["glm"]])
    identical(names(cvm[["glm"]]),names(s)) # TRUE
    identical(nrow(s[[1]]),2L) # TRUE
    # Has two leave-out sample percentages 20%, 40%
    
    # Let's display summarized metrics for the 24-SNP PRS 
    # derived by GLM when leaving out 20% of the samples
    # of the initial dataset
    M <- cvm[["glm"]][["24"]][["0.2"]]
    cvSummary(M,nsnp=24,scr=TRUE)
    
    # Some metric plots (ggplot2 graphs)
    cvo <- cvm[[1]]
    p1 <- plotCvMetrics(cvo,what="r2",pl="prs")
    # p1
    p2 <- plotCvMetrics(cvo,what="rmse",pl="pct")
    # p2
}
\author{
    Panagiotis Moulos
}
