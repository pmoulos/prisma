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
        psig = 0.05, seed = 42, cliopts = NULL, wspace = NULL, 
        rc = NULL)
    
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
    
    \item{glmOpts}{A list of options for GLM regression 
    model fit. See also Details.}
    
    \item{rrblupOpts}{A list of options for ridge regression 
    model fit (rrBLUP). See also Details.}
    
    \item{statgenOpts}{A list of options for regression fit
    with statgenGWAS. See also Details.}
    
    \item{snptestOpts}{A list of options for modeling
    associations with SNPTEST. See also Details.}
    
    \item{plinkOpts}{A list of options for modeling
    associations with PLINK. See also Details.}
    
    \item{usepc}{What to do with PCA and PCs when using
    \code{gwaBlup} and rrBLUP model fitting. Defaults to
    \code{"auto"} and also \code{"rrint"} is preferred.
    See also Details.}
    
    \item{npcs}{Number of PCs to use with rrBLUP. Defaults
    to \code{NULL} for auto-determining what to do.}
    
    \item{family}{The model family for GLM, e.g. 
    \code{"gaussian"} for classic regression with continuous
    response variable or \code{"binomial"} for logistic
    regression. By default (\code{NULL}), it is guessed from
    the form of the \code{response} variable.}
    
    \item{wspace}{A path to create workspaces and temporary
    files for SNPTEST and PLINK which operate outside R.}
    
    \item{test}{One of supported test types for SNPTEST. It
    can be \code{"frequentist"} (default) or 
    \code{"bayesian"}. See also SNPTEST documentation.}
    
    \item{seed}{A random seed for PLINK resampling 
    operations. Defaults to \code{42}.}
    
    \item{cliopts}{A string with additional command-line
    parameters to pass to the PLINK executable. Defaults
    to \code{NULL}.}
    
    \item{model}{One of the supported genotype modeling
    types for SNPTEST. It can be one of \code{"additive"}
    (default), \code{"dominant"}, \code{"recessive"},
    \code{"general"} or \code{"heterozygote"}. See also
    SNPTEST documentation.}
    
    \item{rc}{The fraction of available cores to use for 
    parallel calculations. Default is \code{NULL} for no
    parallelization.}
}
\value{
    Depending on the function used, either a data frame
    with summary statistics (always with at least SNP
    name, effects/ORs and p-values)
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
    Regarding the \code{methods} argument, the following
    briefly explain how the associations are modeled in
    each case:
    \itemize{
        \item \code{glm}: In this case, a SNP is checked
        for association with the \code{response} trait
        using a General Linear Regression model of the
        form: \code{response ~ covariates + [SNP]} or
        \code{response ~ covariates + pcs + [SNP]}
        depending on whether PCs are added as covariates
        or not. \code{[SNP]} is simply the recorded genotype
        (\code{0} for homozygous minor allele, \code{1}
        for heterozygous allele and \code{2} for homozygous
        major allele). The p-value and the regression
        coefficient (effect/OR) of the [SNP] is recorded
        and used to determine if the SNP is associated
        with the covariate or not. The function can be
        used independently by directly using \code{gwaGlm}.
        In this case, a data frame of summary statistics
        is returned instead of the input 
        \code{GWASExperiment}.
        \item \code{rrblup}: In this case, the 
        \code{\link{GWAS}} function from the package
        rrBLUP is used to model the associations, the
        effects and the significance of each SNP using
        Ridge Regression. The function can be used
        independently by directly using \code{gwaBlup}.
        In this case, a data frame of ridge regression 
        summary statistics is returned instead of the 
        input \code{GWASExperiment}.
        \item \code{statgen}: In this case, the
        \code{\link{runSingleTraitGwas}} function from
        the statgenGWAS package is used to model the 
        associations, the effects and the significance 
        of each SNP using also linear regression. The 
        function can be used independently by directly 
        using \code{gwaStatgen}. In this case, a data 
        frame of summary statistics is returned instead 
        of the input \code{GWASExperiment}.
        \item \code{snptest}: In this case, the SNPTEST
        software is used to check for associations
        and summary statistics. It should be noted that
        SNPTEST operates with at least 100 samples. The
        function can be used directly as the previous 
        ones by \code{gwaSnptest}.
        \item \code{plink}: In this case, the PLINK
        software is used to check for associations and
        return summary statistics. The function can be used
        directly as the previous ones by \code{gwaPlink}.
    }
    
    Regarding the \code{combine} argument, this is used to
    control how the p-values returned by the usage of many
    association tests will be combined to produce a 
    summarized p-value. Currently \code{"whitlock"}, 
    \code{"pandora"} and \code{"harmonic"} are placeholders
    since the estimation of p-value weights are work in
    progress. For the rest, please see Moulos and Hatzis,
    2015 (PMID: 25452340) and Fanidis and Moulos, 2021 
    (PMID: 32778872).
    
    Regarding the \code{glmOpts} argument, this is a named
    list with members:
    \itemize{
        \item \code{family}: The regression family, which
        defaults to \code{NULL} for auto-detection.
        \item \code{size}: The number of sequential SNPs to
        include in a lasso regression model (currently in
        testing mode).
    }
    
    Regarding the \code{rrblupOpts} argument, this is a 
    named list with members:
    \itemize{
        \item \code{pcblup}: How should the PCs included
        in ridge regression be calculated? Defaults to
        \code{"auto"}, see the main function arguments
        for the rest of the options.
        \item \code{npcs}: The number of PCs to include.
    }
    
    Regarding the \code{statgenOpts} argument, this is a 
    named list with members:
    \itemize{
        \item \code{kinship}: The kinship matrix estimation
        method, to be included in the statgenGWAS model
        fitting. Defaults to \code{"astle"}, see the
        statgenGWAS package for other options.
        \item \code{reml}: The method for calculated the
        Restricted Maximum Likelihood in the statgenGWAS
        model fitting. Defaults to \code{"EMMA"}, see the
        statgenGWAS package for other options.
        \item \code{gls}: The method used to estimate the
        SNP effects, \code{"signle"} (default) for using
        a single kinship matrix, or \code{"multi"} for 
        using chromosome specific kinship matrices. See
        the statgenGWAS package for further explanations.
    }
    
    Regarding the \code{snptestOpts} argument, this is a 
    named list with members:
    \itemize{
        \item \code{test}: The type of test to be used (see
        also main arguments above). Defaults to 
        \code{"frequentist"}. See also SNPTEST documentation.
        \item \code{model}: The genotype effect estimation
        mode, defaults to \code{"additive"}, See also main
        arguments and SNPTEST documentation.
        \item \code{workspace}: The path to SNPTEST
        temporary workspace. Defaults to \code{NULL} for
        using \code{\link{tempdir}. It should be provided
        for maximum control.}
    }
    
    Regarding the \code{plinkOpts} argument, this is a
    named list with members:
    \itemize{
        \item \code{effect}: How should the effects be
        modeled for each SNP (similar to the \code{model}
        of SNPTEST). Please see PLINK documentation for a
        full list of options.
        \item \code{seed}: A seed for PLINK randomized 
        operations such as resampling and permuting.
        Defaults to \code{42}.
        \item \code{cliopts}: Further options to be passed
        to \code{plink} command, as these would be given
        to the command line. Please note that these should
        be switched and options accepting numbers or srings
        but _not_ input/output files.
    }
    
    Regarding the \code{usepc} argument of \code{gwaBlup},
    this can be one of \code{"auto"}, \code{"estim"},
    \code{"rrint"}, \code{"fixed"} or \code{"none"}. It
    controls how rrBLUP models use the PCs in the generic
    ridge regression model. The following are available:
    \itemize{
        \item \code{"auto"}: When in automated mode, if
        the input object contains calculated PCs to be
        used as covariates, these are used, otherwise
        the optimal number is estimate (switch to 
        \code{"estim"}.)
        \item \code{"estim"}: Automatic iterative estimation.
        \item \code{"rrint"}: In this case, the number of PCs
        to include are taken from existing PCs in the input
        \code{GWASExperiment} object, or from the \code{npcs}
        argument if there are no PCs in the input. Then,
        this number is used internally to the \code{GWAS}
        rrBLUP function to include PCs in the model. Although
        not the default, it seems that this is the more
        stable method and it is recommended to be used.
        \item \code{"fixed"}: This uses \code{npcs} only
        (in a similar way as \code{"auto"}.)
        \item \code{"none"}: Do not include PCs in the
        model.
    }
}
\examples{
    data(toy,package="prisma")
    
    # Only 1 PC will be returned because of dataset size
    toy <- suppressWarnings(calcPcaCovar(toy,method="snprel"))
    
    toyG <- gwa(toy,response="BMI",covariates=c("Age","Gender"),
        pcs=TRUE,psig=0.05,methods=c("glm","statgen"))
    
    glmSumstat <- gwaGlm(toy,response="BMI",
        covariates=c("Age","Gender"),pcs=TRUE)
    
    rrSumstat <- gwaBlup(toy,response="BMI",
        covariates=c("Age","Gender"),usepc="rrint",npcs=1)
    
    sgSumstat <- gwaStatgen(toy,response="BMI",
        covariates=c("Age","Gender"),pcs=TRUE)
    
    if (Sys.which("snptest") != "")
        snpSumstat <- gwaSnptest(toy,response="BMI",
            covariates=c("Age","Gender"),pcs=TRUE)
    
    if (Sys.which("plink") != "")
        plSumstat <- gwaPlink(toy,response="BMI",
            covariates=c("Age","Gender"),pcs=TRUE)
}
\author{
    Panagiotis Moulos
}
