\name{prismaFuns}
\alias{prisma}
\alias{prsSelection}
\alias{prsPipeline}
\title{Main prisma PRS extraction and validation pipelines}
\usage{
    prisma(
        gwe,
        phenotype,
        covariates,
        pcs = FALSE,
        npcs = 0,
        trainSize = 0.8,
        niter = 10,
        resolution = c("frequency", "quantile"),
        step = if (resolution == "frequency") 1 else 
            c(0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 0.8, 0.9, 0.95, 0.99),
        minFreq = 2,
        minSnps = 5,
        dropSameQuantiles = TRUE,
        aggregation = c("intersection", "union"),
        effectWeight = c("mean", "median", "weight"),
        prsSelectMethod = c("maxima", "elbow"),
        prsSelectCrit = c("prs_r2", "prs_pvalue", "prs_aic"),
        prsSelectStat = c("mean", "median", "none"),
        prsSelectR2 = c("adjusted", "raw"),
        sigTest = c("ttest", "wilcoxon", "empirical"),
        filters = getDefaults("filters"),
        pcaMethod = c("auto", "snprel", "grid", "hubert"),
        imputeMissing = FALSE,
        imputeMethod = c("single", "split"),
        gwaMethods = c("glm", "rrblup", "statgen", "snptest",
            "plink"),
        gwaCombine = c("fisher", "simes", "max", "min", "harmonic", 
            "whitlock", "pandora"),
        glmOpts = getDefaults("glm"),
        rrblupOpts = getDefaults("rrblup"),
        statgenOpts = getDefaults("statgen"),
        snptestOpts = getDefaults("snptest"),
        plinkOpts = getDefaults("plink"),
        prsMethods = c("lassosum","prsice"),
        lassosumOpts = getDefaults("lassosum"),
        prsiceOpts = getDefaults("prsice"),
        prsWorkspace = NULL,
        cleanup = c("none", "intermediate", "all"),
        logging = c("screen", "file"),
        dnOutput = c("gwaslist", "summaries"),
        output = c("normal", "full"),
        continue = FALSE,
        useDenovoWorkspace = NULL,
        runId = NULL,
        evalOnSplit = c("original", "new"),
        evalWith = c("vanilla", "prscice"),
        rc = NULL
    )
    
    prsPipeline(
        gwe,
        phenotype,
        covariates,
        pcs = FALSE,
        npcs = 0,
        snpSelection = NULL,
        trainSize = 0.8,
        niter = 10,
        filters = getDefaults("filters"),
        pcaMethod = c("auto", "snprel", "grid", "hubert"),
        imputeMissing = FALSE,
        imputeMethod = c("single", "split"),
        gwaMethods = c("glm", "rrblup", "statgen", "snptest",
            "plink"),
        gwaCombine = c("fisher", "simes", "max", "min", "harmonic", 
            "whitlock", "pandora"),
        glmOpts = getDefaults("glm"),
        rrblupOpts = getDefaults("rrblup"),
        statgenOpts = getDefaults("statgen"),
        snptestOpts = getDefaults("snptest"),
        plinkOpts = getDefaults("plink"),
        prsMethods = c("lassosum","prsice"),
        lassosumOpts = getDefaults("lassosum"),
        prsiceOpts = getDefaults("prsice"),
        prsWorkspace = NULL,
        cleanup = c("none", "intermediate", "all"),
        logging = c("screen", "file"),
        output = c("gwaslist", "summaries"),
        continue = FALSE,
        useDenovoWorkspace = NULL,
        runId = NULL,
        evalWith = c("vanilla", "prscice"),
        rc = NULL
    )
    
    prsSelection(
        dnList,
        gwe,
        phenotype,
        covariates,
        pcs = FALSE,
        npcs = 0,
        trainSize = 0.8,
        niter = 10,
        resolution = c("frequency", "quantile"),
        step = if (resolution == "frequency") 1 else 
            c(0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 0.8, 0.9, 0.95, 0.99),
        minFreq = 2,
        minSnps = 5,
        dropSameQuantiles = TRUE,
        aggregation = c("intersection", "union"),
        effectWeight = c("mean", "median", "weight"),
        filters = getDefaults("filters"),
        pcaMethod = c("auto", "snprel", "grid", "hubert"),
        imputeMissing = FALSE,
        imputeMethod = c("single", "split"),
        gwaMethods = c("glm", "rrblup", "statgen", "snptest",
            "plink"), # lasso later
        gwaCombine = c("fisher", "simes", "max", "min", "harmonic", 
            "whitlock", "pandora"),
        glmOpts = getDefaults("glm"),
        rrblupOpts = getDefaults("rrblup"),
        statgenOpts = getDefaults("statgen"),
        snptestOpts = getDefaults("snptest"),
        plinkOpts = getDefaults("plink"),
        prsMethods = c("lassosum","prsice"),
        lassosumOpts = getDefaults("lassosum"),
        prsiceOpts = getDefaults("prsice"),
        prsWorkspace = NULL,
        cleanup = c("none", "intermediate", "all"),
        logging = c("screen", "file"),
        continue = FALSE,
        useDenovoWorkspace = NULL,
        runId = NULL,
        evalWith = c("vanilla", "prscice"),
        output = c("summary", "full"),
        rc = NULL
    )
}
\arguments{
    \item{gwe}{A GWASExperiment object with GWA tests performed.}
    
    \item{dnList}{A list resulting from a \emph{de novo} PRS 
    extraction pipeline. See also Details.}
    
    \item{phenotype}{}
    
    \item{covariates}{}
    
    \item{pcs}{}
    
    \item{npcs}{}
    
    \item{trainSize}{}
    
    \item{niter}{}
    
    \item{resolution}{}
    
    \item{step}{}
    
    \item{minFreq}{}
    
    \item{minSnps}{}
    
    \item{dropSameQuantiles}{}
    
    \item{aggregation}{}
    
    \item{effectWeight}{}
    
    \item{prsSelectMethod}{}
    
    \item{prsSelectCrit}{}
    
    \item{prsSelectStat}{}
    
    \item{prsSelectR2}{}
    
    \item{sigTest}{}
    
    \item{filters}{}
    
    \item{pcaMethod}{}
    
    \item{imputeMissing}{}
    
    \item{imputeMethod}{}
    
    \item{gwaMethods}{}
    
    \item{gwaCombine}{}
    
    \item{glmOpts}{}
    
    \item{rrblupOpts}{}
    
    \item{statgenOpts}{}
    
    \item{snptestOpts}{}
    
    \item{plinkOpts}{}
    
    \item{prsMethods}{}
    
    \item{lassosumOpts}{}
    
    \item{prsiceOpts}{}
    
    \item{prsWorkspace}{}
    
    \item{cleanup}{}
    
    \item{logging}{}
    
    \item{dnOutput}{}
    
    \item{output}{}
    
    \item{continue}{}
    
    \item{useDenovoWorkspace}{}
    
    \item{runId}{}
    
    \item{evalWith}{}
    
    \item{rc}{The fraction of available cores to use for 
    parallel calculations. Default is \code{NULL} for no
    parallelization.}
}
\value{
    
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

