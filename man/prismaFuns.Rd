\name{prismaFuns}
\alias{prisma}
\alias{prsSelection}
\alias{prsPipeline}
\alias{prismaPipeline}
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
        sitadelaDb = NULL,
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
    
    prismaPipeline(..., retry = 10, outPath = NULL)
}
\arguments{
    \item{gwe}{A GWASExperiment object with GWA tests performed.
    In other words, \code{\link{gwa}} must have been applied.}
    
    \item{dnList}{A list resulting from a \emph{de novo} PRS 
    extraction pipeline. See also Details.}
    
    \item{phenotype}{The trait for which a PRS is sought to be
    extracted. It should be a character of length 1 and exist
    in the phenotypes of \code{gwe}.}
    
    \item{covariates}{The covariates to include in regression
    models. It should be a character vector of traits all 
    existing in the phenotypes of \code{gwe}.}
    
    \item{pcs}{Inlcude PCs in the regression model? Defaults
    to \code{FALSE}.}
    
    \item{npcs}{The number of PCs to use if \code{pcs=TRUE}.
    If \code{npcs=0}, the PCs available in \code{gwe} will
    be used.}
    
    \item{trainSize}{The fraction of the training set used
    for summary statistics calculation in each iteration. 
    See also Details. Defaults to \code{0.8}.}
    
    \item{niter}{Number of PRS extraction iterations. 
    Defaults to \code{100}.}
    
    \item{resolution}{How to resolve the best number of
    SNPs in a PRS candidate, either by validating a PRS
    with as many SNPs as the number of appearances at each
    frequency cutoff (\code{frequency}, default) or as many
    SNPs as those included at different quantile cutoffs of 
    the frequency of appearance distribution 
    (\code{"quantile"}). See also Details.}
    
    \item{step}{The number of steps or step increment for
    \code{resolution}. Defaults to \code{1} for 
    \code{resolution="frequency"} or a set of predefined
    quantiles when \code{resolution="quantile"}. See also
    Details.}
    
    \item{minFreq}{The minimum accepted frequency of 
    appearance of a SNP in \code{niter} iterations to be 
    considered as member of a PRS candidate. Defaults to 
    \code{2}.}
    
    \item{minSnps}{The minimum number of SNPs included in
    a PRS candidate. Defaults to \code{5}.}
    
    \item{dropSameQuantiles}{When 
    \code{resolution="quantile"} drop quantiles with the 
    same SNPs? Defaults to \code{TRUE}. See also Details.}
    
    \item{aggregation}{How should the SNPs from different 
    PRS algorithms should be aggregated at the end of the 
    iterations in order to create the distribution of 
    appearances. One of \code{"intersection"} (default)
    or \code{"union"}. See also Details.}
    
    \item{effectWeight}{How to summarize SNP effects 
    (regression coefficients) across the \code{niter}
    iterations. One of \code{"mean"} (default),
    \code{"mean"} or \code{"weight"}. See also Details.}
    
    \item{prsSelectMethod}{How to select the number of SNPs
    for the PRS after \code{niter} iterations. One of
    \code{"maxima"} (default) or \code{"elbow"}. See also
    Details.}
    
    \item{prsSelectCrit}{One of \code{"prs_r2"} (default), 
    \code{"prs_pvalue"} or \code{"prs_aic"}. Determines what
    criterion is used for PRS selection. See also Details.}
    
    \item{prsSelectStat}{One of \code{"mean"} (default),
    \code{"median"} or \code{"none"}. Determines which
    calculated statistic (mean, median) is used to summarize
    the \code{prsSelectCrit} in the evaluation metrics and
    when applicable (e.g. not applicable in p-values).}
    
    \item{prsSelectR2}{One of \code{"adjusted"} (default) or 
    \code{"raw"}. The first one uses the PRS contribution to
    the evaluation regression model and the second one the 
    R^2 of the full model.}
    
    \item{sigTest}{How to assess statistical significance
    between initial PRS R^2 calculations and the respective
    R^2 (or other \code{prsSelectCrit} metrics) from the 
    prisma aggreagated result. One or of \code{"ttest"} 
    (default), \code{"wilcoxon"} or \code{"empirical"}. See
    also Details.}
    
    \item{filters}{Filters to be applied to the \code{gwe},
    during the workflow executed in each prisma iteration.
    For details, see the \code{\link{filterGWAS}} function.}
    
    \item{pcaMethod}{Which method to use for PCA. It can be
    one of \code{"auto"} (default), \code{"snprel"} for the
    parallel PCA method implemented in the package SNPRelate, 
    \code{"grid"} or \code{"hubert"} for robust (but very 
    slow) PCA. See also the \code{\link{gwa}} function.}
    
    \item{imputeMissing}{Impute missing values in each
    prisma iteration? Defaults to \code{TRUE}, see also
    Details.}
    
    \item{imputeMethod}{One of \code{"single"} (default) or
    \code{"split"}. See the \code{\link{imputeGWAS}} for
    details.}
    
    \item{gwaMethods}{The GWAS SNP effect estimation 
    method(s) to use in order to estimate the effects that
    will be used for PRS weights. It can be one or more from
    \code{"glm"}, \code{"rrblup"}, \code{"snptest"},
    \code{"statgen"} and \code{"plink"}. See the 
    \code{\link{gwa}} function for more details.}
    
    \item{gwaCombine}{How to combine SNP p-values from more
    than one GWAS effect estimation methods. Currently not
    used.}
    
    \item{glmOpts}{Options list specific to GLM effect 
    estimation. See the \code{\link{gwa}} function for 
    more details.}
    
    \item{rrblupOpts}{Options list specific to rrBLUP 
    effect estimation. See the \code{\link{gwa}} function 
    for more details.}
    
    \item{statgenOpts}{Options list specific to statgenGWAS
    effect estimation. See the \code{\link{gwa}} function 
    for more details.}
    
    \item{snptestOpts}{Options list specific to SNPTEST
    effect estimation. See the \code{\link{gwa}} function 
    for more details.}
    
    \item{plinkOpts}{Options list specific to PLINK
    effect estimation. See the \code{\link{gwa}} function 
    for more details.}
    
    \item{prsMethods}{One or more PRS extraction algorithms
    to be included in the prisma workflow. Can be 
    \code{"prcise"}, \code{"lassosum"} or both.}
    
    \item{lassosumOpts}{Options list specific to lassosum
    PRS extraction. See the \code{\link{runPRS}} function
    for more details.}
    
    \item{prsiceOpts}{Options list specific to PRSice2
    PRS extraction. See the \code{\link{runPRS}} function
    for more details.}
    
    \item{prsWorkspace}{A local directory to be used for
    prisma operations. Defaults to \code{NULL} for
    auto-generation. See also Details.}
    
    \item{cleanup}{Cleanup level after the finish of the
    prisma pipeline. \code{"none"} (default) leaves all
    intermediate output intact. \code{"intermediate"} 
    deletes the files produced in effect estimation as
    well as BED, BIM and FAM files required by PRSice2
    and lassosum linkage files. \code{"all"} deletes
    the entire \code{prsWorkspace}.}
    
    \item{logging}{One of \code{"screen"} (default) or
    \code{"file"}. Where to write log messages. If prisma
    runs in parallel, it is \code{"file"} all the time.}
    
    \item{dnOutput}{The output of \code{prsPipeline} called
    from this function. Can be \code{"gwaslist"} (default) 
    for a list of GWASExperiments or \code{"summaries"} for
    a summary required for next steps of the prisma pipelie.
    The latter has obviously a smaller size.}
    
    \item{output}{The level of output of this function. It
    can be \code{"normal"} for a list with the necessary
    information to create the prisma report or \code{"full"}
    which also contains the itermediate output of 
    \code{prsPipeline} called within prisma.}
    
    \item{continue}{Continue a crashed pipeline? Defaults
    to \code{FALSE}, see Details.}
    
    \item{useDenovoWorkspace}{Validate the prisma combined
    PRS in the initial dataset splits or generate new ones.
    Defaults to \code{NULL}, see also Details.}
    
    \item{runId}{A unique run id for the run. If not given,
    it will be auto-generated.}
    
    \item{evalOnSplit}{Which data split to use to evaluate the
    PRS. \code{"original"} (default) uses the original split 
    and SNP effect calculations, while \code{"new"} creates
    a new split taking much more time.}
    
    \item{evalWith}{How to perform the PRS evaluation 
    regressions. It can be \code{"vanilla"} (default) to do
    it with R or \code{"prsice"} to do it with PRSice2. The
    latter takes much more time as there is intermediate
    writing of PLINK files.}
    
    \item{sitadelaDb}{Either \code{NULL} or a sitadela SQLite 
    database with at least human RefSeq annotations for hg19
    and hg38. It is used to add gene names to PRS SNP loci.}
    
    \item{snpSelection}{A data frame with candidate SNPs, the
    output, or part of the output of 
    \code{\link{aggregatePrsMarkers}}, or one of the list 
    members of the output of \code{\link{selectPrs}}.}
    
    \item{rc}{The fraction of available cores to use for 
    parallel calculations. Default is \code{NULL} for no
    parallelization.}
    
    \item{...}{Arguments passed to \code{prisma} from the
    \code{prismaPipeline} wrapper..}
    
    \item{retry}{The number of pipeline relaunches in case
    of crash for any reason. Defaults to \code{10}.}
    
    \item{outPath}{The path to the directory where 
    checkpoint files and the report will be written. It
    can be different from \code{prsWorkspace} and
    defaults to \code{NULL} for auto-generation.}
}
\value{
    A list with the members \code{params} and \code{results}.
    Both are lists.
    
    The \code{params} member has the following contents:
    \itemize{
        \item \code{call} which contains the function call.
        \item \code{args} which contains a list of all
        arguments used, default and user-specific, 
        irrespective of \code{call}.
    }
    
    The \code{results} member is a named list with length as
    the vector \code{gweMethods} and each sublist corresponds
    to results for each SNP effect calculation method. Each 
    one of them has the following contents:
    \itemize{
        \item \code{candidates}, a named list of data frames
        containing the SNPs within each PRS candidate. The
        names correspond to the number of SNPs in the PRS
        candidates and the candidates are determined based
        on \code{prsSelectMethod}.
        \item \code{iterations} which is \code{NULL} or
        the output of each prisma iteration (multiple
        outputs of \code{prsSelection}) if 
        \code{output="full"}.
        \item \code{reportData} which is an internal 
        structure containing PRS evaluation metrics and
        ggplot2 objects used to make the prisma report.
    }
    
    This object can be saved in and .rda file and later used
    with cross-validation metrics and lookup in the GWAS 
    catalog (\code{\link{prismaLookup}} function). These 
    three elements are required to build a prisma report
    (\code{\link{prismaReport}} function).
}
\description{
    These functions comprises the main PRISMA workhorse as 
    they perform all the steps required to extract a PRS 
    potentially chatacterizing a phenotypic trait, given a 
    set of samples, their genetic polymorphims and the 
    phenotypes to participate in the subsequent GWAS 
    associations in a GWASExperiment object. Out of these,
    the \code{prisma} function combines all three and
    performs:
    \itemize{
        \item Input checking and subsequent management
        according to the input details
        \item Iterative calls to \code{prsPipeline} for
        each SNP effect estimation method and is the main 
        prisma iterator for extracting multiple PRSs to be
        later combined.
        \item Call to \code{prsSelection} which evaluates
        the output of \code{prsPipeline} and selects the
        best PRS candidates for each SNP effect estimation
        method.
        \item Wrap up of the data to feed the later PRISMA
        report and basic evaluation graphics.
        \item Output of an object containing all the results 
        and input to the PRISMA reporting function.
    }
    
    \code{prsPipeline} is the iterator of the PRISMA 
    workflow. It will produce a list of GWASExperiment
    objects or their summaries, up to the filling of the
    \code{prsbetas} slot of the GWASExperiment.
    
    \code{prsSelection} accepts as input the output of
    \code{prsPipeline} and evaluate and select the best
    PRS candidates according to the input arguments.
    
    \code{prismaPipeline} seves as a wrapper of 
    \code{prisma} and makes sure that the pipeline will
    resume if crashed at any point as well as write
    checkpoints and create the report.
}
\details{
    Regarding the \code{dnList} argument, this is a list
    of either \code{GWASExperiment} objects or a list of
    lists, each summarizing a \code{GWASExperiment} object
    to save space/memory as not all the elements of a
    \code{GWASExperiment} object are required for downstream
    analyses. The content of \code{dnList} is determined by
    the \code{output} argument of \code{prsPipeline}. When
    the latter is \code{"summaries"}, \code{dnList} sublists
    have the following members:
    \itemize{
        \item baseIndex: the sample indices that were used
        for the extraction of summary statistics via
        \code{gwa} (training).
        \item targetIndex: the sample indices that were
        used for PRS extraction and validation in the
        supported pipelines (testing).
        \item effects: the effects (regression coefficients
        or odds ratios) for each SNP (summary statistics).
        \item prsbetas: a subset of regression coefficients
        and/or adjusted ones in the case of lassosum, which
        comprise the PRS weights.
        \item npcs: the number of principal components 
        included in the regression model.
        \item pr2: a data frame with metrics returned from
        PRSice2, such as PRS p-value and R^2.
    }
    \code{dnList} can essentially be the outcome only of a
    \emph{de novo} PRS extraction pipeline.
    
    Regarding the \code{trainSize} argument, this controls
    the fraction of the input GWASExperiment object (i.e.
    the fraction of the input dataset) that will be used
    for the calculation of summary statistics using the
    \code{link{gwa}} function (usually). These summary
    statistics are then used with the rest of the samples
    (fraction 1 - \code{trainSize}) for PRS extraction,
    that is a input to the supported PRS algorithms 
    (PRSice2, lassosum). Therefore, in each prisma
    iteration (out of \code{niter}), 
    \code{trainSize*ncol(gwe)} samples are used to calculate
    summary statistics and \code{(1-trainSize)*ncol(gwe)}
    are used as targets in the PRS extrction process.
    
    Regarding the \code{dropSameQuantiles} argument, when
    \code{resolution="quantile"}, it is possible that the
    number of SNPs at different quantiles of the frequency
    of appearance distribution is the same. By default,
    the PRS is not evaluated at the highest quantile with
    the same number of SNPs as the next lower one
    (\code{dropSameQuantiles=FALSE}). For example, if
    the number of PRS SNPs at Q=0.5 and Q=0.75 is 1000,
    then there will be no evaluation at Q=0.75. To
    force such an evaluation, set 
    \code{dropSameQuantiles=TRUE}.
    
    Regarding the \code{"aggregation"} argument, this
    controls how should PRS SNP candidates from different
    PRS algorithms should be consolidated. When
    \code{"intersection"}, SNPs included in the PRS
    suggested by all algorithms synthesize the PRS
    candidate set. When \code{"union"}, it is the union
    of SNPs present in all PRSs. The final SNP set
    is then used to create a distribution of appearances
    (frequency of each SNP across \code{"niter"}
    iterations) which is the basis for later selection.
    
    Regarding the \code{effectWeight} argument, this
    controls how the SNP effects in PRS candidates are
    summarized over \code{niter} iterations to produce the
    final effects for the calculation of weighted PRS. When
    \code{"mean"} (default), the final effect for a SNP is
    the mean off all its effects across \code{niter}
    iterations, when it appears. When \code{"median"}, it's
    the respective median while when \code{"weight"},
    the SNP effects are weighted by 1 minus the recorded
    PRS R^2 so that the "best" runs can be given more
    value.
    
    Regarding the \code{prsSelectMethod} argument, this
    controls how the final PRS is selected after the end
    of \code{niter} iterations and based on the requested
    evaluation metrics (typically the PRS R^2). When
    \code{"maxima"}(default), a distribution of PRS R^2
    values is created with the PRS R^2 values (or rarer, 
    other chosen metric) derived from \code{niter} 
    iterations. This distribution (if R^2-related) may
    further be adjusted according to the number of SNPs
    comrpising in the PRS, according to the 
    \code{prsSelectR2} argument. Then, from this 
    distribution, local maxima are calculated. The local
    maximum with the highest R^2 is considered the best
    candidate and provides the number of SNPs. The other 
    maxima are returned also as further candidates. When
    \code{"elbow"} the, elbow point method is used for
    determining the optimal number of SNPs. The elbow method
    is a heuristic used in determining the number of 
    features in a data set with respect to some optimization
    criterion. In this case the criterion is the value of
    PRS R^2 and the method consists of plotting the PRS
    explained variation (as depicted by R^2) as a function 
    of the number of SNPs in the PRS.
    
    Regarding the \code{prsSelectCrit} argument, this
    controls the criterion used to select the final PRS
    after the end of the iterative process and according
    to the evaluation metrics. When \code{"prs_r2"} 
    (default), the decision making plots and the 
    \code{prsSelectMethod} operate on the PRS R^2 calculated
    in each PRS extraction iteration. When 
    \code{"prs_pvalue"}, the -log_10(PRS p-valye) is used
    and when \code{"prs_aic"}, the AIC (Akaike's Information
    Criterion) is used.
    
    Regarding the \code{sigTest} argument, this controls
    how (any) statistical significance is assessed in the
    comparison of the initial PRS R^2 calculations (as this
    is returned by PRSice2) and the PRS R^2 of the prisma
    combined approach. Essentially, the two PRS R^2 (or
    full model R^2 as determined by \code{prsSelectCrit})
    distributions are compared using a t-test (default),
    a Wilcoxon rank test or an empirical test. In the latter
    case, the 'p-value' is the number of times that the
    R^2 of the prisma combination method is smaller than
    the initial R^2 divided by the number of iterations.
    The outcome is presented in the prisma report.
    
    Regarding the \code{imputeMissing} argument, it is
    practically not used with a pre-processed dataset which
    should not contain missing values. If the initial 
    dataset has missing values, imputation will be performed
    in each step, slowing down the procedure.
    
    Regarding the \code{prsWorkspace} argument, this
    controls where the prisma intermediate files and
    results as well as configuration files will be stored.
    If \code{NULL} (not recommended!), a temporary
    directory is used. This directory hosts:
    \itemize{
        \item a JSON file with run parameters called
        \code{params.json}.
        \item one subdirectory for each iteration which 
        contains intermediate files used e.g. by SNPTEST,
        PLINK, PRSice2 etc. These directories are named
        based on the \code{runId} and the iteration number.
        \item log file for each iteration, named similarly
        as the run subdirectory.
    }
    
    Regarding the \code{continue} argument, this can be
    set to \code{TRUE} to continue a crashed pipeline, 
    given the initial \code{gwe} argument and the prisma
    workspace. If the workspace is not found, a new run
    will be started. Generally, as the prisma pipeline
    combines several often disjoint components (including
    several system commands), it is possible that something
    can go wrong. We have done our best to modularize the
    pipeline which can be restarted from various 
    automatically detected checkpoints.
    
    Regarding the \code{useDenovoWorkspace} argument, if
    \code{useDenovoWorkspace} is a valid prisma workspace or
    a list output from \code{prsPipeline}, the PRS validation
    is performed using the original data split. If 
    \code{NULL}, then a completely new split is generated.
    The later can take very long time, so it is advised to 
    use an existing split if available. This argument is 
    available in the \code{prsPipeline} and 
    \code{prsSelection} functions when used directly. When
    executing the complete pipeline through \code{prisma},
    this is controlled by the \code{evalOnSplit} argument.
    When the latter is \code{evalOnSplit="original"}, the 
    workspace generated in \code{prisma} initiation is used. 
    When \code{evalOnSplit="new"}, the initial split is
    ignored and a new split is created and evaluated,
    taking much more time as the SNP effect calculation
    takes place again.
}
\examples{
data(toy,package="prisma")
gdsfile(toy) <- system.file(package="prisma","extdata","toy.gds")

# Lose filters for the toy dataset
filts <- getDefaults("filters")
filts$IBD <- NA
filts$hwe <- 1e-3
filts$pcaOut <- FALSE
filts$inbreed <- NA

# Workspace
d1 <- tempfile()
wspace1 <- file.path(tempdir(),d1)

# Traits
response <- "BMI"
covariates <- c("Case_Control","Gender","Age")

# Complete 4-iteration PRS pipeline - ~1 min to run with 4 cores - no PCA
if (Sys.which("PRSice_linux") != "") {

prismaOut <- prisma(
    gwe=toy,
    phenotype=response,
    covariates=covariates,
    pcs=FALSE,
    trainSize=0.8,
    niter=4,
    resolution="frequency",
    minSnps=2,
    filters=filts,
    pcaMethod="snprel",
    imputeMissing=FALSE,
    gwaMethods=c("glm","statgen"),
    prsiceOpts=list(clump_p=0.75,perm=100),
    prsWorkspace=wspace1,
    cleanup="intermediate",
    logging="file",
    output="normal",
    rc=0.25
)

# Below, individual wrapper steps implemented in prisma as whole
# Output only of the discovery PRS pipeline with statgenGWAS
d2 <- tempfile()
wspace2 <- file.path(tempdir(),d2)
dnList <- prsPipeline(
    gwe=toy,
    phenotype=response,
    covariates=covariates,
    pcs=FALSE,
    trainSize=0.8,
    niter=4,
    filters=filts,
    pcaMethod="snprel",
    imputeMissing=FALSE,
    gwaMethods="statgen",
    prsiceOpts=list(clump_p=0.75,perm=100),
    prsWorkspace=wspace2,
    cleanup="intermediate",
    logging="file",
    rc=0.25
)

# After, calculate the evaluation metrics
d3 <- tempfile()
wspace3 <- file.path(tempdir(),d3)
evalMetrics <- prsSelection(
    dnList=dnList,
    gwe=toy,
    phenotype=response,
    covariates=covariates,
    pcs=FALSE,
    trainSize=0.8,
    niter=4,
    resolution="frequency",
    minSnps=2,
    filters=filts,
    pcaMethod="snprel",
    imputeMissing=FALSE,
    gwaMethods="statgen",
    prsiceOpts=list(clump_p=0.75,perm=100),
    prsWorkspace=wspace3,
    cleanup="intermediate",
    logging="file",
    output="summary",
    useDenovoWorkspace=dnList,
    rc=0.25
)

candidates <- selectPrs(
    metrics=evalMetrics$metrics,
    snpSelection=evalMetrics$pgs,
    gwe=toy,
    base=evalMetrics$baseline
)

}

# Complete pipeline with report and checkpoints
# Complete 4-iteration PRS pipeline - ~1 min to run with 4 cores - no PCA
if (Sys.which("PRSice_linux") != "") {

outPath <- file.path(tempdir(),tempfile())
d4 <- tempfile()
wspace4 <- file.path(outPath,d4)
prismaPipeline(
    gwe=toy,
    phenotype=response,
    covariates=covariates,
    pcs=FALSE,
    trainSize=0.8,
    niter=4,
    resolution="frequency",
    minSnps=2,
    filters=filts,
    pcaMethod="snprel",
    imputeMissing=FALSE,
    gwaMethods=c("glm","statgen"),
    prsiceOpts=list(clump_p=0.75,perm=100),
    prsWorkspace=wspace4,
    cleanup="intermediate",
    logging="file",
    output="normal",
    rc=0.25,
    outPath=outPath
)

file.exists(file.path(outPath,"prismaOut.RData"))
# TRUE

}
}
\author{
    Panagiotis Moulos
}

