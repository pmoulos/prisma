\name{prsFuns}
\alias{lassosumPRS}
\alias{PRS}
\alias{prsicePRS}
\alias{runPRS}
\title{Functions to extract PRS candidates with various algorithms}
\usage{
    lassosumPRS(base, target = base, response, covariates = NULL,
        pcs = FALSE, wspace = NULL, anc = c("eur", "asn", "afr"),
        valid = c("auto", "stack", "split"), rc = NULL)
    
    PRS(gwe, snpSelection, type = c("avg", "sum", "std"))
    
    prsicePRS(base, target = base, response, covariates = NULL,
        pcs = FALSE, mode = c("calculate", "apply"), wspace = NULL,
        clump_kb = 250, clump_r2 = 0.1, clump_p = 1,
        score = c("avg", "sum", "std", "con-std"), perm = 10000,
        seed = 42, rc = NULL)
    
    runPRS(base, target = base, response, covariates = NULL,
        pcs = FALSE, methods = c("lassosum", "prsice"),
        lassosumOpts = getDefaults("lassosum"),
        prsiceOpts = getDefaults("prsice"), wspace = NULL, 
        rc = NULL)
}
\arguments{
    \item{base}{A \code{GWASExperiment} object with GWA
    tests performed and containing summary statistics apart
    from genotypes.}
    
    \item{target}{A \code{GWASExperiment} object with or
    without GWA tests, they are ignored. The PRS is 
    calculated based on the target dataset.}
    
    \item{response}{The trait for which a PRS is sought to be
    extracted. It should be a character of length 1 and exist
    in the phenotypes of \code{base} and \code{target}.}
    
    \item{covariates}{The covariates to include in regression
    models. It should be a character vector of traits all 
    existing in the phenotypes of \code{base} and 
    \code{target}.}
    
    \item{pcs}{Inlcude PCs in the regression model? Defaults
    to \code{FALSE}.}
    
    \item{mode}{The mode to run PRSice2. One of calculating and
    extracting a PRS (\code{"calculate"}, default) or applying
    a PRS given a GWASExperiment with existing \code{prsbetas}
    (\code{"apply"}).}
    
    \item{wspace}{A path to a workspace directory for 
    PRSice. If \code{NULL}, a temporary folder is used.}
    
    \item{anc}{One of the precalculated LD structures
    on lassosum, to be taken into account while extracting
    PRS. It can be \code{"eur"} (default) for European
    ancestry, \code{"asn"} for Asian or \code{"afr"} for
    African.}
    
    \item{valid}{Validation mode for lassosum PRS extraction. 
    One of \code{"auto"} (default), \code{"split"} or 
    \code{"stack"}. See also the lassosum package.}
    
    \item{gwe}{A code{GWASExperiment} object with the
    slot \code{effects} and/or \code{prsbetas} filled.}
    
    \item{snpSelection}{A data frame with candidate SNPs, the
    output, or part of the output of 
    \code{\link{aggregatePrsMarkers}}, or one of the list 
    members of the output of \code{\link{selectPrs}}.}
    
    \item{type}{Same as \code{score} in \code{prsicePRS}
    function (first three cases).}
    
    \item{clump_kb}{The genomic interval (in kbs) within 
    which PRSice2 performs clumping. Defaults to 
    \code{250}.}
    
    \item{clump_r2}{The R^2 threshold to keep SNPs after
    clumping in PRSice2. Defaults to \code{0.1}.}
    
    \item{clump_p}{The p-value threshold for SNPs to 
    include in clumping. Defaults to \code{1}.}
    
    \item{score}{How should PRSice2 calculate the PRS?
    One of \code{"avg"} (default), \code{"sum"}, 
    \code{"std"}, \code{"con-std"}. See also the manual
    of PRSice2 and Details.}
    
    \item{perm}{How many permutations should PRSice2
    perform? Defaults to 10000.}
    
    \item{seed}{Seed for PRSice2 permutations, for 
    reproducibility. Defaults to \code{42}.}
    
    \item{methods}{The methods to use for PRS extraction,
    one or more of \code{"prsice"} and \code{"lassosum"}.}
    
    \item{lassosumOpts}{List of options passed to lassosum
    (\code{lassosumPRS}). See also Details.}
    
    \item{prsiceOpts}{List of options passed to PRSice2
    (\code{prsicePRS}). See also Details.}
    
    \item{rc}{The fraction of available cores to use for 
    parallel calculations. Default is \code{NULL} for no
    parallelization.}
}
\value{
    For \code{lassosumPRS} and \code{prsicePRS}, a data
    frame with the SNPs comprising the extracted PRS,
    compatible with the score files returned by the PGS
    Catalog and the candidates returned by the complete
    \code{\link{prisma}} pipeline.
    
    For \code{runPRS}, the input \code{target} 
    GWASExperiment object with the \code{prsbetas}
    slots filled. Essentially the effect weights or
    ORs for the SNPs comprising the PRS.
}
\description{
    These functions perform the PRS extraction given one
    or two (base and target) \code{GWASExperiment} objects.
    \code{lassosumPRS} is a wrapper around the lassosum
    PRS extraction algorithm while \code{prsicePRS} is a
    wrapper around the PRSice2 PRS extraction algorithm.
    The \code{runPRS} is a wrapper around the two others
    and makes sure that the proper objects are retured.
}
\details{
    Regarding the \code{score} argument, it controls how
    PRSice2 calculates the PRS and creates regression models
    as well as how it reports it. Although typically it 
    should be the same as the \code{prsiceOpts$score}
    argument in the main \code{\link{prisma}} pipeline
    (and it is when \code{\link{prisma}} is used), it is
    not obligatory as PRISMA recalculates the final PRS
    before reporting and validating. However, if not the
    same, the actual score distribution will differ much
    between PRSice2 results and PRISMA results.
}
\examples{
    data(toy,package="prisma")
    gdsfile(toy) <- system.file(package="prisma","extdata","toy.gds")
    
    # Only 1 PC will be returned because of dataset size
    toy <- suppressWarnings(calcPcaCovar(toy,method="snprel"))
    
    toy <- gwa(toy,response="BMI",covariates=c("Age","Gender"),
        pcs=TRUE,psig=0.05,methods=c("glm","statgen"))
    
    prsiceOpts <- getDefaults("prsice")
    prsiceOpts$perm <- 100
    
    if (Sys.which("PRSice_linux") != "")
        toy <- runPRS(base=toy,response="BMI",
            covariates=c("Age","Gender"),pcs=TRUE,
            prsiceOpts=prsiceOpts)
    
    lprs <- lassosumPRS(base=toy,response="BMI",
        covariates=c("Age","Gender"),pcs=TRUE)
    
    if (Sys.which("PRSice_linux") != "")
        pprs <- prsicePRS(base=toy,response="BMI",
            covariates=c("Age","Gender"),pcs=TRUE,perm=100)
}
\author{
    Panagiotis Moulos
}
