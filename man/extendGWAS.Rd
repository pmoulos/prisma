\name{extendGWAS}
\alias{extendGWAS}
\title{External imputation of GWAS data using IMPUTE2}
\usage{
    extendGWAS(obj, intSize = 1e+6, wspace = NULL, 
        refSpace = NULL, continue = FALSE,
        cleanup = c("none", "intermediate", "all"),
        runId = NULL, convTool = c("gtool", "qctool"),
        rc = NULL)
}
\arguments{
    \item{obj}{A \code{GWASEXperiment} object.}
    
    \item{intSize}{The interval size where IMPUTE2
    will operate. Defaults to \code{1e+6}. See also
    Details.}
    
    \item{wspace}{The path to a directory that will be used
    as IMPUTE2 workspace. If \code{NULL} a temporary 
    directory is used, however no intermediate files can
    be saved in this way. It is therefore advised to
    explicitly provide a path.}
    
    \item{refSpace}{The local directory where the IMPUTE2
    reference panel files are located. See also the 
    \code{link{download1000GP3}} function. If \code{NULL}
    it is assumed that the files have been downloaded in
    the default path used by the aforementioned function.}
    
    \item{continue}{Continue a crashed imputation session?
    Defaults to \code{FALSE} for the initial run. See also
    Details.}
    
    \item{cleanup}{The cleanup level after the imputation 
    process. \code{"none"} (default) leaves all the 
    intermediate files intact while \code{"intermediate"}
    removes interval GEN files and \code{"all"} removes the
    entire \code{wspace}.}
    
    \item{runId}{An identification string for the run. It
    is used to create a subdirectory inside \code{wspace}.
    If not provided (\code{NULL}, default), a random string
    is generated.}
    
    \item{convTool}{The conversion tool to use between GEN
    and PLINK file formats. Defaults to \code{"gtool"}. See
    also Details.}
    
    \item{rc}{The fraction of available cores to use. These
    are used to impute SNPs in multiple genomic intervals
    at the same time.}
}
\value{
    A (larger) GWASExperiment object with the imputed SNPs
    based on the IMPUTE2 reference panel. Bare in mind that
    the SNP ids that were not present in the initial
    GWASExperiment object do not follow the same naming 
    format but they have the names of the SNPs present
    in the reference panel used by IMPUTE2. The 
    \code{gfeatures} of the output will be augmented by
    one column with the INFO score returned by IMPUTE2 and
    denoting the "trustworthiness" of the imputed position.
}
\description{
    This function is a convenient wrapper around the IMPUTE2
    program for imputing genotypes based on an initial 
    dataset of genotyped SNPs and a reference panel based on
    population studies such as the 1000 genomes project.
}
\details{
    Regarding the \code{intSize} argument, this controls
    the genomic interval (on the same chromosome) where
    IMPUTE2 will impute reference panel variants based on
    the variants that exist on the same interval in
    \code{obj}. By default IMPUTE2 does not operate on
    whole chromosomes as the process will not be as accurate
    as operating on intervals and will take longer. Using
    intervals also IMPUTE2 is parallelized according to the
    \code{rc} parameter.
    
    Regarding the \code{continue} argument, this can be used
    to continue a crashed imputation session, provided that
    the path to the IMPUTE2 working space is given 
    (\code{wspace}). The function then automatically 
    determines what has been already done and restarts from
    a specific checkpoint. This funcionality is provided as 
    the imputation process may take a long time, especially
    for \code{GWASExperiment}s containing a lot of samples,
    and something could go wrong resulting in the stop of
    the process.
    
    Regarding the \code{"cleanup"} argument, the intermediate
    files produced are PLINK files and GEN files produced and
    and used by IMPUTE2 and also intermediate GEN files per
    imputation interval. These are removed in the case of
    \code{cleanup="intermediate"}, where also the final 
    imputed PLINK files are left on the workspace. If the
    workspace is not deleted, the output of the function can
    be recostructed using also the \code{runId}, otherwise, 
    the resulting \code{GWASExperiment} should be stored in 
    an RDA file for later usage.
    
    Regarding the \code{convTool} argument, this controls
    which external tool will be used for the conversion of
    PLINK files to GEN files and vice-versa. The default is
    \code{"gtool"} but we have experienced crashed therefore
    support for \code{"qctool"} was added. See also the
    \code{\link{downloadGtool}} and 
    \code{\link{downloadQctool}} functions.
}
\examples{
\dontrun{
    # The process may take long or crash because of the
    # small sample size.
    data(toy,package="prisma")
    toyEx <- extendGWAS(toy)
}
}
\author{
    Panagiotis Moulos
}
