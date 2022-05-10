\name{ebiQueryFuns}
\alias{getGWASVariants}
\alias{getPGSScores}
\alias{prismaLookup}
\title{Functions to query EBI GWAS and PGS Catalogs}
\usage{
    getGWASVariants(efoId = NULL, efoTrait = NULL, 
        removeUnknownRisk = TRUE, retries = 5)
    
    getPGSScores(pgsId = NULL, efoId = NULL, pubmedId = NULL,
        base = NULL, retries = 5, validateLoc = FALSE)
    
    prismaLookup(prismaOut)
}
\arguments{
    \item{efoId}{One or more EFO ids in proper format. See
    also Details.}
    
    \item{efoTrait}{One or more EFO traits in free text 
    format, which should be existing. See also Details.}
    
    \item{pgsId}{One or more PGS ids in proper format. See
    also Details.}
    
    \item{pubmedId}{One or more PubMed ids in proper format.}
    
    \item{retries}{The number of retries in case of API call
    failures.}
    
    \item{removeUnknownRisk}{Should variants of unknown risk
    be removed in a score from GWAS Catalog? Defaults to
    \code{TRUE}.}
    
    \item{base}{The base path for the PGS Catalog location.
    It can be a local path, see also Details.}
    
    \item{validateLoc}{Validate SNP locations in retrieved
    PGS scores? Defaults to \code{FALSE}, see also Details.}
    
    \item{prismaOut}{An output object from \code{\link{prisma}}
    or a list of PRS candidate data frames.}
}
\value{
    A data frame such as the one from 
    \code{\link{aggregatePrsMarkers}} when used with a
    \code{GWASExperiment} object as accompanying input,
    with as many values filled as possible (e.g. GWAS
    Catalog does not record variant regression effects.)
    
    For \code{prismaLookup}, a list of lists. The outer
    level is named with the GWA methods in \code{prismaOut}.
    The inner level is named according to the PRS candidate
    names for each GWA method in \code{prismaOut}.
}
\description{
    These functions can be used to retrieve PRS candidates
    in the format of data frames (compatible with the general
    \code{\link{prisma}} output as well as the one from 
    \code{\link{aggregatePrsMarkers}}). This can be done from
    two sources, the EBI PGS Catalog and the EBI GWAS Catalog
    given one or more identifiers. These functions build upon
    the gwasrapidd and quincunx packages which leverage the
    EBI APIs to fetch the required data. An additional 
    feature regarding PGS scores is that while quincunx only
    fetches metadata, \code{getPGSScores} actually fetches
    and combines (if possible) the PGS scores, with alleles,
    effects etc. Then, the resulting data frames can be used
    with a \code{GWASExperiment} object and either a PRISMA
    external validation pipeline with PRScice2 (takes longer)
    or with vanilla R regression and cross-validation using
    e.g. \code{\link{prsCrossValidate}} to check their
    performance against one or more traits present in the
    \code{GWASExperiment} object.
}
\details{
    Regarding the \code{efoId} argument, this has to be one
    or more valid EFO ids from the GWAS Catalog in proper
    format (e.g. \code{"EFO_0004566"}). Improper format is 
    detected and rejected by the respective function. The 
    same applies to the \code{pgsId} argument which should be 
    like \code{"PGS000034"} and the \code{pubmedId} format 
    which should be like \code{"32527150"}.
    
    Regarding the \code{base} argument, this can be either
    \code{NULL} (the default) to retrieve the PGS score from
    the link present in the quincunx API call metadata, or
    a local path following the same structure as the PGS
    catalog FTP storage. The first case is more guaranteed
    but prone to connection errors, depending also on the
    operation system, and especially when multiple score
    files must be retrieved and processed. The second case
    assumes that the scores of the PGS Catalog are stored
    somewhere locally (they do not occupy much space).
    This can be done for example with Linux \code{wget}
    command:
    
    \preformatted{
        mkdir -p LOCAL_PATH
        cd LOCAL_PATH
        wget -m ftp://anonymous:@ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/
        
        # Then in R:
        # base <- LOCAL_PATH/ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/
        # If the long path is not desired:
        
        cd LOCAL_PATH/ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/
        mv * ../../../../../../
        cd LOCAL_PATH
        rm -r ./ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/
        
        # Then in R: base <- LOCAL_PATH
    }
    
    Regarding the \code{validateLoc} argument, in case of 
    duplicates remaining after basic position/id/risk allele 
    deduplication, when \code{TRUE}, performs location 
    validation, as in some cases, duplicate rs ids have 
    different locations (especially when the genome assembly 
    is not reported (\code{"nr"}) and potentially indicating 
    very old genome assemblies (<hg18). The process relies on 
    the package biomaRt and uses it to query the duplicate rs 
    ids and cross-validate the locations retrieved from 
    Ensembl variation databases. SNPs with locations that 
    cannot be validated are dropped and deduplicated. As the 
    process involves calling Ensembl APIs, it ma be slow. 
    Therefore, if the intended use of the retrieved PRS is 
    to plugin the SNPs to other GWAS data, \code{validateLoc} 
    can be \code{FALSE} to speed up the process and reduce 
    failures, as the unique SNP ids are used with local 
    genotypes instead of location. In other cases where 
    strict validation is required (e.g. when ORs must be 
    used or reported, then it should be \code{TRUE}.
}
\examples{
\donttest{
    # Assume an internet connection!
    efoId <- c("EFO_0004566","EFO_0005937")
    efoTrait <- c("body weight gain",
        "longitudinal BMI measurement") 
    
    test1 <- getGWASVariants(efoId=efoId,
        removeUnknownRisk=FALSE)
    test2 <- getGWASVariants(efoTrait=efoTrait)
    
    pgsId <- "PGS000034"
    test3 <- getPGSScores(pgsId=pgsId)
    
    data(prisma_out,package="prisma")
    lookup <- prismaLookup(prismaOut)
}
}
\author{
    Panagiotis Moulos
}
