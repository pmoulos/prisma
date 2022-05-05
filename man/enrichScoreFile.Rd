\name{enrichScoreFile}
\alias{enrichScoreFile}
\title{Enrichment of PGS Catalog retrieved scores}
\usage{
    enrichScoreFile(sf, gb = c("hg19", "hg38", "nr"),
        clean = FALSE)
}
\arguments{
    \item{sf}{A data frame created from a file retrieved 
    from PGS Catalog.}
    
    \item{gb}{human genome build, \code{"hg19"} (default),
    for hg19/GRCh37, \code{"hg38"} for hg38/GRCh38 or 
    \code{"nr"} if unavailable.}
    
    \item{clean}{keep only matched SNPs with dbSNP and a
    genome version with respect to location. Default to
    \code{FALSE}.}
}
\value{
    The same data frame with whatever information (such
    as location, rs ids etc.) was missing and was 
    retrieved.
}
\description{
    This function will try to enrich the PRS SNPs with 
    information from the respective genome version if
    such info is not already there (e.g. some PGS files
    are missing locations or rs ids)
}
\details{
    When \code{"gb"} is provided (i.e. retrieved from PGS
    Catalog with the respective \code{quincunx} API call),
    the function will try to enrich the PRS SNPs with 
    information from the respective genome version if
    such info is not already there (e.g. some PGS files
    are missing locations or rs ids). If not available
    (\code{"nr"}) or is simply is not known, the function
    will try to use enrichment information from hg38.
}
\examples{
    scoreFile <- system.file(package="prisma","extdata","PGS000299.txt.gz")
    pgs <- read.delim(scoreFile,comment.char="#")
    pgse <- enrichScoreFile(pgs)
    is(pgse,"GPos") # TRUE
}
\author{
    Panagiotis Moulos
}
