\name{enrichScoreFile}
\title{Enrichment of PGS Catalog retrieved scores}
\usage{
    enrichScoreFile(sf, gb = c("hg19", "hg38", "nr"),
        clean = FALSE)
}
\arguments{
    \item{sf}{}
    
    \item{gb}{human genome build, \code{"hg19"} (default),
    for hg19/GRCh37, \code{"hg38"} for hg38/GRCh38 or 
    \code{"nr"} if unavailable.}
    
    \item{clean}{keep only matched SNPs with dbSNP and a
    genome version with respect to location. Default to
    \code{FALSE}.}
}
\value{
    
}
\description{
    
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
    
}
\author{
    Panagiotis Moulos
}

