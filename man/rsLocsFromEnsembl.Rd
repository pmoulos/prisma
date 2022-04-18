\name{rsLocsFromEnsembl}
\title{Fill up missing SNP locations using biomaRt}
\usage{
    rsLocsFromEnsembl(rs, gv = c("hg38", "hg19"), 
        canonical = TRUE)
}
\arguments{
    \item{rs}{a character vector of SNP (rs) ids.}
    
    \item{gv}{human genome version (\code{"hg38"}, 
    \code{"hg19"}).}
    
    \item{canonical}{check only canonical chromosomes? 
    Defaults to \code{TRUE}.}
}
\value{
    A data frame with the following columns:
    \itemize{
        \item \code{chr_name}: The chromosome name
        \item \code{chrom_start}: The SNP position
        \item \code{refsnp_id}: The reference SNP id 
        (usually dbSNP)
        \item \code{allele}: The major allele
        \item \code{minor_allele}: The minor allele
        \item \code{minor_allele_freq}: The recorded minor 
        allele frequency
        \item \code{chrom_strand}: The SNP strand
    }
}
\description{
    This function queries Ensembl to retrieve SNP info from
    a set of rs ids. It can be used not only to retrieve
    information but to salvage ambiguous or previously 
    unmapped SNPs in older datasets. Also heavily used from
    \code{\link{mergeGWAS}} function.
}
\examples{
\dontrun{
    # Requires active internet connection!
    rs <- c("rs1009238","rs2212890")
    info <- rsLocsFromEnsembl(rs)
}
}
\author{
    Panagiotis Moulos
}
