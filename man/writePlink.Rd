\name{writePlink}
\title{PLINK files from GWASExperiment}
\usage{
    writePlink(obj, pheno = NULL, outBase = NULL, salvage = FALSE,
        reverse = FALSE)
}
\arguments{
    \item{obj}{an object of GWASExperiment class}
    
    \item{pheno}{the name of the phenotype to be written to
    the respective FAM column.}
    
    \item{outBase}{}
    
    \item{salvage}{try to salvage SNPs with missing locations?
    Defaults to \code{FALSE}. See also Details.}
    
    \item{reverse}{}
}
\value{
    
}
\description{
    
}
\details{
    
}
\examples{
    x <- rnorm(20)
    
    o1 <- normalityCheck(x)
    o2 <- normalityCheck(x,combine="simes")
    # Two instead of 3 tests
    o3 <- normalityCheck(x,tests=c("sw","jb"),combine="simes")
}
\author{
    Panagiotis Moulos
}
