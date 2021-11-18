\name{twTest}
\title{Enrichment of PGS Catalog retrieved scores}
\usage{
    normalityCheck(x, pval = 0.05, tests = c("sw", "ks", "jb"),
        lower = 30, forceTest = TRUE, combine = c("fisher",
        "simes", "max", "min"), nsample = 1000)
}
\arguments{
    \item{x}{}
    
    \item{pval}{p-value cutoff for the test significance.}
    
    \item{tests}{}
    
    \item{lower}{}
    
    \item{forceTest}{}
    
    \item{combine}{}
    
    \item{nsample}{}
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
