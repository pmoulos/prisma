\name{twTest}
\title{Tracy-Widom test for selection of PCs}
\usage{
    twTest(eigv, p = 0.05, tol = 1e-8)
}
\arguments{
    \item{eigv}{}
    
    \item{p}{p-value cutoff for the test significance.}
    
    \item{tol}{}
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
