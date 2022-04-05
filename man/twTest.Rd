\name{twTest}
\title{Tracy-Widom test for selection of PCs}
\usage{
    twTest(eigv, p = 0.05, tol = 1e-8)
}
\arguments{
    \item{eigv}{The eigen values calculated from a PCA.
    See also Details.}
    
    \item{p}{p-value cutoff for the test significance.}
    
    \item{tol}{A value that will replace eigen values which
    are smaller than this value. Required for internal
    calculations.}
}
\value{
    A named list with members: \code{tw} containing the
    TW statistic for each eigen value, \code{pval} 
    containing the p-values and \code{index} containing
    the index of the significant eigen values, aka PCs.
}
\description{
    This function implements the Tracy-Widom statistic
    for the assessment of the most significant principal
    components of a Principal Component Analysis. It is
    implemented based on assumptions described in Zhao
    et al., PMID: 30507552
}
\details{
    Regarding the \code{eigv} input, it can be retrieved
    according to which method or package has been used for
    PCA. For example, if robust PCA has been used (package
    rrcov), the \code{getEigenvalues} function can be used
    to get the eigen values. If \code{\link{prcomp}} has
    been used, then \code{eigv} is the squared \code{sdev} 
    component of the \code{prcomp} output.
}
\examples{
# Matrix borrowed from prcomp man page
C <- chol(S <- toeplitz(0.9^(0:31)))
X <- matrix(rnorm(32000),1000,32)
Z <- X %*% C

P <- prcomp(Z)
tw <- twTest(P$sdev^2)
# tw$index provides the significant PCs
identical(tw$index,which(tw$pval<0.05)) # TRUE
}
\author{
    Panagiotis Moulos
}
