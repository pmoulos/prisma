\name{normalityCheck}
\title{Enrichment of PGS Catalog retrieved scores}
\usage{
    normalityCheck(x, pval = 0.05, tests = c("sw", "ks", "jb"),
        lower = 30, combine = c("fisher", "simes", "max", "min"), 
        nsample = 1000)
}
\arguments{
    \item{x}{a numeric vector to be tested for normality.}
    
    \item{pval}{p-value cutoff for the test significance.}
    
    \item{tests}{one or more of \code{"sw"} (Shapiro-Wilk),
    \code{"ks"} (Kolmogorov-Smirnov) and \code{"jb"}
    (Jarque-Bera).}
    
    \item{lower}{a lower threshold for sample size, below
    which the selected tests will be executed without 
    resampling. See also Details.}
    
    \item{combine}{p-value combination if multiple tests are
    selected. One of \code{"fisher"}, \code{"simes"}, 
    \code{"max"}, \code{"min"}. See also Details.}
    
    \item{nsample}{Number of resampling iterations when the
    sample size is larger than \code{lower}. Defaults to 
    \code{1000}.}
}
\value{
    A named list with members \code{statistics}, the actual
    values of the statistics used to assess normality, 
    \code{pvalues}, the respective p-values, \code{"pval"},
    the p-value threshold of the tests (same as input \code{p},
    \code{normal}, a logical indicating the output of the test.)
}
\description{
    This function uses three distinct tests to assess normality
    of a numeric vector, namely the Shapiro-Wilk test, the
    Kolmogorov Smirnov test against the normal distribution
    (\code{y="pnorm"} in the \code{\link{ks}} function) and the
    Jarque-Bera test implemented in the package tseries. If the
    latter is not present, the test is not available.
}
\details{
    Regarding the \code{lower} argument, it represents a 
    length threshold of the input vector (aka number of
    observations) above which the normality is assessed
    using resampling of \code{x} where each sample has 
    length \code{lower}. The final outcome is then 
    determined by the \code{p = 100 * (1 - 2*pval)} quantile 
    of the \code{nsample} length Shapiro-Wilk p-value 
    distribution created after \code{nsample} resampling
    iterations. If \code{p >= pval}, \code{x} can be 
    considered to follow the normal distribution.
    
    Regarding the \code{combine} argument, it determines how
    the outcomes of multiple tests will be combined. When
    \code{"fisher"} or \code{"simes"}, then the Fisher or
    Simes method for p-value combination is used. For
    details on these two methods please refer to 
    PMID: 25452340 and the Bioconductor package metaseqR2.
    When \code{"min"} or \code{"max"}, then the minimum
    (more lose) or the maximum (stricter) p-value out of
    three tests is considered to determine the outcome.
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
