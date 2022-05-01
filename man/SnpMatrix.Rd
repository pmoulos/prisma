\name{SnpMatrix}
\alias{SnpMatrix}
\title{Simple SnpMatrix object constructor}
\usage{
    SnpMatrix(snp) 
}
\arguments{
    \item{snp}{A matrix of genotypes with integers
    or numerics. See also Details.}
}
\value{
    A SnpMatrix object.
}
\description{
    This funtion serves as a constructor of a SnpMatrix
    object from a matrix of coded genotypes.
}
\details{
    Regarding the main \code{snp} argument, if this is an
    integer matrix, genotypes should be coded as \code{1}
    (homozygous minor allele), \code{2} (heterozygous),
    \code{3} (homozygous major allele). Missing values
    should be represented as \code{0}s. If this is a
    numeric matrix, genotypes should be coded as \code{0}
    (homozygous minor allele), \code{1} (heterozygous),
    \code{2} (homozygous major allele). Missing values
    should be represented as \code{NA}s.
}
\examples{
    g <- z <- cbind(c(1L,2L,3L),c(0L,1L,2L),c(3L,2L,1L))
    S <- SnpMatrix(g)
    as(S,"numeric")
}
\author{
    Panagiotis Moulos
}
