\name{getDefaults}
\alias{getDefaults}
\title{Get default parameters for several prisma functions}
\usage{
    getDefaults(what)
}
\arguments{
    \item{what}{get parameters for what. It can be one of
    \code{"filters"}, \code{"glm"}, \code{"rrblup"},
    \code{"statgen"}, \code{"snptest"}, \code{"plink"},
    \code{"externalTools"}, \code{"lassosum"}, 
    \code{"prsice"}.}
}
\value{
    A named list with parameters corresponding to 
    \code{what} was asked.
}
\description{
    This function hard-codes default options for many of 
    prisma functionalities, such as genotype filtering,
    genome-wide association testing with various algorithms
    and polygenic risk score extraction algorithms.
}
\examples{
    filterDefaults <- getDefaults("filters")
    
    snptestDefaults <- getDefaults("snptest")
}
\author{
    Panagiotis Moulos
}
