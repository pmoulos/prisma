\name{importGWAS}
\title{Create GWASExperiments from PLINK files}
\usage{
    importGWAS(input, phenos = NULL, backend = c("snpStats",
        "bigsnpr"), selection = NULL, genome = NA_character_,
        alleleOrder = c("plink", "reverse"),
        gdsfile = ifelse(backend == "snpStats", tempfile(), NA),
        gdsOverwrite = TRUE)
}
\arguments{
    \item{input}{}
    
    \item{phenos}{}
    
    \item{backend}{the R backend to use for object management.
    See Details.}
    
    \item{selection}{}
    
    \item{genome}{}
    
    \item{alleleOrder}{}
    
    \item{gdsfile}{}
    
    \item{gdsOverwrite}{overwrite the provided GDS file if
    exists? Defaults to \code{TRUE}.}
}
\value{
    
}
\description{
    
}
\details{
    
}
\examples{
    
}
\author{
    Panagiotis Moulos
}

