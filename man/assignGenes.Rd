\name{assignGenes}
\title{Assign genes to a PRISMA output or list of PRS candidates}
\usage{
    assignGenes(prismaOut, gv = c("hg19","hg38"), sitDb = NULL)
}
\arguments{
    \item{prismaOut}{An output object from \code{\link{prisma}}
    or a list of PRS candidate data frames.}
    
    \item{gv}{Human genome version, one of \code{"hg19"} (default)
    or \code{"hg38"}.}
    
    \item{sitDb}{A sitadela SQLite database with at least the
    two genomes described above from RefSeq.}
}
\value{
    Same as input, with HUGO gene names assigned to the 
    \code{locus_name} column of PRS candidates.
}
\description{
    
}
\details{
    If the sitadela database is not found in the provided location,
    the required information will be downloaded on-the-fly. This is
    not advised as the annotation is done iteratively. When 
    \code{sitDb=NULL}, no annotation will be performed.
}
\examples{
    
}
\author{
    Panagiotis Moulos
}
