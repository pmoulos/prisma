\name{harvestWorkspace}
\alias{harvestWorkspace}
\title{Reconstruct various prisma outputs from a workspace}
\usage{
    harvestWorkspace(wspace, rid, denovo = TRUE, fast = FALSE)
}
\arguments{
    \item{wspace}{path to a prisma workspace.}
    
    \item{rid}{a prisma run id.}
    
    \item{denovo}{does the workspace contain the results of
    a \emph{denovo} PRS pipeline execution?}
    
    \item{fast}{does the workspace contain the results of
    a PRS pipeline fast validation? (see \code{\link{prisma}}
    for further details.)}
}
\value{
    A list or data frame or prisma pipeline outcomes,
    according to the (auto-detected) run type.
}
\description{
    This function checks and traverses a path supposed to
    be a valid prisma workspace, either provided by the user
    or auto-generated, reads the existing results and returns
    a structure suitable for pipeline continuation or usable
    by other prisma functions. This function can be directly
    used by the user to break prisma steps for better 
    monitoring or continue a previously interrupted pipeline.
    It is extensively used when \code{continue=TRUE} in the
    respective pipelines that support continuation after a
    possible failure.
}
\details{
    It is important to note that prisma creates workspace
    paths with certain attributes that differentiate them
    from other random paths. One such attribute is that a
    prisma workspace is expected to contain at least:
    \itemize{
        \item a JSON file \code{params.json} with stored
        parameters for analysis continuation
        \item prisma iteration directories named with a
        certain convention according to the run type.
    }
}
\examples{
    wspace <- "/my/local/path"
    
    #output <- harvestWorkspace(wspace)
}
\author{
    Panagiotis Moulos
}
