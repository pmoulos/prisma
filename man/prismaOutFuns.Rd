\name{prismaOutFuns}
\alias{getPrsCandidates}
\alias{gwaTests}
\alias{prismaReport}
\alias{prsEvalReport}
\title{Functions to get data from PRISMA outputs}
\usage{
    getPrsCandidates(prismaOut, method, index = NULL)
    
    gwaTests(prismaOut)
    
    prismaReport(gwe, prismaOut, cvMetricsOut, lookupOut,
        path = NULL)
    
    prsEvalReport(evalList, path = NULL)
}
\arguments{
    \item{prismaOut}{An output object from \code{\link{prisma}}
    or a list of PRS candidate data frames.}
    
    \item{method}{A GWA method supported, see also the 
    \code{\link{gwa}} function.}
    
    \item{index}{An integer denoting which candidate to get
    out of the available ones. If \code{NULL}, a list of all
    candidates is returned.}
    
    \item{gwe}{The GWASExperiment object used to generate the
    \code{prismaOut}.}
    
    \item{cvMetricsOut}{The output from 
    \code{\link{prismaCrossValidate}} function.}
    
    \item{lookupOut}{The output from \code{\link{prismaLookup}} 
    function.}
    
    \item{evalList}{A list of evaluation metrics of the PRS 
    when applied to a new dataset, essentially a (named) list 
    of outputs from one or multiple runs of the 
    \code{\link{applyPRS}} function.}
    
    \item{path}{The path to write the interactive prisma 
    report. Defaults to \code{NULL} for auto-generation in
    the current path (\code{getwd()}).}
}
\value{
    A data frame or a list of data frame with SNPs comprising
    the PRS candidate(s). In the case of \code{prismaReport}
    and \code{prsEvalReport}, an interactive HTML report.
}
\description{
    The \code{getPrsCandidates} function facilitates the
    retrieval of PRS candidates (data frames) from the more
    complex output of \code{link{prisma}}.
    
    The \code{gwaTests} function lists the GWA methods that
    have been used in a \code{\link{prisma}} call.
    
    The \code{prismaReport} function produces the main PRISMA
    interactive report of the PRS and evaluations with the
    original dataset.
    
    The \code{prsEvalReport} function produces a short HTML
    report of evaluation metrics for PRS candidates when 
    applied to another dataset.
}
\examples{
    data(prisma_out,package="prisma")
    
    can <- getPrsCandidates(prismaOut,"glm",1)
    gwaTests(prismaOut)
}
\author{
    Panagiotis Moulos
}
