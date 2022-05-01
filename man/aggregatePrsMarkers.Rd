\name{aggregatePrsMarkers}
\alias{aggregatePrsMarkers}
\title{Aggregate PRS candidate SNPs from iterative extraction}
\usage{
    aggregatePrsMarkers(gwaList, mode = c("intersect", 
        "union"), qcut = 0.9, assoc = c("auto", "glm",
        "rrblup", "statgen", "snptest", "plink", "lasso"),
        avgfun = c("mean", "median", "weight"), gwe = NULL)
}
\arguments{
    \item{gwaList}{A list of GWASExperiments or their 
    summaries as derived from \code{\link{prsPipeline}}.}
    
    \item{mode}{Aggregation mode, can be \code{"intersect"}
    (default) or \code{"union"}.}
    
    \item{qcut}{A quantile cutoff for the frequency of SNP
    appearances in \code{gwaList}. SNPs with frequencies 
    below the one corresponding to this quantile are ignored
    from the output.}
    
    \item{assoc}{Which effects to use out of the performed
    GWA tests as the SNP weights in PRS. Can be \code{"auto"}
    (default) or any of \code{"glm"}, \code{"rrblup"}, 
    \code{"statgen"}, \code{"snptest"}, \code{"plink"}, as
    long as it exists. See also Details.}
    
    \item{avgfun}{How to summarize SNP effects (regression 
    coefficients) after the iterative extraction process.
    One of \code{"mean"} (default), \code{"mean"} or 
    \code{"weight"}. See also Details.}
    
    \item{gwe}{The GWASExperiment object used to derive the
    PRS. If given, the output is a more complete data frame,
    harmonized with the scores from EBI PGS Catalog.}
}
\value{
    If \code{gwe} is not \code{NULL}, a data frame of SNP
    markers summarized a \code{prsPipeline} run with the
    following columns:
    \itemize{
        \item \code{chromosome}: the chromosome where the 
        marker is located
        \item \code{position}: the marker position in
        chromosome
        \item \code{variant_id}: SNP rs id or other 
        identification where available
        \item \code{risk_allele}: The risk allele
        \item \code{reference_allele}: The non-risk allele
        \item \code{locus_name}: The SNP's closest gene
        (\code{NA} here, assigned later through 
        \code{\link{assignGenes}})
        \item \code{effect_weight}: The SNP effect as 
        derived from regression functions and \code{avgfun}
        \item \code{OR}: The \code{effect_weight} 
        transformed to Odds Ratio
        \item \code{asm}: The genome assembly (\code{"hg19"}
        or \code{"hg38"} or \code{"nr"} if unknown)
        \item \code{freq}: The number of appearances in the
        \code{prsPipeline} iterative process
    }
    
    If \code{gwe} is \code{NULL}, a simpler data frame of
    SNP marker with the following columns:
    \itemize{
        \item \code{snp} SNP rs id or other 
        identification where available
        \item \code{freq} The number of appearances in the
        \code{prsPipeline iterative process}
        \item \code{effect} The SNP effect as derived
        from regression functions and \code{avgfun}
    }
}
\description{
    This function aggregates the results of an iterative
    PRS extraction derived using \code{\link{prsPipeline}}.
    The aggregation is done by intersecting or unifying the
    SNPs extracted by more than one PRS algorithms. The
    aggregated result can be then subjected to a PRS 
    selection process according to the frequency of 
    appearance of each marker and result in a subset of
    the aggregated markers as the final PRS.
}
\details{
    Regarding the \code{assoc} argument, this controls which
    effects should be included in the final PRS data frame,
    out of the available ones. If \code{"auto"} (the 
    default), the final effects are selected based on the
    GWA algorithms that have been executed and an internal
    prioritization scheme according to our perception 
    regarding algorithm popularity. Specifically, the
    first available effects are used based on the following
    prioritization:
    \enumerate{
        \item SNPTEST (\code{"snptest"})
        \item PLINK (\code{"plink"})
        \item statgenGWAS (\code{"statgen"})
        \item GLM (\code{"glm"})
        \item rrBLUP (\code{"rrblup"})
    }
    In the future we intend to combine all available
    effects based on linear optimization methods (work in
    progress)
    
    Regarding the \code{avgfun} argument, this controls how
    the SNP effects in PRS candidates are summarized after 
    the iterative PRS extraction process in order to produce 
    the final effects for the calculation of weighted PRS. 
    When \code{"mean"} (default), the final effect for a SNP 
    is the mean off all its effects across the iterations
    and when it appears. When \code{"median"}, it is the 
    respective median while when \code{"weight"}, the SNP 
    effects are weighted by 1 minus the recorded PRS R^2 so 
    that the "best" runs can be given more value.
}
\examples{
    data(toy,package="prisma")
    data(prisma_out,package="prisma")
    
    iters <- prismaOut$results$glm$iterations
    #iters <- getRawIterations(prismaOut,"glm")
    snps <- aggregatePrsMarkers(iters,qcut=0.25,gwe=toy)
}
\author{
    Panagiotis Moulos
}

