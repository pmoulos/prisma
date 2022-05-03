\name{calcPcaCovar}
\alias{calcPcaCovar}
\title{Calculate PCA covariates for GWASExperiments}
\usage{
    calcPcaCovar(obj, ld = 0.2, method = c("auto", 
        "snprel", "grid", "hubert"), npc = 0, pval = 0.05, 
        rc = NULL)
}
\arguments{
    \item{obj}{A \code{GWASExperiment} object without 
    missing values.}
    
    \item{ld}{Linkage Disequilibrium cutoff (defaults to
    \code{0.2}) for SNPs to include in PCA. See also
    Details.}
    
    \item{method}{Method used to calculate PCs. One of
    \code{"auto"} (default), \code{"snprel"}, \code{"grid"}
    or \code{"hubert"}. See also Details.}
    
    \item{npc}{The number of PCs to attach to \code{obj}.
    When \code{0} (default), the Tracy-Widom test is used
    to determine the best number of PCs to include.}
    
    \item{pval}{p-value for the Tracy-Widom test. Defaults
    to \code{0.05}. PCs with p-value below this threshold
    are ignored.}
    
    \item{rc}{Fraction of available cores to use when
    \code{method="snprel"}.}
}
\value{
    The input GWASExperiment object with attached PC
    covariates corresponding to the data hosted in the
    object. They can be retrieved with 
    \code{pcaCovariates(obj).}
}
\description{
    This function performs PCA using three available methods
    to the genetic data included in the iput object in order
    to capture potential underlying population structure or
    outliers, and derive PCs to be used as covariates in
    downstream GWA. Two of the three methods can be used to
    perform robust PCA with the intention of better 
    capturing potential outliers.
}
\details{
    Regarding the \code{ld} argument, it is generally
    recommended in the literature to use LD-pruned SNPs
    for PCA in GWAS Experiments. This argument controls
    the threshold of LD for LD-pruning. It is based on
    recommendations from Reed at al., 2013.
    
    Regarding the \code{method} argument, when 
    \code{"auto"}, if PCs already exist in the object, LD
    pruning is redone following by robust PCA if the 
    attached filters in the object have robust PCA marked.
    If not, PCA with the package SNPRelate is performed.
    It should be noted that for larger datasets, robust
    PCA can take a lot of time, while much less for the
    optimized method in SNPRelate.
}
\examples{
    data(toy,package="prisma")
    gdsfile(toy) <- system.file(package="prisma","extdata","toy.gds")
    
    # Only 1 PC will be returned because of dataset size
    toy <- suppressWarnings(calcPcaCovar(toy,method="snprel"))
}
\author{
    Panagiotis Moulos
}
