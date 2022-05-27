\name{mergeGWAS}
\alias{mergeGWAS}
\title{Merge GWASExperiments from different populations}
\usage{
    mergeGWAS(gwe1, gwe2, output = c("common", "all"), 
        mafResolve = NULL, gdsfile = NULL, writegds = FALSE)
}
\arguments{
    \item{gwe1}{A \code{GWASExperiment} object from the first
    dataset.}
    
    \item{gwe2}{A \code{GWASExperiment} object from the second
    dataset.}
    
    \item{output}{One of \code{"common"} (default) to operate
    on common SNPs between \code{gwe1} and \code{gwe2} or 
    \code{"all"} to operate on all SNPs. SNPs unique to each
    dataset are returned with missing genotypes in the
    merged set.}
    
    \item{mafResolve}{How would alleles that cannot be 
    resolved just by strand fliping be resolved? Based on
    the frequency in one of the two datasets (\code{1} or
    \code{2}) or from online resources (dbSNP, population
    studies etc.). The later is the default (\code{NULL}).}
    
    \item{gdsfile}{path to the new GDS file. If \code{NULL}, a 
    temporary file will be written (unadvised if 
    \code{writegds=TRUE}!).}
    
    \item{writegds}{write the new GDS file? Defaults to 
    \code{FALSE}. It can be written anytime with the
    \code{\link{GWASExperiment2GDS}} function.}
}
\value{
    A new GWASExperiment object with merged genotypes,
    samples and phenotypes.
}
\description{
    This function attempts to merge two GWASExperiments with
    various goals, such as merging similar populations with
    common phenotypic and/or anthropometric data to increase
    statistical power. See Details for a more rigorous 
    explanation of the merging process.
}
\details{
    The merging of \code{GWASExperiment} objects is performed
    based on the following steps and assumptions. Firstly, the
    two objects are checked for genome version compatibility.
    If they are not compatible the process does not start.
    You should consider using the 
    \code{\link{GWASExperimentLiftOver}} function that can
    potentially render the two datasets compatible. Also,
    both datasets should have a common reference SNP naming.
    It is not possible to merge a dataset where the markers
    are named e.g. by Illumina IDs and a dataset where the
    markers are named with dbSNP (rs) IDs.
    
    If the two datasets are compatible in terms of genome 
    version and coordinates, the common SNPs are then
    identified between the two datasets. As these must
    point to the same information (location, minor and major
    alleles, same strand), the process then checks these
    attributes and tries to resolve any discrepancies.
    
    The first check is whether non identical SNPs in terms
    of alleles can be solved with strand-fliping. Those that
    cannot be resolved with strand-fliping may not point to
    the same risk allele. This is resolved either by using
    one dataset as a reference (\code{mafResolve=1 or 2})
    or even better, by querying online resources (e.g.
    Ensembl and dbSNP) to identify the true minor and major
    alleles based on summarized frequencies from population
    studies (there should be an internet connection for 
    this). SNPs whose alleles cannot be resolved by any
    means (strand-flip, online MAF resolve) are completely
    dropped.
    
    If the common SNPs between the two objects continue not
    being identical, it probably means that there are 
    location/position mismatches potentially derived from
    older annotations. These can only be resolved online.
    If after this there are any remaining non-resolvable
    SNPs, they are dropped. Then, the genotypes of common
    SNPs are merged. The genotypes of the alleles that were 
    flipped are also switched.
    
    Finally, the common annotation elements and genotypes
    are combined. If \code{output=="common"} the process
    stops. If \code{output="all"}, the unique elements of
    each dataset are also combined with the common. The
    genotypes of SNPs unique to each dataset are set to
    missing.
    
    With respect to the allele resolving of the common
    SNPs that cannot be resolved with strand-flipping 
    online:
    \enumerate{
        \item SNPs are queried with biomaRt to resolve 
        allele frequencies. Genome version is not important 
        at this point as only alleles are used and IDs.
        \item If minor allele is available from Ensembl, 
        then the flip decision is based on it. If not (e.g. 
        outdated/merged SNPs, dbSNP is queried with package
        rsnps and fill the Ensembl result.
        \item Some SNPs may have been removed from dbSNP. 
        Nothing can be done. Flip for these will be based 
        on the majority of the rest.
        \item Some found SNPs may not have information on 
        MAF and minor alleles. dbSNP is requeried with 
        rsnps and filling missing info is retried. Again, 
        some will not be resolved and will be flipped 
        according to majority.
    }
}
\examples{
    data(toy,package="prisma")
    
    toy1 <- toy[1:800,1:50]
    toy2 <- toy[201:1000,51:100]
    
    toym <- mergeGWAS(toy1,toy2)
    dim(toym) # 600 100
    
    toym <- mergeGWAS(toy1,toy2,output="all")
    dim(toym) # 1000 100
}
\author{
    Panagiotis Moulos
}
