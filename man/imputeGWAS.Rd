\name{imputeGWAS}\title{Genotype imputation for GWASExperiment objects}
\usage{
    imputeGWAS(obj, mode = c("single", "split"),
        failed = c("scrime", "none"), rc = NULL)
}
\arguments{
    \item{obj}{A \code{GWASExperiment} object with missing
    genotype values to be imputed.}
    
    \item{mode}{Perform the imputation within the same
    dataset (\code{"single"}, default) or split the dataset
    into a part containing missing values and a part with
    no missing values (\code{"split"}). See also Details and
    Description for further explanation of the internal 
    imputation process.}
    
    \item{failed}{What to do with genotypes that failed to 
    be imputed with the iterative regression process? Impute
    with kNN from the scrime package (\code{"scrime"}, 
    default) or no further imputation (\code{"none"}).}
    
    \item{rc}{The fraction of available cores to use for 
    parallel calculations. Default is \code{NULL} for no
    parallelization.}
}
\value{
    A \code{GWASExpriment} object with no 
    (\code{failed="scrime"}) missing values or fewer
    (\code{failed="none"}) missing values in the genotypes.
}
\description{
    This function performs internal imputation in a
    \code{GWASExperiment} object using a combined approach
    with fucntionalities from the snpStats packages and the
    scrime package. The internal imputation in PRISMA refers
    to missing genotype imputation within the same dataset
    and using existing genotypes _within_ the same dataset
    to impute missing genotypes. It should not be confused
    with imputation based on an external reference panel
    such as the 1000 Genomes Project which can be used to
    both impute and extend a dataset to include 
    non-genotyped SNPs/markers. This process is also 
    available with PRISMA and is refered as "extension"
    (see also the \code{\link{extendGWAS}} function.)
}
\details{
    Regarding the \code{mode} argument, this controls how
    the imputation based on snpStats functions 
    (\code{\link{snp.imputation}}) will be performed.
    When \code{mode="single"}, the imputation rules are
    created within the total initial dataset without any
    splitting to a "better" and a "worse" dataset in terms
    of missing value content. When \code{mode="split"},
    the dataset is split in a part containing samples with
    missing values and a part with samples without missing
    values. While this is not realistic for smaller datasets
    it may happen in large dataset and could be more accurate.
    If splitting is chosen for small datasets, it's likely
    that the imputation rules will be created with a very
    small number of samples, resulting in fewer rules (i.e.
    fewer SNPs can be imputed) and with poorer accuracy.
    
    The internal imputation process in PRISMA is split
    in two phases:
    \enumerate{
        \item Phase 1: The genotypes of the 
        \code{GWASExperiment} input is used to create
        imputation rules (with the respective function of
        snpStats). These rules are essentially linear
        regression models which model potential predictions
        of SNP genotypes according to nearby genotypes and
        present high R^2 values and therefore can be used
        to provide missing cases with relatively good 
        accuracy. This process is iteratively deployed in 
        PRISMA, until no more imputation rules can be
        produced, i.e. the result of each imputation
        iteration is used to create new rules to impute
        more SNPs if possible, until no more new rules
        can be produced.
        
        \item Phase 2: In this phase, the missing genotypes
        that could not be imputed with the iterative process
        of Phase 1, are imputed with the k-nearest neighbor
        algorithm for genotype imputation from the R package
        scrime. This step is optional.
        
        \item Phase 3: Missing genotypes that could not be
        imputed with the scrime package in Phase 2 are
        replaced with the global average genotype of the
        rest of the genotypes. This phase is optional and
        should be quite rare to be needed.
    }
}
\examples{
    data(toy,package="prisma")
    
    # Make some genotypes missing
    set.seed(42)
    toMiss <- lapply(seq_len(ncol(toy)),function(j) {
        n <- sample(20,1)
        miss <- logical(nrow(toy))
        miss[sample(nrow(toy),n)] <- TRUE
        return(miss)
    })
    
    G <- as(genotypes(toy),"numeric") + 1
    mode(G) <- "integer"
    for (j in seq_len(ncol(toy)))
        G[toMiss[[j]],j] <- 0L
    G <- SnpMatrix(G)
    genotypes(toy) <- G
    
    toyI <- imputeGWAS(toy,fail="none")
    before <- length(which(is.na(genotypes(toy))))
    after <- length(which(is.na(genotypes(toyI))))
    before > after # TRUE
}
\author{
    Panagiotis Moulos
}

