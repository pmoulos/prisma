\name{GWASExperiment-class}
\docType{class}

% Classes
\alias{GWASExperiment}
\alias{class:GWASExperiment}
\alias{GWASExperiment-class}

% Accessors
\alias{genotypes}
\alias{genotypes,GWASExperiment-method}
\alias{genotypes<-}
\alias{genotypes<-,GWASExperiment-method}
\alias{genotypes<-,SnpMatrix-method}

\alias{gfeatures}
\alias{gfeatures,GWASExperiment-method}
\alias{gfeatures<-}
\alias{gfeatures<-,DataFrame-method}
\alias{gfeatures<-,data.frame-method}

\alias{gsamples}
\alias{gsamples,GWASExperiment-method}
\alias{gsamples<-}
\alias{gsamples<-,DataFrame-method}
\alias{gsamples<-,data.frame-method}

\alias{phenotypes}
\alias{phenotypes,GWASExperiment-method}
\alias{phenotypes<-}
\alias{phenotypes<-,DataFrame-method}
\alias{phenotypes<-,data.frame-method}
\alias{phenotypes<-,NULL-method}

\alias{pvalues}
\alias{pvalues,GWASExperiment-method}
\alias{pvalues<-}
\alias{pvalues<-,SimpleList-method}

\alias{effects}
\alias{effects,GWASExperiment-method}
\alias{effects<-}
\alias{effects<-,SimpleList-method}

\alias{prsbetas}
\alias{prsbetas,GWASExperiment-method}
\alias{prsbetas<-}
\alias{prsbetas<-,SimpleList-method}

\alias{genome}
\alias{genome,GWASExperiment-method}
\alias{genome<-}
\alias{genome<-,character-method}

\alias{gdsfile}
\alias{gdsfile,GWASExperiment-method}
\alias{gdsfile<-}
\alias{gdsfile<-,character-method}

\alias{filterRecord}
\alias{filterRecord,GWASExperiment-method}
\alias{filterRecord<-}
\alias{filterRecord<-,data.frame-method}

\alias{pcaCovariates}
\alias{pcaCovariates,GWASExperiment-method}
\alias{pcaCovariates<-}
\alias{pcaCovariates<-,matrix-method}

\alias{allPvalues}
\alias{allPvalues,SimpleList-method}
\alias{allEffects}
\alias{allEffects,SimpleList-method}
\alias{allPrsbetas}
\alias{allPrsbetas,SimpleList-method}

% Subsetting
\alias{[,GWASExperiment-method}
\alias{[,GWASExperiment,ANY-method}
\alias{[,GWASExperiment,ANY,ANY,ANY-method}
\alias{[<-,GWASExperiment,ANY,ANY,GWASExperiment-method}
\alias{subset,GWASExperiment-method}

% Quick colData access
\alias{[[,GWASExperiment,ANY,missing-method}
\alias{[[<-,GWASExperiment,ANY,missing-method}

% Display
\alias{show,GWASExperiment-method}

% Combine
\alias{rbind,GWASExperiment-method}
\alias{cbind,GWASExperiment-method}

\title{GWASExperiment objects}

\description{
    The GWASExperiment class is a matrix-like container which
    directly inherits all the properties from the 
    \code{\link{SummarizedExperiment}} class and extends it with
    respect to holding GWAS data in terms of accomodating not only 
    genotypic and phenotypic data but also association scores and
    SNP effects from Genome-Wide Association tests. In addition,
    as PRISMA's main purpose is the construction and evaluation
    of Polygenic Risk Scores, it offers facilities to store the
    results of respective analyses.
    
    Specifically, the assays of a GWASExperiment hold a matrix
    object (SnpMatrix or bigSNP) where rows represent SNPs or
    other variants and columns represent samples. The values of
    the assays represent genotypes and they are accessed and
    manipulated through the respective functions of the packages
    \code{snpStats} and \code{bigSNP}. Apart from the assays,
    the \code{rowData} of the GWASExperiment hold information
    about the markers and the \code{colData} hold information
    about the samples. Such information are typically derived
    from the initial input data used to construct the 
    GWASExperiment. The easiest way to construct it is from a
    set of PLINK files in BED, BIM and FAM formats. Direct
    construction is also possible, as long as the provided
    information about markers and samples is complete, for
    example, markers should be accompanied at least by
    chromosome and position data as well as major/minor 
    allele information.
    
    Apart from the basic GWASExperiment data which make direct
    use of the \code{\link{SummarizedExperiment}} slots, the
    class offers slots to store associated phenotypic and
    anthropometric data for each sample which are used for
    GWA tests, and slots to store the results of these tests
    in terms of p-values and reggression coefficients (effects).
    For the SNPs that may then become part of a PRS, the
    GWASExperiment also offers a slot for storing the weights
    of the PRS SNPs. Furthermore, GWASExperiment takes
    advantage of the \code{metadata} slot of the
    \code{\link{SummarizedExperiment}} to store additional
    information on GWAS and subsequent PRS extraction through
    the GWA summary statistics, such as Principal Component
    Analysis results used as covariates to capture and
    account for potential underlying population stratification,
    reference genome information, track of the filters applied
    to the genotypes and samples and other information. 
    Finally, it offers the ability to avoid genotypic and
    phenotypic data duplication when conducting multiple GWA
    tests by allowing the association and storage of many test
    results (e.g. with different covariates) within the same
    object.
}

\usage{
    ## Constructor

    GWASExperiment(
        genotypes = SnpMatrix(), 
        features = DataFrame(),
        samples = DataFrame(),
        phenotypes = DataFrame(),
        pvalues = SimpleList(),
        effects = SimpleList(),
        prsbetas = SimpleList(),
        metadata = list(
            genome = NA_character_,
            backend = NA_character_,
            filters = setNames(
                data.frame(matrix(ncol = 5, nrow = 0)),
                c("parameter", "name", "value", "type", "filtered")
            )
        )
    )


    ## Accessors - apart from the ones of SummarizedExperiment

    genotypes(x, ...)
    genotypes(x, ...) <- value
    gfeatures(x, withDimnames=TRUE, ...)
    gfeatures(x, withDimnames=TRUE, ...) <- value
    gsamples(x, withDimnames=TRUE, ...)
    gsamples(x, withDimnames=TRUE, ...) <- value

    phenotypes(x, withDimnames=TRUE, ...)
    phenotypes(x, withDimnames=TRUE, ...) <- value

    pvalues(x, response = 1L, covariates = NULL, npcs = 0,
        withDimnames = TRUE, ...)
    pvalues(x, response = 1L, covariates = NULL, npcs = 0,
        withDimnames = TRUE, ...) <- value
    allPvalues(x, i = seq_len(ncol(x)))
    effects(x, response = 1L, covariates = NULL, npcs = 0,
        withDimnames = TRUE, ...)
    effects(x, response = 1L, covariates = NULL, npcs = 0,
        withDimnames = TRUE, ...) <- value
    allEffects(x, i = seq_len(ncol(x)))
    prsbetas(x, response = 1L, covariates = NULL, npcs = 0,
        withDimnames = TRUE, ...)
    prsbetas(x, response = 1L, covariates = NULL, npcs = 0,
        withDimnames = TRUE, ...) <- value
    allPrsbetas(x, i = seq_len(ncol(x)))

    genome(x)
    genome(x) <- value
    gdsfile(x)
    gdsfile(x) <- value
    filterRecord(x)
    filterRecord(x) <- value
    pcaCovariates(x)
    pcaCovariates(x) <- value


    ## Quick colData access

    \S4method{$}{GWASExperiment}(x, name)
    \S4method{$}{GWASExperiment}(x, name) <- value
    \S4method{[[}{GWASExperiment,ANY,missing}(x, i, j, ...)
    \S4method{[[}{GWASExperiment,ANY,missing}(x, i, j, ...) <- value


    ## Subsetting

    \S4method{[}{GWASExperiment}(x, i, j, ..., drop=TRUE)
    \S4method{[}{GWASExperiment,ANY,ANY,GWASExperiment}(x, i, j) <- value


    ## Combining

    \S4method{cbind}{GWASExperiment}(..., deparse.level=1)
    \S4method{rbind}{GWASExperiment}(..., deparse.level=1)
}

\arguments{
    \item{x}{a GWASExperiment object.}
  
    \item{...}{for \code{cbind}, \code{rbind}, \code{...} 
    contains GWASExperiment objects to be combined. For other 
    accessors, ignored.}
  
    \item{value}{an object of a class specified in the S4 method
    signature or as outlined in Details.}

    \item{i, j}{for \code{genotypes}, \code{genotypes<-}, 
    \code{i} is an integer or numeric scalar, see Details 
    for additional constraints.

    For \code{[,GWASExperiment}, \code{[,GWASExperiment<-}, 
    \code{i}, \code{j} are subscripts that can act to subset the 
    rows and columns of \code{x}, that is the \code{matrix} 
    elements of \code{genotypes}.

    For \code{[[,GWASExperiment}, \code{[[<-,GWASExperiment}, 
    \code{i} is a scalar index (e.g., \code{character(1)} or 
    \code{integer(1)}) into a column of \code{gsamples}.}

    \item{withDimnames}{a \code{logical(1)}, indicating whether 
    the dimnames of the GWASExperiment object should be applied 
    (i.e. copied) to the extracted assays (genotypes or other
    data such as samples, features, phenotypes). More precisely, 
    setting \code{withDimnames=FALSE} in the \emph{getter} 
    returns the assays \emph{as-is} whereas setting 
    \code{withDimnames=FALSE} return them with possibly modified 
    dimnames.

    Setting \code{withDimnames=FALSE} in the \emph{setter} 
    (\code{assays<-}) is required when the dimnames on the supplied 
    assays are not identical to the dimnames on the GWASExperiment 
    object. It does not influence actual assignment of dimnames to 
    assays (they're always stored as-is). For more details, see
    \code{\link{SummarizedExperiment}}.}
    
    \item{drop}{A \code{logical(1)}, ignored by these methods.}

    \item{deparse.level}{See \code{?base::\link[base]{cbind}} for a 
    description of this argument.}
    
    \item{genotypes}{a \code{SnpMatrix} or \code{bigSNP} 
    matrix-like object that holds the genotypes for the 
    SNPs/markers in each row and each sample in each column.}
    
    \item{features}{a \code{\link{DataFrame}} or 
    \code{\link{data.frame}} holding detailed information about
    each genotyped SNP. See also Details.}
    
    \item{samples}{a \code{\link{DataFrame}} or 
    \code{\link{data.frame}} holding information about each
    genotyped sample. See also Details.}
    
    \item{phenotypes}{a \code{\link{DataFrame}} or 
    \code{\link{data.frame}} holding phenotypic records about 
    each genotyped sample, typically used for GWA tests. See 
    also Details.}
    
    \item{pvalues}{a \code{\link{SimpleList}} or 
    \code{\link{list}} with p-values from GWA tests for each
    SNP. Each member corresponds to a GWA test with specific
    covariates. See also Details.}
    
    \item{effects}{a \code{\link{SimpleList}} or 
    \code{\link{list}} with regression coefficients (or converted
    odds ratios) from GWA tests for each SNP (effects). Each 
    member corresponds to a GWA test with specific covariates. 
    See also Details.}
    
    \item{prsbetas}{a \code{\link{SimpleList}} or 
    \code{\link{list}} with weights for each SNP participating
    in a PRS. Each member corresponds to a GWA test with specific
    covariates. Typically, subset of effects, unless an adjustment
    algorithms (e.g. lassosum) has been applied. See also Details.}
    
    \item{metadata}{a list similar to \code{\link{metadata}}
    of \code{\link{SummarizedExperiment}}. It should contain at
    least \code{genome}, \code{backend} and \code{filters}
    members. See also Details.}
    
    \item{response}{a numeric or integer or character scalar.
    If a number, it should correspond to the column number of
    the response phenotype in a GWA test. If character, it
    should be one of the phenotypes associated with \code{x}.
    See also Details.}
    
    \item{covariates}{a character vector or \code{NULL}. The
    covariates to be used in regression models and GWA tests.
    They should exist in the phenotypes of \code{x}. See also 
    Details.}
    
    \item{npcs}{the number of PCs included in the covariates
    in a GWA test. See also Details.}
    
    \item{i}{a scalar or vector of indices for which to retrieve
    statistical scores, corresponding to a number of models fit
    to \code{x}. See also Details.}
}

\details{
    The GWASExperiment class is a straighforward expansion of the
    \code{\link{SummarizedExperiment}} class which is meant for 
    numeric and other data types derived typically from a sequencing 
    experiment but is suitable to host genotypic data as the main 
    features of the latter are Single Nucleotide Polymorphisms 
    (SNPs) which share a lot of common information with sequencing
    experiments such as genomic coordinates and certain types of
    abundance. The structure is rectangular like a \code{matrix}, 
    but with additional annotations on the rows and columns, and 
    with the possibility to manage several assays simultaneously 
    so long as they be of the same dimensions. More generic
    information can be found on the manual pages of the
    \link{SummarizedExperiment} class. The following focus on the
    expansions introduced in the GWASExperiment class.

    The rows of a GWASExperiment object represent SNPs/markers of
    interest. Information about these markers is stored in a
    \link{DataFrame} or \link{data.frame} object, accessible 
    using the function \code{gfeatures} (or the native 
    \code{rowData} as \code{gfeatures} is an alias to it). The 
    object must have as many rows as there are rows in the 
    GWASExperiment object, with each row of the object providing 
    information on the SNP in the corresponding row of the 
    GWASExperiment object. Columns of the object represent 
    different attributes of the features of interest. Specifically, 
    for the GWASExperiment case, it must have at least genomic
    coordinates, that is chromosome and position, a unique id (it
    can be a dbSNP id) and the two alleles comprising the genotype.
    A typical case is when the GWASExperiment is created from PLINK
    files with the backend provided by \code{snpStats} package
    and the columns of the object are \code{chromosome},
    \code{snp.name}, \code{cM}, \code{position}, \code{allele.1}
    and \code{allele.2}. These are also the preferred column names
    of the object. It must also have assigned rownames, preferrably
    the same as the \code{snp.names} column.

    Each column of a GWASExperiment object represents a sample
    (e.g. a person contributing to a GWAS project). Information 
    about the samples are stored in a \link{DataFrame} or ,
    \link{data.frame} object, accessible using the function 
    \code{gsamples} (or the native \code{colData} as 
    \code{gsamples} is an alias to it). The object must have as 
    many rows as there are columns in the GWASExperiment object, 
    with each row of the object providing information on the sample
    in the corresponding column of the GWASExperiment object.
    Columns of the object represent different sample attributes, 
    e.g., pedigree, sex etc. Column names typically provide a short 
    identifier unique to each sample. A typical case is when the 
    GWASExperiment is created from PLINK files with the backend 
    provided by \code{snpStats} package and the columns of the
    object are \code{pedigree}, \code{member}, \code{father},
    \code{mother}, \code{sex} and \code{affected}. These are also 
    the preferred column names of the object. It must also have 
    assigned rownames, preferrably the same as the \code{member} 
    column.

    A GWASExperiment object, as a SummarizedExperiment, can also 
    contain information about the overall experiment, for instance 
    the lab in which it was conducted, the publications with which 
    it is associated, etc. This information is stored as a 
    \code{list} object, accessible using the \code{metadata} 
    function. The form of the data associated with the experiment
    is left to the discretion of the user. However, the 
    \code{metadata} of a GWASExperiment object must have at least
    three elements. These are:
    
    \itemize{
        \item \code{genome}: the reference genome used in the
        genotyping technology (usually \code{hg38} or \code{hg19}).
        It can be accessed and set using the \code{genome} function.
        
        \item \code{filters}: the set of filters applied to the
        GWASExperiment to filter out noisy genotyping data (if any)
        and a summary of their results in a \code{data.frame}. A
        placeholder is automatically added on construction. It can
        be accessed and set with the \code{filterRecord} fuction.
        
        \item \code{backend}: the backend package used to read and
        store the genotype data (one of \code{snpStats} or
        \code{bigsnpr}). This is crucial for downstream functions
        and manipulation and must not be changed.
    }

    The GWASExperiment container is appropriate for matrix-like
    genotypic data, such as those derived from gentype microarray
    technologies. Apart from the related SummarizedExperiment data
    accessor functions, the data are accessed using the 
    \code{genotypes} function, described below. This returns a
    \link{SnpMatrix} or \link{bigSNP} object. A \link{SimpleList} 
    object with genotypic data from potentially different 
    experiment can be accessed with \code{\link{assays}}. 
    \code{genotypes} returns the first one, but typically, this is 
    what is being used. Each element of the list must itself either 
    a \code{SnpMatrix} or a {bigSNP} object and must have 
    dimensions that are the same as the dimensions of the
    GWASExperiment in which they are stored. Row and column
    names of each object must either be \code{NULL} or match those 
    of the GWASExperiment during construction.
    
    Furthermore, the GWASExperiment offers facilities for storing
    other data typical for subsequent GWAS analyses and PRS
    constructions, such as storing phenotypic and anthropometric
    data used for association analysis, the results of GWA tests
    and the results of PRS analysis. Instead of storing these
    data to the metadata slot of a SummarizedExperiment, they
    are bound to slots following the dimensions and subsetting
    of the main genotypic matrix object. Therefore, phenotypic
    data are accessed through the \code{phenotypes} function,
    the results of GWA tests (p-values, regression coefficients)
    are accessed through the \code{pvalues} and \code{effects}
    functions respectively, and associated PRS weights are
    accessed throught the function \code{prsbetas}. All these
    functions are described below and all return 
    \code{SimpleList} objects, as a GWASExperiment object may
    have multiple GWA and PRS results attached to it to avoid
    genotype repetition. A listing of associated analyses
    with a GWASExperiment object may be retrieved with the
    \code{listGwa} and \code{listPrs} functions and the 
    functions \code{allPvalues}, \code{allEffects} and
    \code{allPrsbetas} can be used to retrieve a list of
    all the analyses results.
}

\section{Constructor} {
    GWASExperiment instances are constructed using the
    \code{GWASExperiment} function with arguments outlined 
    above.
}

\section{Accessors} {
    In the following code snippets, \code{x} is a GWASExperiment
    object. We do not describe standard \link{SummarizedExperiment}
    accessors such as \code{rowData}, \code{colData}, \code{dim}
    etc. Details can be found in the respective man pages. General
    patterns such as that row data must have \code{NULL} or identical
    row names as with the main matrix object (here, the genotypes)
    are assumed.

\describe{
    \item{\code{genotypes(x)}, \code{genotypes(x) <- value}:}{Get 
    or set the genotypes. \code{value} is a \code{SnpMatrix} or 
    \code{bigSNP} with the same dimensions as \code{x}.}

    \item{\code{gfeatures(x, withDimnames=TRUE)}, 
    \code{gfeatures(x) <- value}:}{Get or set the feature data (similar
    to row data in a SummarizedExperiment). \code{value} is a 
    \code{\link{DataFrame}} or a \code{\link{data.frame}} object.}
    
    \item{\code{gsamples(x, withDimnames=TRUE)}, 
    \code{gsamples(x) <- value}:}{Get or set the sample data (similar
    to column data in a SummarizedExperiment). \code{value} is a 
    \code{\link{DataFrame}} or a \code{\link{data.frame}} object.}
    
    \item{\code{phenotypes(x, withDimnames=TRUE)}, 
    \code{phenotypes(x) <- value}:}{Get or set the phenotypic data
    for each sample. \code{value} is a \code{\link{data.frame}} 
    object with as many rows as the number of samples and associated
    phenotypes in the experiment.}
    
    \item{\code{pvalues(x, response = 1L, covariates = NULL, 
    npcs = 0, withDimnames = TRUE)}, pvalues(x, response = 1L, 
    covariates = NULL, npcs = 0, withDimnames = TRUE, ...) <- value}:{
    Get or set the p-values after the application of one or more
    genotype-phenotype association tests. \code{response} is the
    phenotype for which the association test is performed. It can be
    a numeric or integer value denoting the column number in the 
    phenotypes data frame which is associated with x \code{value},
    or a character in which case it's the name of the phenotype and
    corresponds to a column name of the aforementioned data frame. 
    It defaults to the first available phenotype. \code{covariates}
    are the phenotypic covariates to be included in an association 
    model and must be available in the phenotypes data frame. It can
    be either a numeric/integer vector denoting column numbers or
    a character vector denoting column names of the phenotypes data
    frame. It defaults to \code{NULL} for no addition of covariates,
    other than the main phenotype (\code{response}). \code{npcs}
    denotes the number of PCs to include in the model as additional
    covariates to capture potential population substructures. It
    defaults to \code{0} for no PCs as additional covariates. For
    \code{withDimnames}, see above e.g. \code{gfeatures}. The
    reason for this accessor to have additional arguments
    (\code{response}, \code{covariates}, \code{npcs}) is that
    \code{x} may have multiple association test results attached
    to it, and these variable offer a way of unique access to
    them. Finally, \code{value} is a numeric matrix with p-values
    for each SNP. This matrix has as many columns as the different
    test types performed (e.g. SNPTEST and statgenGWAS) and are
    named accordingly.}
    
    \item{\code{allPvalues(x, i = seq_len(ncol(x)))} will
    retrieve (access-only) either all the p-value matrices
    attached to \code{x} or as many as requested by \code{i},
    assuming at least \code{i} tests have been performed.}
    
    \item{\code{effects(x, response = 1L, covariates = NULL, 
    npcs = 0, withDimnames = TRUE)}, effects(x, response = 1L, 
    covariates = NULL, npcs = 0, withDimnames = TRUE, ...) <- value}:{
    Get or set the SNP effects with respect to the phenotype after 
    the application of one or more genotype-phenotype association 
    tests. The effects are essentially regression coefficients or
    odds ratios for binomial regression models. \code{response}, 
    \code{covariates} and \code{npcs} are as in the \code{pvalues}
    accessor. \code{value} is a numeric matrix with effects for
    each SNP. This matrix has as many columns as the different test
    types performed (e.g. SNPTEST and statgenGWAS) and are named
    accordingly.}
    
    \item{\code{allEffects(x, i = seq_len(ncol(x)))} will
    retrieve (access-only) either all the effect matrices
    attached to \code{x} or as many as requested by \code{i},
    assuming at least \code{i} tests have been performed.}
    
    \item{\code{prsbetas(x, response = 1L, covariates = NULL, 
    npcs = 0, withDimnames = TRUE)}, prsbetas(x, response = 1L, 
    covariates = NULL, npcs = 0, withDimnames = TRUE, ...) <- value}:{
    Get or set the beta coefficients (or odd ratios) effects 
    derived after the application of a Polygenic Risk Scores
    algorithm, with respect to the phenotype tested with one or
    more of the association tests. Essentially, the PRS betas are
    the regression coefficients for the SNPs participating in a
    PRS either corrected in the case of effect shrinkage algorithms
    (e.g. lassosum), or as returned by the association test.
    \code{response}, \code{covariates} and \code{npcs} are as in 
    the \code{pvalues} accessor. \code{value} is a numeric matrix 
    with effects for each SNP participating in a PRS candidate and
    0 otherwise. This matrix has as many columns as the different 
    PRS algorithms executed on \code{x} and are named accordingly.}
    
    \item{\code{allPrsbetas(x, i = seq_len(ncol(x)))} will
    retrieve (access-only) either all the PRS beta matrices
    attached to \code{x} or as many as requested by \code{i},
    assuming at least \code{i} tests have been performed.}
    
    \item{\code{genome(x), \code{genome(x) <- value}:}{Get or set
    the genome build where the features (SNPs, other variations)
    are mapped to. Typically, one of \code{"hg19"} or \code{"hg38"}
    as the class is designed mainly for human data. Other genomes
    may be provided, but for these human builts, other 
    functionalities are also supported such as SNP information
    retrieval from dbSNP and/or Ensembl, and also enrichment of
    PGS score files.}
    
    \item{\code{gdsfile(x), \code{gdsfile(x) <- value}:}{Get or set
    the supporting GDS file location. This is required by the 
    package SNPRelate to perform certain functionalities such as
    Linkage Disequilibrium (LD) pruning, Identity By Descent (IBD)
    filtering and others. If not provided, these functionalities
    will not be available.}
    
    \item{\code{filterRecord(x), \code{filterRecord(x) <- value}:}{
    Get or set a data frame holding a record of the filters applied
    to \code{x}, the cutoff values and the filtering results (i.e.
    how many SNPs where filtered). This data frame has a specific
    format and violating will cause an error. Typically, users will
    only access the filter records. To see the structure of this
    data frame, call \code{filterRecord(x)}.}
    
    \item{\code{pcaCovariates(x), \code{pcaCovariates(x) <- value}:}{
    Get or set the PCA covariates that capture underlying population
    stratification structures. The PCA covariates are typically
    calculated using the \code{\link{calcPcaCovars}} function with
    one of the supporting methods within the package. Users will
    usually only access the PCA covariates while setting them using
    the \code{\link{calcPcaCovars}} function.}
}

\section{Subsetting}{
    The subsetting operations are as described in the respective
    section of the \link{SummarizedExperiment} help page. There 
    are restrictions in the case of a GWASExperiment as in the
    case of replacement, the data types stored in its slots must
    be compatible (i.e. SnpMatrix, bigSNP objects for genotypes 
    and DataFrame, data.frame objects for associated genotype and
    sample data).
}

\section{Combining}{
    The combining functions of GWASExperiment objects follow the
    same principles as those of SummarizedExperiment objects.
    Special attention is given regarding the integrity of the main
    SummarizedExperiment class extensions. These are briefly the 
    phenotypes which follow the samples (columns/colData) of
    the GWASExperiment object and statistical scores parallel to
    the features (rowData) of the object (pvalues, effects etc.).
    
    When combining objects, the following attributes and values are
    not preserved, as they become meaningless with the combination
    of samples and/or SNPs:
    
    \itemize{
        \item \code{pvalues, effects, prsbetas}: these scores are
        dropped as the initial population in potentially both 
        samples and SNPs from which they were derived is changed.
        
        \item \code{filters}: any filtered SNPs/samples statistics
        are re-initialized as the population of SNPs/samples is
        changed.
        
        \item any values related to PCA, as again, these statistics
        depend on the population which changes with object 
        combination.
    }
}

\author{
    Panagiotis Moulos
}

\seealso{
    \link{SummarizedExperiment} and everything therein
}

\examples{
nrows <- 200; ncols <- 6
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
colData <- DataFrame(Treatment=rep(c("ChIP", "Input"), 3),
                     row.names=LETTERS[1:6])
se0 <- SummarizedExperiment(assays=SimpleList(counts=counts),
                            colData=colData)
se0
dim(se0)
dimnames(se0)
assayNames(se0)
head(assay(se0))
assays(se0) <- endoapply(assays(se0), asinh)
head(assay(se0))

rowData(se0)
colData(se0)

se0[, se0$Treatment == "ChIP"]
subset(se0, select = Treatment == "ChIP")

## cbind() combines objects with the same features of interest
## but different samples:
se1 <- se0
se2 <- se1[,1:3]
colnames(se2) <- letters[seq_len(ncol(se2))]
cmb1 <- cbind(se1, se2)
dim(cmb1)
dimnames(cmb1)

## rbind() combines objects with the same samples but different
## features of interest:
se1 <- se0
se2 <- se1[1:50,]
rownames(se2) <- letters[seq_len(nrow(se2))]
cmb2 <- rbind(se1, se2)
dim(cmb2)
dimnames(cmb2)
}
