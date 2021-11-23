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

    \S4method{$}{SummarizedExperiment}(x, name)
    \S4method{$}{SummarizedExperiment}(x, name) <- value
    \S4method{[[}{SummarizedExperiment,ANY,missing}(x, i, j, ...)
    \S4method{[[}{SummarizedExperiment,ANY,missing}(x, i, j, ...) <- value


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
    object.

\describe{
    \item{\code{genotypes(x)}, \code{genotypes(x) <- value}:}{Get 
    or set the genotypes. \code{value} is a \code{SnpMatrix} or 
    \code{bigSNP} with the same dimensions as \code{x}.}

    \item{\code{rowData(x, use.names=TRUE)}, \code{rowData(x) <- value}:}{
      Get or set the row data. \code{value} is a \link{DataFrame} object.}

    \item{\code{colData(x)}, \code{colData(x) <- value}:}{Get or set the
      column data. \code{value} is a \link{DataFrame} object. Row
      names of \code{value} must be NULL or consistent with the existing
      column names of \code{x}.}

    \item{\code{metadata(x)}, \code{metadata(x) <- value}:}{Get or set
      the experiment data. \code{value} is a \code{list} with arbitrary
      content.}

    \item{\code{dim(x)}:}{Get the dimensions (features of interest x samples)
      of the SummarizedExperiment.}

    \item{\code{dimnames(x)}, \code{dimnames(x) <- value}:}{Get or set
      the dimension names. \code{value} is usually a list of length 2,
      containing elements that are either \code{NULL} or vectors of
      appropriate length for the corresponding dimension. \code{value}
      can be \code{NULL}, which removes dimension names. This method
      implies that \code{rownames}, \code{rownames<-}, \code{colnames},
      and \code{colnames<-} are all available.}

  }
}

\section{Subsetting}{

  In the code snippets below, \code{x} is a SummarizedExperiment object.

  \describe{

    \item{\code{x[i,j]}, \code{x[i,j] <- value}:}{Create or replace a
      subset of \code{x}. \code{i}, \code{j} can be \code{numeric},
      \code{logical}, \code{character}, or \code{missing}. \code{value}
      must be a SummarizedExperiment object with dimensions,
      dimension names, and assay elements consistent with the subset
      \code{x[i,j]} being replaced.}

    \item{\code{subset(x, subset, select)}:}{Create a subset of \code{x}
      using an expression \code{subset} referring to columns of
      \code{rowData(x)} and / or \code{select} referring to column names
      of \code{colData(x)}.}

  }

  Additional subsetting accessors provide convenient access to
  \code{colData} columns

  \describe{

    \item{\code{x$name}, \code{x$name <- value}}{Access or replace
    column \code{name} in \code{x}.}

    \item{\code{x[[i, ...]]}, \code{x[[i, ...]] <- value}}{Access or
    replace column \code{i} in \code{x}.}

  }

}

\section{Combining}{

  In the code snippets below, \code{x}, \code{y} and \code{...} are
  SummarizedExperiment objects to be combined.

  \describe{

    \item{\code{cbind(...)}:}{
      \code{cbind} combines objects with the same features of interest
      but different samples (columns in \code{assays}).
      The colnames in \code{colData(SummarizedExperiment)} must match or
      an error is thrown.
      Duplicate columns of \code{rowData(SummarizedExperiment)} must
      contain the same data.

      Data in \code{assays} are combined by name matching; if all assay
      names are NULL matching is by position. A mixture of names and NULL
      throws an error.

      \code{metadata} from all objects are combined into a \code{list}
      with no name checking.
    }

    \item{\code{rbind(...)}:}{
      \code{rbind} combines objects with the same samples
      but different features of interest (rows in \code{assays}).
      The colnames in \code{rowData(SummarizedExperiment)} must match or
      an error is thrown.
      Duplicate columns of \code{colData(SummarizedExperiment)} must
      contain the same data.

      Data in \code{assays} are combined by name matching; if all assay
      names are NULL matching is by position. A mixture of names and NULL
      throws an error.

      \code{metadata} from all objects are combined into a \code{list}
      with no name checking.
    }

    \item{\code{combineRows(x, ..., use.names=TRUE, delayed=TRUE, fill=NA)}:}{
      \code{combineRows} acts like more flexible \code{rbind}, returning a
      SummarizedExperiment with features equal to the concatenation of features
      across all input objects. Unlike \code{rbind}, it permits differences in
      the number and identity of the columns, differences in the available
      \code{\link{rowData}} fields, and even differences in the available
      \code{\link{assays}} among the objects being combined.

      If \code{use.names=TRUE}, each input object must have non-\code{NULL},
      non-duplicated column names. These names do not have to be the same, or
      even shared, across the input objects. The column names of the returned
      \code{SummarizedExperiment} will be a union of the column names across
      all input objects. If a column is not present in an input, the
      corresponding assay and \code{colData} entries will be filled with
      \code{fill} and \code{NA}s, respectively, in the combined
      SummarizedExperiment.

      If \code{use.names=FALSE}, all objects must have the same number of
      columns. The column names of the returned object is set to
      \code{colnames(x)}. Any differences in the column names between input
      objects are ignored.

      Data in \code{assays} are combined by matching the names of the assays.
      If one input object does not contain a named assay present in other input
      objects, the corresponding assay entries in the returned object will be
      set to \code{fill}. If all assay names are NULL, matching is done by
      position. A mixture of named and unnamed assays will throw an error.

      If \code{delayed=TRUE}, assay matrices are wrapped in
      \code{\link{DelayedArray}}s to avoid any extra memory allocation during
      the matrix \code{rbind}ing. Otherwise, the matrices are combined as-is;
      note that this may still return \code{DelayedMatrix}s if the inputs were
      also \code{DelayedMatrix} objects.

      If any input is a \code{RangedSummarizedExperiment}, the returned object
      will also be a \code{RangedSummarizedExperiment}. The \code{rowRanges} of
      the returned object is set to the concatenation of the \code{rowRanges}
      of all inputs. If any input is a \code{SummarizedExperiment}, the
      returned \code{rowRanges} is converted into a \code{GRangesList} and the
      entries corresponding to the rows of the \code{SummarizedExperiment} are
      set to zero-length \code{GRanges}. If all inputs are
      \code{SummarizedExperiment} objects, a \code{SummarizedExperiment} is
      also returned.

      \code{rowData} are combined using \code{\link{combineRows}} for
      \code{DataFrame} objects. It is not necessary for all input objects to
      have the same fields in their \code{rowData}; missing fields are filled
      with \code{NA}s for the corresponding rows in the returned object.

      \code{metadata} from all objects are combined into a \code{list}
      with no name checking.
    }

    \item{\code{combineCols(x, ..., use.names=TRUE, delayed=TRUE, fill=NA)}:}{
      \code{combineCols} acts like more flexible \code{cbind}, returning a
      SummarizedExperiment with columns equal to the concatenation of columns 
      across all input objects. Unlike \code{cbind}, it permits differences in
      the number and identity of the rows, differences in the available
      \code{\link{colData}} fields, and even differences in the available
      \code{\link{assays}} among the objects being combined.

      If \code{use.names=TRUE}, each input object must have non-\code{NULL},
      non-duplicated row names. These names do not have to be the same, or
      even shared, across the input objects. The row names of the returned
      \code{SummarizedExperiment} will be a union of the row names across
      all input objects. If a row is not present in an input, the
      corresponding assay and \code{rowData} entries will be filled with
      \code{fill} and \code{NA}s, respectively, in the combined
      SummarizedExperiment.

      If \code{use.names=FALSE}, all objects must have the same number of rows.
      The row names of the returned object is set to \code{rownames(x)}. Any
      differences in the row names between input objects are ignored.

      Data in \code{assays} are combined by matching the names of the assays.
      If one input object does not contain a named assay present in other input
      objects, the corresponding assay entries in the returned object will be
      set to \code{fill}. If all assay names are NULL, matching is done by
      position. A mixture of named and unnamed assays will throw an error.

      If \code{delayed=TRUE}, assay matrices are wrapped in
      \code{\link{DelayedArray}}s to avoid any extra memory allocation during
      the matrix \code{rbind}ing. Otherwise, the matrices are combined as-is;
      note that this may still return \code{DelayedMatrix}s if the inputs were
      also \code{DelayedMatrix} objects.

      If any input is a \code{RangedSummarizedExperiment}, the returned object
      will also be a \code{RangedSummarizedExperiment}. The \code{rowRanges} of
      the returned object is set to a merge of the \code{rowRanges} of all
      inputs, where the coordinates for each row are taken from the input
      object that contains that row. Any conflicting ranges for shared rows
      will raise a warning and all \code{rowRanges} information from the
      offending \code{RangedSummarizedExperiment} will be ignored. If any
      input is a \code{SummarizedExperiment}, the returned \code{rowRanges} is
      converted into a \code{GRangesList} and the entries corresponding to the
      unique rows of the \code{SummarizedExperiment} are set to zero-length
      \code{GRanges}. If all inputs are \code{SummarizedExperiment} objects, a
      \code{SummarizedExperiment} is also returned.

      \code{colData} are combined using \code{\link{combineRows}} for
      \code{DataFrame} objects. It is not necessary for all input objects to
      have the same fields in their \code{colData}; missing fields are filled
      with \code{NA}s for the corresponding columns in the returned object.

      \code{metadata} from all objects are combined into a \code{list}
      with no name checking.
    }
  }

}

\section{Implementation and Extension}{

  This section contains advanced material meant for package developers.

  SummarizedExperiment is implemented as an S4 class, and can be extended in
  the usual way, using \code{contains="SummarizedExperiment"} in the new
  class definition.

  In addition, the representation of the \code{assays} slot of
  SummarizedExperiment is as a virtual class Assays. This
  allows derived classes (\code{contains="Assays"}) to implement
  alternative requirements for the assays, e.g., backed by file-based
  storage like NetCDF or the \code{ff} package, while re-using the existing
  SummarizedExperiment class without modification.
  See \link{Assays} for more information.

}

\author{Martin Morgan; \code{combineRows} and \code{combineCols} by Aaron Lun}

\seealso{
  \itemize{
    \item \link{RangedSummarizedExperiment} objects.

    \item \link[S4Vectors]{DataFrame}, \link[S4Vectors]{SimpleList}, and
          \link[S4Vectors]{Annotated} objects in the \pkg{S4Vectors} package.

    \item The \code{\link[S4Vectors]{metadata}} and
          \code{\link[S4Vectors]{mcols}} accessors in the \pkg{S4Vectors}
          package.

    \item \code{\link[HDF5Array]{saveHDF5SummarizedExperiment}} and
          \code{\link[HDF5Array]{loadHDF5SummarizedExperiment}} in the
          \pkg{HDF5Array} package for saving/loading an HDF5-based
          SummarizedExperiment object to/from disk.

    \item The \code{\link[DelayedArray]{realize}} generic function in the
          \pkg{DelayedArray} package for more information about on-disk
          realization of objects carrying delayed operations.
  }
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

## ---------------------------------------------------------------------
## ON-DISK REALIZATION
## ---------------------------------------------------------------------
library(DelayedArray)
setAutoRealizationBackend("HDF5Array")
cmb3 <- realize(cmb2)
assay(cmb3, withDimnames=FALSE)  # an HDF5Matrix object
}
