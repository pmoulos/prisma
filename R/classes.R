# There will probably be another class for storing retrieved from public APIs

setClassUnion("data.frame_OR_DataFrame_OR_NULL",
    c("data.frame","DataFrame","NULL"))

# Definition
.GWASExperiment <- setClass("GWASExperiment",
    contains="SummarizedExperiment",
    slots=c(
        phenotypes="data.frame_OR_DataFrame_OR_NULL"
    )
)

# Constructor
GWASExperiment <- function(genotypes=SnpMatrix(),features=DataFrame(),
    samples=DataFrame(),phenotypes=DataFrame(),
    metadata=list(genome=NA_character_,backend=NA_character_,
        filters=setNames(data.frame(matrix(ncol=4,nrow=0)),
            c("parameter","name","value","type","filtered"))),...) {
    
    # It MUST have metadata with certain names in the list
    if (!is.list(metadata)) {
        warning("The provided metadata argument is not a list! ",
            "Reinitializing... ")
        metadata <- list(genome=NA_character_,backend=NA_character_,
            filters=setNames(data.frame(matrix(ncol=4,nrow=0)),
            c("parameter","name","value","type","filtered")))
    }
        
    if (!("genome" %in% names(metadata))) {
        warning("Required 'genome' member was not found in metadata list and ",
            "was be added automatically.")
        metadata$genome <- NA_character_
    }
    if (!("backend" %in% names(metadata))) {
        warning("Required 'backend' member was not found in metadata list and ",
            "was be added automatically.")
        m$backend <- NA_character_
    }
    if (!("filters" %in% names(metadata))) {
        warning("Required 'filters' member was not found in metadata list and ",
            "was be added automatically.")
        metadata$filters <- .initFilterInfo()
    }
    
    se <- SummarizedExperiment(assays=genotypes,rowData=features,
        colData=samples,metadata=metadata)
    return(.GWASExperiment(se,phenotypes=phenotypes))
}

setGeneric("gfeatures",function(x,...) standardGeneric("gfeatures"))

setGeneric("gfeatures<-",function(x,...,value) standardGeneric("gfeatures<-"))

setGeneric("gsamples",function(x,...) standardGeneric("gsamples"))

setGeneric("gsamples<-",function(x,...,value) standardGeneric("gsamples<-"))

setGeneric("phenotypes",function(x,...) standardGeneric("phenotypes"))

setGeneric("phenotypes<-",function(x,...,value) standardGeneric("phenotypes<-"))

setGeneric("filters",function(x,...) standardGeneric("filters"))

setGeneric("filters<-",function(x,...,value) standardGeneric("filters<-"))

setMethod("gfeatures","GWASExperiment",function(x,withDimnames=TRUE) {
    return(rowData(x))
})

setReplaceMethod("gfeatures","GWASExperiment",function(x,...,value) {
    rowData(x) <- value
    return(x)
})

setMethod("gsamples","GWASExperiment",function(x,withDimnames=TRUE) {
    return(colData(x))
})

setReplaceMethod("gsamples","GWASExperiment",function(x,...,value) {
    if (!is(value,"DataFrame"))
        value <- DataFrame(value)
    colData(x) <- value
    return(x)
})

setMethod("phenotypes","GWASExperiment",function(x,withDimnames=TRUE) {
    p <- x@phenotypes
    if (!is.null(p) && ncol(p) > 0) {
        if (withDimnames)
            rownames(p) <- colnames(x)
        else
            rownames(p) <- NULL
    }
    return(p)
})

setReplaceMethod("phenotypes","GWASExperiment",function(x,...,value) {
    x@phenotypes <- value
    if (validObject(x))
        return(x)
})

setMethod("genome","GWASExperiment",function(x) {
    m <- metadata(x)
    return(m$genome)
})

# NOTE: If the 2nd argument is not 'value' but is 'y', it violates the generic
# which is defined with value as 2nd argument
setReplaceMethod("genome","GWASExperiment",function(x,value) {
    m <- metadata(x)
    m$genome <- value
    metadata(x) <- m
    if (validObject(x))
        return(x)
})

setMethod("filters","GWASExperiment",function(x) {
    m <- metadata(x)
    return(m$filters)
})

setReplaceMethod("filters","GWASExperiment",function(x,...,value) {
    m <- metadata(x)
    m$filters <- value
    metadata(x) <- m
    if (validObject(x))
        return(x)
})

setMethod("[","GWASExperiment",function(x,i,j,drop=TRUE) {
    p <- phenotypes(x,withDimnames=FALSE)
    
    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<",class(x),">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(i,
                rownames(x),fmt)
        }
        i <- as.vector(i)
    }

    if (!missing(j)) {
        if (is.character(j)) {
            fmt <- paste0("<",class(x),">[,j] index out of bounds: %s")
            j <- SummarizedExperiment:::.SummarizedExperiment.charbound(j,
                colnames(x),fmt)
        }
        j <- as.vector(j)
        p <- p[j,,drop=FALSE]
    }

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out,phenotypes=p,check=FALSE)
})

setReplaceMethod("[",c("GWASExperiment","ANY","ANY","GWASExperiment"),
    function(x,i,j,...,value) {
    p <- phenotypes(x,withDimnames=FALSE)
    
    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<",class(x),">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(i,
                rownames(x),fmt)
        }
        i <- as.vector(i)
    }

    if (!missing(j)) {
        if (is.character(j)) {
            fmt <- paste0("<",class(x),">[,j] index out of bounds: %s")
            j <- SummarizedExperiment:::.SummarizedExperiment.charbound(j,
                colnames(x),fmt)
        }
        j <- as.vector(j)
        p[j,] <- phenotypes(value,withDimname=FALSE)
    }

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out,phenotypes=p,check=FALSE)
})

# cbind means adding samples - phenotypes may not be identical
setMethod("cbind","GWASExperiment",function(...,deparse.level=1) {
    args <- list(...)
    p <- lapply(args,phenotypes,withDimnames=FALSE)
    p <- do.call(rbind,p)

    # Checks for identical column state.
    ref <- args[[1]]
    rp <- phenotypes(ref,withDimnames=FALSE)
    for (x in args[-1]) {
        if (!identical(colnames(rp),colnames(phenotypes(x,withDimnames=FALSE))))
            stop("Combining phenotype column names are not compatible")
    }

    oldValidity <- S4Vectors:::disableValidity()
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(oldValidity))

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out,phenotypes=p,check=FALSE)
})

# rbind means adding SNPs - phenotypes must be identical
setMethod("rbind","GWASExperiment",function(...,deparse.level=1) {
    args <- list(...)
    
    # Checks for identical column state.
    ref <- args[[1]]
    rp <- phenotypes(ref,withDimnames=FALSE)
    for (x in args[-1]) {
        if (!identical(rp,phenotypes(x,withDimnames=FALSE)))
            stop("Combining phenotypes must be identical")
    }

    oldValidity <- S4Vectors:::disableValidity()
    S4Vectors:::disableValidity(TRUE)
    on.exit(S4Vectors:::disableValidity(oldValidity))

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out,phenotypes=rp,check=FALSE)
})

setAs("SummarizedExperiment","GWASExperiment",function(from) {
    new("GWASExperiment",from, 
        phenotypes=DataFrame(),
        metadata=list(genome=NA_character_,backend=NA_character_,
        filters=setNames(data.frame(matrix(ncol=4,nrow=0)),
            c("name","value","type","filtered")),metadata(from)))
})

setValidity2("GWASExperiment",function(obj) {
    msg <- NULL
    
    # assays must be a transposed SnpMatrix
    if (!is(assays(obj)[[1]],"SnpMatrix") && !is(assays(obj)[[1]],"bigSNP"))
        msg <- c(msg,"The assays must be either SnpMatrix of bigSNP class!")
    
    # rowData must have the same #rows as rows of assays and same rownames
    if (nrow(assays(obj)[[1]]) != nrow(rowData(obj)))
        msg <- c(msg,paste0("Supporting feature (row) data must have the same ",
            "number of rows as the rows of assay data!"))
    no <- rownames(assays(obj)[[1]])
    nr <- rownames(rowData(obj))
    if (!is.null(no) && !is.null(nr)) {
        if (!identical(no,nr))
            msg <- c(msg,paste0("When provided, assay and feature data must ",
                "have identical rownames!"))
    }
    
    # colData must have the same #rows as columns of assays and colnames of
    # assays must be identical to rownames of colData
    if (ncol(assays(obj)[[1]]) != nrow(colData(obj)))
        msg <- c(msg,paste0("Supporting sample (column) data must have the ",
            "same number of rows as the columns of assay data!"))
    no <- colnames(assays(obj)[[1]])
    nr <- rownames(colData(obj))
    if (!is.null(no) && !is.null(nr)) {
        if (!identical(no,nr))
            msg <- c(msg,paste0("When provided, assay colnames and sample ",
                "rownames must be identical!"))
    }
    
    # If phenotypes are given, they must have same #rows as assay columns
    # and at least one column and colnames of assays must be identical to 
    # rownames of colData
    p <- phenotypes(obj)
    if (!is.null(p) && nrow(p) > 0) {
        if (nrow(p) != ncol(assays(obj)[[1]]) || nrow(p) < 1)
            msg <- c(msg,paste0("When provided, phenotype data must have the ",
                "same number of rows as the number of samples and at least ",
                "one phenotype!"))
    }
    no <- colnames(assays(obj)[[1]])
    nr <- rownames(p)
    if (!is.null(no) && !is.null(nr)) {
        if (!identical(no,nr))
            msg <- c(msg,paste0("When provided, assay colnames and phenotype ",
                "data rownames must be identical!"))
    }
    
    # It MUST have metadata with certain names in the list
    m <- metadata(obj)
    if (!is.character(m$genome))
        msg < c(msg,"The 'genome' member of metadata must be a character!")
    if (!is.character(m$backend))
        msg < c(msg,"The 'backend' member of metadata must be a character!")
    if (!.checkFilterInfo(m$filters))
        msg <- c(msg,paste0("The 'filters' member of metadata is not properly ",
            "formatted!"))

    if (length(msg) > 0)
        return(msg)
    else
        return(TRUE)
})

# Essentially the same as SummarizedExperiment but 
setMethod("show","GWASExperiment",function(object) {
    cat("class:",class(object),"\n")
    cat("dim:",dim(object),"\n")

    # metadata()
    expt <- names(metadata(object))
    if (is.null(expt))
        expt <- character(length(metadata(object)))
    coolcat("metadata(%d): %s\n",expt)

    # assays()
    # Assay type
    cls <- class(assays(object)[[1]])[1]
    nr <- nrow(assays(object)[[1]])
    nc <- ncol(assays(object)[[1]])
    cat("assays class:",cls,"with",nr,"rows (SNPs) and",nc,
        "columns (samples)\n")
    nms <- assayNames(object)
    if (is.null(nms))
        nms <- character(length(assays(object,withDimnames=FALSE)))
    coolcat("assays(%d): %s\n",nms)

    # rownames()
    rownames <- rownames(object)
    if (!is.null(rownames)) 
        coolcat("rownames(%d): %s\n",rownames)
    else 
        cat("rownames: NULL\n")

    # rowData()
    coolcat("feature names(%d): %s\n",names(rowData(object,use.names=FALSE)))

    # colnames()
    colnames <- colnames(object)
    if (!is.null(colnames)) 
        coolcat("colnames(%d): %s\n",colnames)
    else 
        cat("colnames: NULL\n")

    # colData()
    coolcat("sample names(%d): %s\n",names(colData(object)))
    
    # phenotypes()
    coolcat("phenotype names(%d): %s\n",names(phenotypes(object)))
})

# Constructor for SnpMatrix as it's missing from the snpStats package
SnpMatrix <- function(snp) {
    if (missing(snp))
        return(new("SnpMatrix",matrix(as.raw(0),0,0)))
    else
        return(new("SnpMatrix",snp))
}

# test: all(names(out) %in% ...)
.initFilterInfo <- function() {
    return(setNames(data.frame(matrix(ncol=5,nrow=0)),
        c("parameter","name","value","type","filtered")))
}

.checkFilterInfo <- function(d) {
    return(is.data.frame(d) && ncol(d) == 5
        && all(names(d) %in% c("parameter","name","value","type","filtered")))
}


#~ .GWASExperiment <- setClass("GWASExperiment",
#~  contains="SummarizedExperiment",
#~     slots=c(
#~         genotypes="SnpMatrix",
#~         fam="data.frame",
#~         map="data.frame",
#~         pheno="data.frame",
#~         backend="character",
#~         filters="list",
#~         genome="character_OR_NULL" # From S4Vectors
#~     ),
#~     prototype=(
#~         genotypes=matrix(as.raw(0),0,0),
#~         fam=data.frame(),
#~         map=data.frame(),
#~         pheno=data.frame(),
#~         backend=NA_character_,
#~         filters=setNames(data.frame(matrix(ncol=4,nrow=0)),
#~             c("name","value","type","filtered")),
#~         genome=NA_character_
#~     )
#~ )

#rowToColMat
# phenotypes is a table with as many rows as the samples and as many columns
# as the recorded phenotypes

#~ gfeatures <- function(x,y) {
#~     .checkX(x)
#~     if (!missing(y)) {
#~         y <- .checkY(y)
#~         rowData(x) <- y
#~         return(x)
#~     }
#~     return(rowData(x))
#~ }

#~ gsamples <- function(x,y) {
#~     .checkX(x)
#~     if (!missing(y)) {
#~         y <- .checkY(y)
#~         colData(x) <- y
#~         return(x)
#~     }
#~     return(colData(x))
#~ }

#~ .checkX <- function(x) {
#~     if (!is(x,"GWASExperiment"))
#~         stop("x must be a GWASExperiment object!")
#~ }

#~ .checkY <- function(y) {
#~     if (!is.data.frame(y) && !is(y,"DataFrame"))
#~         stop("y must be a data.frame or a DataFrame object!")
#~     if (is.data.frame(y))
#~         y <- DataFrame(y)
#~     return(y)
#~ }
