# There will probably be another class for storing retrieved from public APIs

setClassUnion("data.frame_OR_matrix",c("data.frame","matrix"))

setClassUnion("data.frame_OR_matrix_OR_NULL",c("data.frame","matrix","NULL"))

setClassUnion("data.frame_OR_DataFrame_OR_NULL",
    c("data.frame","DataFrame","NULL"))

# Definition
.GWASExperiment <- setClass("GWASExperiment",
    contains="SummarizedExperiment",
    slots=c(
        phenotypes="data.frame_OR_DataFrame_OR_NULL",
        pvalues="data.frame_OR_matrix_OR_NULL",
        effects="data.frame_OR_matrix_OR_NULL",
        prsnps="data.frame_OR_matrix_OR_NULL"
    )
)

# Constructor
GWASExperiment <- function(genotypes=SnpMatrix(),features=DataFrame(),
    samples=DataFrame(),phenotypes=DataFrame(),pvalues=NULL,effects=NULL,
    prsnps=NULL,metadata=list(genome=NA_character_,backend=NA_character_,
        filters=setNames(data.frame(matrix(ncol=5,nrow=0)),
            c("parameter","name","value","type","filtered"))),...) {
    
    # It MUST have metadata with certain names in the list
    if (!is.list(metadata)) {
        warning("The provided metadata argument is not a list! ",
            "Reinitializing... ")
        metadata <- list(genome=NA_character_,backend=NA_character_,
            filters=setNames(data.frame(matrix(ncol=5,nrow=0)),
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
        metadata$backend <- NA_character_
    }
    if (!("filters" %in% names(metadata))) {
        warning("Required 'filters' member was not found in metadata list and ",
            "was be added automatically.")
        metadata$filters <- .initFilterInfo()
    }
    
    se <- SummarizedExperiment(assays=genotypes,rowData=features,
        colData=samples,metadata=metadata)
    return(.GWASExperiment(se,phenotypes=phenotypes,pvalues=pvalues,
        effects=effects,prsnps=prsnps))
}

setGeneric("genotypes",function(x,...) standardGeneric("genotypes"))

setGeneric("genotypes<-",function(x,...,value) standardGeneric("genotypes<-"))

setGeneric("gfeatures",function(x,...) standardGeneric("gfeatures"))

setGeneric("gfeatures<-",function(x,...,value) standardGeneric("gfeatures<-"))

setGeneric("gsamples",function(x,...) standardGeneric("gsamples"))

setGeneric("gsamples<-",function(x,...,value) standardGeneric("gsamples<-"))

setGeneric("phenotypes",function(x,...) standardGeneric("phenotypes"))

setGeneric("phenotypes<-",function(x,...,value) standardGeneric("phenotypes<-"))

setGeneric("pvalues",function(x,...) standardGeneric("pvalues"))

setGeneric("pvalues<-",function(x,...,value) standardGeneric("pvalues<-"))

setGeneric("effects",function(x,...) standardGeneric("effects"))

setGeneric("effects<-",function(x,...,value) standardGeneric("effects<-"))

setGeneric("prsnps",function(x,...) standardGeneric("prsnps"))

setGeneric("prsnps<-",function(x,...,value) standardGeneric("prsnps<-"))

setGeneric("filterRecord",function(x,...) standardGeneric("filterRecord"))

setGeneric("filterRecord<-",function(x,...,value) standardGeneric("filterRecord<-"))

setGeneric("pcaCovariates",function(x,...) standardGeneric("pcaCovariates"))

setGeneric("pcaCovariates<-",
    function(x,...,value) standardGeneric("pcaCovariates<-"))

setMethod("genotypes","GWASExperiment",function(x,i=1) {
    return(assay(x,i))
})

setReplaceMethod("genotypes","GWASExperiment",function(x,i=1,...,value) {
    assay(x,i) <- value
    return(x)
})

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

setMethod("pvalues","GWASExperiment",function(x,withDimnames=TRUE) {
    p <- x@pvalues
    if (!is.null(p)) {
        if (withDimnames)
            rownames(p) <- rownames(x)
        else
            rownames(p) <- NULL
    }
    return(p)
})

setReplaceMethod("pvalues","GWASExperiment",function(x,...,value) {
    x@pvalues <- value
    if (validObject(x))
        return(x)
})

setMethod("effects","GWASExperiment",function(x,withDimnames=TRUE) {
    e <- x@effects
    if (!is.null(e)) {
        if (withDimnames)
            rownames(e) <- rownames(x)
        else
            rownames(e) <- NULL
    }
    return(e)
})

setReplaceMethod("effects","GWASExperiment",function(x,...,value) {
    x@effects <- value
    if (validObject(x))
        return(x)
})

setMethod("prsnps","GWASExperiment",function(x,withDimnames=TRUE) {
    s <- x@prsnps
    if (!is.null(s)) {
        if (withDimnames)
            rownames(s) <- rownames(x)
        else
            rownames(s) <- NULL
    }
    return(s)
})

setReplaceMethod("prsnps","GWASExperiment",function(x,...,value) {
    x@prsnps <- value
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

setMethod("filterRecord","GWASExperiment",function(x) {
    m <- metadata(x)
    return(m$filters)
})

setReplaceMethod("filterRecord","GWASExperiment",function(x,...,value) {
    m <- metadata(x)
    m$filters <- value
    metadata(x) <- m
    if (validObject(x))
        return(x)
})

setMethod("pcaCovariates","GWASExperiment",function(x) {
    m <- metadata(x)
    return(m$PCs)
})

setReplaceMethod("pcaCovariates","GWASExperiment",function(x,...,value) {
    m <- metadata(x)
    # Validation of value, e.g. data.frame, must have same number of samples as
    # the object etc.
    m$PCs <- value
    metadata(x) <- m
    if (validObject(x))
        return(x)
})

setMethod("[","GWASExperiment",function(x,i,j,drop=TRUE) {
    p <- phenotypes(x,withDimnames=FALSE)
    v <- pvalues(x,withDimnames=FALSE)
    e <- effects(x,withDimnames=FALSE)
    s <- prsnps(x,withDimnames=FALSE)
    
    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<",class(x),">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(i,
                rownames(x),fmt)
        }
        i <- as.vector(i)
        v <- v[i,,drop=FALSE]
        e <- e[i,,drop=FALSE]
        s <- s[i,,drop=FALSE]
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
    
    # If i,j is not missing, i.e. we are subsetting, LDsnps, pcaCov and PCs are
    # no longer valid, so must be dropped to avoid accidental usage. We keep
    # pcaRob for visualization purposes.
    # There must be an option to recalculate, set with R options. We keep 
    # filters for now. Whether deleted or not, will be controlled by various
    # levels of this R option.
    if (!missing(i) || !missing(j)) {
        m <- metadata(x)
        m$LDsnps <- m$pcaCov <- m$PCs <- NULL
        metadata(x) <- m
    }

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out,phenotypes=p,pvalues=v,effects=e,prsnps=s,
        check=FALSE)
})

setReplaceMethod("[",c("GWASExperiment","ANY","ANY","GWASExperiment"),
    function(x,i,j,...,value) {
    p <- phenotypes(x,withDimnames=FALSE)
    v <- pvalues(x,withDimnames=FALSE)
    e <- effects(x,withDimnames=FALSE)
    s <- prsnps(x,withDimnames=FALSE)
    
    if (!missing(i)) {
        if (is.character(i)) {
            fmt <- paste0("<",class(x),">[i,] index out of bounds: %s")
            i <- SummarizedExperiment:::.SummarizedExperiment.charbound(i,
                rownames(x),fmt)
        }
        i <- as.vector(i)
        v[i,] <- pvalues(value,withDimnames=FALSE)
        e[i,] <- effects(value,withDimnames=FALSE)
        s[i,] <- prsnps(value,withDimnames=FALSE)
    }

    if (!missing(j)) {
        if (is.character(j)) {
            fmt <- paste0("<",class(x),">[,j] index out of bounds: %s")
            j <- SummarizedExperiment:::.SummarizedExperiment.charbound(j,
                colnames(x),fmt)
        }
        j <- as.vector(j)
        p[j,] <- phenotypes(value,withDimnames=FALSE)
    }
    
    # Similarly as above, if i or j not missing, values have been replaced and
    # some metadata values are no longer valid.
    if (!missing(i) || !missing(j)) {
        m <- metadata(x)
        m$LDsnps <- m$pcaCov <- m$PCs <- NULL
        metadata(x) <- m
    }

    out <- callNextMethod()
    BiocGenerics:::replaceSlots(out,phenotypes=p,pvalues=v,effects=e,prsnps=s,
        check=FALSE)
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
    
    # Drop most metadata    
    m <- metadata(out)
    m$LDsnps <- m$pcaRob <- m$pcaCov <- m$PCs <- m$filters <- NULL
    metadata(out) <- m
    
    # pvalues must also be dropped as the population changes
    if (!is.null(pvalues(ref)))
        warning("Previous association tests will be dropped during ",
            "cbind operation as the population changes!")
    
    BiocGenerics:::replaceSlots(out,phenotypes=p,metadata=m,pvalues=NULL,
        effects=NULL,prsnps=NULL,check=FALSE)
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
    
    # Drop most metadata    
    m <- metadata(out)
    m$LDsnps <- m$pcaRob <- m$pcaCov <- m$PCs <- m$filters <- NULL
    metadata(out) <- m
    
    # pvalues must also be dropped as the population changes
    if (!is.null(pvalues(ref)))
        warning("Previous association tests will be dropped during ",
            "rbind operation as the number of markers changes!")
    
    BiocGenerics:::replaceSlots(out,phenotypes=rp,metadata=m,pvalues=NULL,
        effects=NULL,prsnps=NULL,check=FALSE)
})

setAs("SummarizedExperiment","GWASExperiment",function(from) {
    new("GWASExperiment",from, 
        phenotypes=DataFrame(),
        metadata=list(genome=NA_character_,backend=NA_character_,
        filters=setNames(data.frame(matrix(ncol=4,nrow=0)),
            c("name","value","type","filtered")),metadata(from)))
})

.checkGenotypesClass <- function(obj,msg) {
    if (!is(assays(obj)[[1]],"SnpMatrix") && !is(assays(obj)[[1]],"bigSNP"))
        msg <- c(msg,"The assays must be either SnpMatrix of bigSNP class!")
    return(msg)
}

.checkFeaturesValidity <- function(obj,msg) {
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
    return(msg)
}

.checkSamplesValidity <- function(obj,msg) {
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
    return(msg)
}

.checkPhenotypesValidity <- function(obj,msg) {
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
    return(msg)
}

.checkPvaluesValidity <- function(obj,msg) {
    v <- pvalues(obj)
    if (!is.null(v) && nrow(v) > 0) {
        if (nrow(v) != nrow(assays(obj)[[1]]) || ncol(v) < 1)
            msg <- c(msg,paste0("When provided, association pvalues must have ",
                "the same number of rows as the number of markers and at ",
                "least one test performed (1 column)!"))
    }
    no <- rownames(assays(obj)[[1]])
    nr <- rownames(v)
    if (!is.null(no) && !is.null(nr)) {
        if (!identical(no,nr))
            msg <- c(msg,paste0("When provided, assay rownames and associated ",
                "p-value rownames must be identical!"))
    }
    return(msg)
}

.checkEffectsValidity <- function(obj,msg) {
    e <- effects(obj)
    if (!is.null(e) && nrow(e) > 0) {
        if (nrow(e) != nrow(assays(obj)[[1]]) || ncol(e) < 1)
            msg <- c(msg,paste0("When provided, association effects must have ",
                "the same number of rows as the number of markers and at ",
                "least one test performed (1 column)!"))
    }
    no <- rownames(assays(obj)[[1]])
    nr <- rownames(e)
    if (!is.null(no) && !is.null(nr)) {
        if (!identical(no,nr))
            msg <- c(msg,paste0("When provided, assay rownames and associated ",
                "effect rownames must be identical!"))
    }
    return(msg)
}

.checkPrsnpsValidity <- function(obj,msg) {
    s <- prsnps(obj)
    if (!is.logical(s))
        msg <- c(msg,paste0("When provided, SNPs in a PRS must be a matrix of ",
                "logical values (TRUE/FALSE)!"))
    if (!is.null(s) && nrow(s) > 0) {
        if (nrow(s) != nrow(assays(obj)[[1]]) || ncol(s) < 1)
            msg <- c(msg,paste0("When provided, SNPs in a PRS must have the ",
                "same number of rows as the number of markers and at least ",
                "one PRS algorithm associated (1 column)!"))
    }
    no <- rownames(assays(obj)[[1]])
    nr <- rownames(s)
    if (!is.null(no) && !is.null(nr)) {
        if (!identical(no,nr))
            msg <- c(msg,paste0("When provided, assay rownames and associated ",
                "polygenic score SNP rownames must be identical!"))
    }
    return(msg)
}

.checkMetadataValidity <- function(obj,msg) {
    m <- metadata(obj)
    if (!is.character(m$genome))
        msg < c(msg,"The 'genome' member of metadata must be a character!")
    if (!is.character(m$backend))
        msg < c(msg,"The 'backend' member of metadata must be a character!")
    if (!.checkFilterInfo(m$filters))
        msg <- c(msg,paste0("The 'filters' member of metadata is not properly ",
            "formatted!"))
    if (!is.null(m$PCs)) {
        if (!is(m$PCs,"data.frame_OR_matrix"))
            msg <- c(msg,paste0("If provided, the 'PCs' member of metadata ",
                "must be a data.frame or a matrix!"))
        if (is(m$PCs,"data.frame_OR_matrix")) {
            if (!all(rownames(m$PCs) %in% colnames(assays(obj)[[1]])))
                msg <- c(msg,paste0("If provided, the 'PCs' member of ",
                    "metadata must have the same rownames as the column names ",
                    "of the GWASExperiment object!"))
        }
    }
    return(msg)
}

setValidity2("GWASExperiment",function(obj) {
    msg <- NULL
    
    # assays must be a transposed SnpMatrix or bigSNP
    msg <- .checkGenotypesClass(obj,msg)
    
    # rowData must have the same #rows as rows of assays and same rownames
    msg <- .checkFeaturesValidity(obj,msg)
    
    # colData must have the same #rows as columns of assays and colnames of
    # assays must be identical to rownames of colData
    msg <- .checkSamplesValidity(obj,msg)
    
    # If phenotypes are given, they must have same #rows as assay columns
    # and at least one column and colnames of assays must be identical to 
    # rownames of colData
    msg <- .checkPhenotypesValidity(obj,msg)
    
    # If pvalues are given, they must have same #rows as assay rows
    # and at least one column and must be a matrix-like object
    msg <- .checkPvaluesValidity(obj,msg)
    
    # If effects are given, they must have same #rows as assay rows
    # and at least one column and must be a matrix-like object
    msg <- .checkEffectsValidity(obj,msg)
    
    # If prsnps are given, they must have same #rows as assay rows
    # and at least one column and must be a matrix-like object of logicals
    msg <- .checkPrsnpsValidity(obj,msg)
    
    # It MUST have metadata with certain names in the list
    msg <- .checkMetadataValidity(obj,msg)
    
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
    
    # pvalues()
    coolcat("performed tests(%d): %s\n",colnames(pvalues(object)))
    
    # prsnps()
    coolcat("associated PRS(%d): %s\n",colnames(prsnps(object)))
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

.hasPcaCovariates <- function(x) {
    if (!is(x,"GWASExperiment"))
        stop("Input must be an object of class GWASExperiment!")
    return(!.isEmpty(pcaCovariates(x)))
}

#setAs("GWASExperiment","gData",function(from) {
#    if (!requireNamespace("statgenGWAS"))
#        stop("R package statgenGWAS is required!")
#    
#    # Genotypes
#    geno <- as(t(assay(from,1)),"numeric")
#    
#    # Map
#    map <- as.data.frame(gfeatures(from))
#    map <- map[,c("chromosome","position")]
#    names(map) <- c("chr","pos")
#    
#    # Phenotypes - drop those which are not numeric for now...
#    pheno <- phenotypes(from)
#    chn <- vapply(names(pheno),function(x,p) {
#        return(is.character(p[[x]]))
#    },logical(1),pheno)
#    if (!all(chn)) {
#        warning("Some phenotypes are not numeric (",
#            paste0(names(pheno)[which(chn)],collapse=", "),")!\nThese will ",
#            "not be passed to the resulting gData object.",immediate.=TRUE)
#        pheno <- pheno[,!chn,drop=FALSE]
#    }
#    pheno <- data.frame(genotype=rownames(pheno),pheno)
#    
#    return(createGData(geno=geno,map=map,pheno=pheno))
#})

# Convert GWASExperiment to gData (statgenGWAS)
GWASExperiment2gData <- function(obj,covariates=NULL,pcs=FALSE) {
    if (!requireNamespace("statgenGWAS"))
        stop("R package statgenGWAS is required!")
    
    # Phenotypes - drop those which are not numeric for now... and check with
    # covariates... If not everything satisfied, we do not proceed.
    pheno <- phenotypes(obj)
    chn <- vapply(names(pheno),function(x,p) {
        return(is.character(p[[x]]))
    },logical(1),pheno)
    if (any(chn)) {
        warning("Some phenotypes are not numeric (",
            paste0(names(pheno)[which(chn)],collapse=", "),")!\nThese will ",
            "not be passed to the resulting gData object.",immediate.=TRUE)
        pheno <- pheno[,!chn,drop=FALSE]
    }
    
    if (pcs) { # Include precalculated PCs if any
        if (.hasPcaCovariates(obj)) {
            pcov <- pcaCovariates(obj)
            pheno <- cbind(pheno,pcov)
            covariates <- c(covariates,colnames(pcov))
        }
        else
            warning("PC covariates reqeusted in the model, but no calculated ",
                "PC covariates found! Ignoring...",immediate.=TRUE)
    }
    
    covar <- NULL
    if (!is.null(covariates)) {
        if (!(all(covariates %in% names(pheno)))) {
            no <- which(!(covariates %in% names(pheno)))
            stop("Covariate(s) ",paste0(covariates[no],collapse=", "),
                " not found in gData phenotypes!\nHave they been dropped ",
                "because they were not numeric?")
        }
        
        covar <- pheno[,covariates,drop=FALSE]
        pheno <- pheno[,-match(covariates,names(pheno)),drop=FALSE]
    }
    pheno <- data.frame(genotype=rownames(pheno),pheno)
    
    
    # Genotypes
    geno <- as(t(assay(obj,1)),"numeric")
    
    # Map
    map <- as.data.frame(gfeatures(obj))
    map <- map[,c("chromosome","position")]
    names(map) <- c("chr","pos")
    
    return(createGData(geno=geno,map=map,pheno=pheno,covar=covar))
}
