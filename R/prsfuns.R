lassosumPRS <- function(obj,response,covariates=NULL,pcs=FALSE,lsWspace=NULL,
    anc=c("EUR","ASN","AFR")) {
    # Can run PRS? Has necessary fields?
    .canRunPrs(obj,response)
    # Workspace valid?
    lsWspace <- .validateWorkspacePath(lsWspace,"lassosum")
    # Ancestry for LD blocks based on 1000G provided by lassosum
    anc <- toupper(anc[1])
    .checkTextArgs("Ancestry (anc)",anc,c("EUR","ASN","AFR"),multiarg=FALSE)
    
    # If all ok, prepare the run
    prepList <- .prepareLassosumRun(obj,response,lsWspace,anc)
    
    # Prepare the phenotypes for later prediction/validation
    p <- phenotypes(obj)
    chResCov <- .validateResponseAndCovariates(p,response,covariates)
    response <- chResCov$res
    covariates <- chResCov$cvs
    pheno <- p[,c("FID","IID",response)]
    covar <- p[,c("FID","IID",covariates)]
    if (pcs) { # Include robust PCs in the model
        if (.hasPcaCovariates(obj)) {
            pcov <- pcaCovariates(obj)
            covar <- cbind(covar,pcov)
        }
        else
            warning("PC covariates requested in the model, but not calculated ",
                "PC covariates found! Ignoring...",immediate.=TRUE)
    }
    
    # Transform the P-values into correlation
    ss <- prepList$ss
    cor <- p2cor(p=ss$pvalue,n=ncol(obj),sign=ss$effect)
    
    # Run the lassosum pipeline
    cwd <- getwd()
    setwd(lsWspace)
    out <- lassosum.pipeline(
        cor=cor,
        chr=ss$chromosome,
        pos=ss$position,
        A1=ss$allele.1,
        A2=ss$allele.2,
        ref.bfile=prepList$bfile,
        test.bfile=prepList$bfile,
        LDblocks=paste0(anc,".",prepList$gb)
    )
    setwd(cwd)
    
    # PRS?
    targetRes <- validate(out,pheno=pheno,covar=covar)
    # Get the maximum R2
    r2 <- max(target.res$validation.table$value)^2
}

.prepareLassosumRun <- function(obj,pheno,wspace,anc) {
    disp("Preparing lassosum run at workspace directory: ",wspace)
    
    # Prepare the sum stat
    disp("  preparing summary statistics...")
    statsCoords <- gfeatures(obj)
    statsCoords <- as.data.frame(statsCoords[,c("chromosome","position",
        "allele.1","allele.2")])
    p <- pvalues(obj)[,ncol(pvalues(obj))]
    e <- effects(obj)
    
    # Until we find a better way of summarizing effects across methods,
    # prioritize effect selection according to algorithm popularity
    pri <- .getGwaLinArgPrior()
    if (any(pri %in% colnames(e)))
        oe <- e[,pri[which(pri %in% colnames(e))[1]]]
    sumStat <- data.frame(statsCoords,effect=oe,pvalue=p)

    # Prepare PLINK files
    disp("  preparing PLINK temporary files...")
    bfileBase <- paste0("plink_lassosum_",.randomString())
    writePlink(obj,pheno,outBase=file.path(wspace,bfileBase))
    
    # Copy LD block files to workspace directory
    disp("  copying linkage disequilibrium block files...")
    ldHome <- system.file("data",package="lassosum")
    gb <- ifelse(is.null(genome(obj)),"hg19",genome(obj))
    fromLd <- file.path(ldHome,paste0("Berisa.",anc,".",gb,".bed"))
    toLd <- file.path(wspace,paste0("Berisa.",anc,".",gb,".bed"))
    file.copy(from=fromLd,to=toLd,overwrite=TRUE)
    
    disp("Preparation done!")
    return(list(ss=sumStat,bfile=bfileBase,ldb=toLd,gb=gb))
}

.validateWorkspacePath <- function(path,tool="tool") {
    # Various path options
    if (!is.null(path)) {
        if (!is.character(path)) {
            warning("The desired ",tool," workspace must be a valid string! ",
                "Setting to NULL for autocreation",immediate.=TRUE)
            path <- NULL
        }
        else {
            if (!dir.exists(path)) {
                # No disp here
                message("The desired ",tool," workspace does not exist!")
                confirm <- menu(c("Yes","No"),
                    title=sprintf(paste0(tool," workspace is going to be ",
                        "created at %s. Do you agree?"),path))
                if (confirm == 1)
                    dir.create(path,recursive=TRUE)
                else
                    stop(paste0(tool," cannot be executed without a valid ",
                        "workspace.\nPlease agree or use the 'path' argument ",
                        "in the respective function."))
            }
        }
    }
    else {
        path <- file.path(tempdir(),paste0("lassosum_",.randomString()))
        dir.create(path,recursive=TRUE)
    }
    
    # Finally
    disp(tool," operations will be executed at ",path)
    return(path)
}


.getGwaLinArgPrior <- function() {
    return(c("snptest","statgen","glm","rrblup"))
}

.canRunPrs <- function(o,p) {
    # Object class?
    if (!is(o,"GWASExperiment"))
        stop("PRS input must be an object of class GWASExperiment!")
    # Has an association test been performed?
    if (is.null(pvalues(o)) || is.null(effects(o)))
        stop("At least one association test must be performed prior to ",
            "building PRS with lassosum!")
    # The given phenotype exists?
    ph <- phenotypes(o)
    if (is.null(ph))
        stop("PRS cannot be performed without phenotypes!")
    else {
        if (is.numeric(p)) {
            if (p > ncol(ph))
                stop("The phenotype index ",p," is larger than the object's ",
                    "available phenotypes!")
        }
        else if (is.character(p)) {
            if (!(p %in% colnames(ph)))
                stop("The phenotype ",p," cannot be found in the object's ",
                    "available phenotypes (",paste(colnames(ph),collapse=", "),
                    ")!")
        }
    }
}
