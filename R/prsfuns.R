#! Add betas accessor to GWASExperiment (maybe after PRSice)

# More arguments to come, like PRS formula like in PRSice etc.
PRS <- function(obj,snps,response,covariates=NULL,pcs=FALSE,...) {
    # Some validation here
    
    X <- t(as(genotypes(obj),"numeric"))
    prs <- X[,snps$variant_id] %*% snps[snps$variant_id,"effect_weight"]
}

lassosumPRS <- function(obj,response,covariates=NULL,pcs=FALSE,
    lsWspace=NULL,anc=c("EUR","ASN","AFR"),valid=c("auto","stack","split")) {
    # Can run PRS? Has necessary fields?
    .canRunPrs(obj,response)
    # Workspace valid?
    lsWspace <- .validateWorkspacePath(lsWspace,"lassosum")
    # Ancestry for LD blocks based on 1000G provided by lassosum
    anc <- toupper(anc[1])
    valid <- valid[1]
    .checkTextArgs("Ancestry (anc)",anc,c("EUR","ASN","AFR"),multiarg=FALSE)
    .checkTextArgs("Validation type (valid)",valid,c("auto","stack","split"),
        multiarg=FALSE)
    
    if (valid == "auto")
        valid <- ifelse(ncol(obj) >= 1000,"split","stack")
    
    # If all ok, prepare the run
    prepList <- .prepareLassosumRun(obj,response,lsWspace,anc)
    
    # Prepare the phenotypes for later prediction/validation
    p <- phenotypes(obj)
    chResCov <- .validateResponseAndCovariates(p,response,covariates)
    response <- chResCov$res
    covariates <- chResCov$cvs
    pheno <- p[,c("FID","IID",response)]
    covar <- p[,c("FID","IID",covariates)]
    # Make sure FID and IID matches, otherwise cannot be properly handled by
    # snpStats::write.plint
    if (!identical(pheno[,"FID"],pheno[,"IID"])) {
        pheno[,"IID"] <- pheno[,"FID"]
        covar[,"IID"] <- covar[,"FID"]
    }
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
    
    # We need to attach marker ids to ss for later extraction
    marks <- rownames(ss)
    names(marks) <- paste0(ss$chromosome,"_",ss$position,"_",ss$allele.1,"_",
        ss$allele.2)
    
    # Run the lassosum pipeline
    disp("Running lassosum PRS algorithm")
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
    # Adjust the output - actual location of temporary PLINK file
    out$test.bfile <- file.path(lsWspace,prepList$bfile)
    # Marker names
    tmp <- out$sumstats
    rownames(tmp) <- paste0(tmp$chr,"_",tmp$pos,"_",tmp$A1,"_",tmp$A2)
    rownames(tmp) <- marks[rownames(tmp)]
    out$sumstats <- tmp
    setwd(cwd)

    # PRS? If sample size is >1000 then splitvalidate else validate
    disp("Running lassosum validation and PRS calculation")
    if (valid == "stack")
        val <- lassosum::validate(out,pheno=pheno,covar=covar,plot=FALSE)
    else if (valid == "split")
        val <- lassosum::splitvalidate(out,pheno=pheno,covar=covar,plot=FALSE)
    
    disp("Done! R^2 is: ",max(targetRes$validation.table$value)^2)
    return(.extractLassosumPrsComponents(out,val,genome(obj)))
    # or
    # betas(obj) <- .extractLassosumPrsComponents(out,val,genome(obj),full=T)
    # return(obj)
}

.extractLassosumPrsComponents <- function(out,val,gb) {
    # Construct initial data frame (we have matched marker names with locations)
    df <- out$sumstats[,c("chr","pos","A1","A2")]
    colnames(df) <- c("chromosome","position","risk_allele","reference_allele")
    df$variant_id <- rownames(df)
    
    # Assuming that betas are aligned with sumstats (they are)
    nonZero <- val$best.beta!=0
    df <- df[nonZero,,drop=FALSE]
    # Until we get an answer from lassosum authors, we reverse betas so as to be
    # concordant with the classical X %*% betas PRS definition
    df$effect_weight <- -val$best.beta[nonZero]
    # and conversion to odds ratios
    df$OR <- exp(df$effect_weight)
    df$locus_name <- rep(NA,nrow(df))
    
    # Final alignment with the external API fetch outcomes from PGS catalog
    df <- df[,c("chromosome","position","variant_id","risk_allele",
        "reference_allele","locus_name","effect_weight","OR")]
    if (is.null(gb))
        gb <- "nr"
    df$asm <- rep(gb,nrow(df))
    
    return(df)
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
