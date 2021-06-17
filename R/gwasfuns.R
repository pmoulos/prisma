gwaGlm <- function(obj,response,covariates,pcs=TRUE,family=NULL,psig=0.05,
    rc=NULL,...) {
    if (!is(obj,"GWASExperiment"))
        stop("obj must be an object of class GWASExperiment!")
    
    .checkNumArgs("Association p-value",psig,"numeric",c(0,1),"both")
    
    # Construct the glm data frame, check if response and covariate names are
    # valid
    p <- phenotypes(obj)
    chResCov <- .validateResponseAndCovariates(p,response,covariates)
    response <- chResCov$res
    covariates <- chResCov$cvs
    
    # Reduce the phenotypes for testing and check GLM call integrity
    p <- p[,c(response,covariates)]
    #glmArgs <- list(...)
    # We must ensure that if binomial requested, res is 0-1 and binary
    #if ("family" %in% names(glmArgs) && glmArgs[["family"]] == "binomial")
    #    p[,response] <- .validateBinaryForBinomial(p[,response])
    if (!is.null(family)) {
        family <- family[1]
        .checkTextArgs("Regression family",family,
            c("gaussian","binomial","poisson"),multiarg=FALSE)
        if (family == "binomial")
            p[,response] <- .validateBinaryForBinomial(p[,response])
        fpredHelper <- "provided"
    }
    else {
        # If family not provided, we must decide whether to perform logistic
        # of gaussian regression. One way, tryCatch .validateBinaryForBinomial.
        # If fail, then gaussian else binomial
        fam <- NULL
        fam <- tryCatch(.validateBinaryForBinomial(p[,response]),
            error=function(e) { return("gaussian") },finally="")
        family <- ifelse(fam!="gaussian","binomial","gaussian")
        fpredHelper <- "predicted"
    }
    
    if (pcs) { # Include robust PCs in the model
        if (.hasPcaCovariates(obj)) {
            pcov <- pcaCovariates(obj)
            p <- cbind(p,pcov)
        }
        else
            warning("PC covariates reqeusted in the model, but no calculated ",
                "PC covariates found! Ignoring...",immediate.=TRUE)
    }
    
    # We now need to split the SnpMatrix, as glm is quite fast, parallelizng
    # per snp will most probably cause more overhead than chunking and iterating
    # We split according to the number of available cores for analysis.
    snps <- assay(obj,1)
    splitFactor <- .splitFactorForParallel(nrow(snps),rc)
    # We now need to split the SnpMatrix, as glm is quite fast, parallelizng
    # per snp will most probably cause more overhead than chunking and iterating
    # We split according to the number of available cores for analysis.
    #splits <- split(seq_len(nrow(snps)),splitFactor)
    splits <- split(rownames(snps),splitFactor)
    
    disp("Performing GWAS with GLM over ",length(splits)," chunks")
    disp("Trait(s)    : ",response)
    disp("Covariate(s): ",paste0(covariates,collapse=", "))
    disp("Regression  : ",family," (",fpredHelper,")")
    disp("Use PCs     : ",ifelse(pcs,"Yes","No"))
    assocList <- cmclapply(splits,function(n,p,r,f,...) {
        disp("  testing SNPs from ",n[1]," to ",n[length(n)])
        batch <- lapply(n,function(s,p,r,f,...) {
            disp("    testing ",s,level="full")
            dat <- cbind(p,t(as(snps[s,],"numeric")))
            tmp <- coefficients(.gwaGlmWorker(dat,r,f,...))
            rownames(tmp)[nrow(tmp)] <- s
            return(tmp[s,,drop=FALSE])
        },p,r,f,...)
        batch <- do.call("rbind",batch)
        rownames(batch) <- n
        return(batch)
    },p,response,family,...,rc=rc)
    
    #disp("Performing GWAS with GLM over ",length(splits)," chunks")
    #disp("Trait(s)    : ",response)
    #disp("Covariate(s): ",paste0(covariates,collapse=", "))
    #disp("Use PCs     : ",ifelse(pcs,"Yes","No"))
    #assocList <- cmclapply(splits,function(n,p,r,g) {
    #    disp("  testing SNPs from ",n[1]," to ",n[length(n)])
    #    batch <- lapply(n,function(s,p,r,g) {
    #        disp("    testing ",s,level="full")
    #        dat <- cbind(p,t(as(snps[s,],"numeric")))
    #        tmp <- coefficients(do.call(.gwaGlmWorker,dat,r,glmArgs))
    #        rownames(tmp)[nrow(tmp)] <- s
    #        return(tmp[s,,drop=FALSE])
    #    },p,r,g)
    #    batch <- do.call("rbind",batch)
    #    rownames(batch) <- n
    #    return(batch)
    #},p,response,glmArgs,rc=rc)
    
    assocResult <- do.call("rbind",assocList)
    disp("Found ",length(which(assocResult[,4]<psig))," associations ",
        "(uncorrected p-values).")
    
    return(assocResult)
}

# ... is other options passed to GLM
.gwaGlmWorker <- function(dat,res,fam,...) {
    # We assume that dat contains only the deserved covariates
    ii <- which(colnames(dat)==res)
    # Before constructing formula, all the covariates must have proper names
    colnames(dat) <- make.names(colnames(dat))
    covs <- colnames(dat)[-ii]
    cres <- colnames(dat)[ii]
    f <- as.formula(paste(cres,paste0(covs,collapse="+"),sep="~"))
    return(summary(glm(f,data=dat,family=fam,...)))
}

gwaBlup <- function(obj,response,covariates,usepc=c("auto","estim","fixed",
    "none"),npcs=NULL,psig=0.05,rc=NULL) {
    if (!is(obj,"GWASExperiment"))
        stop("obj must be an object of class GWASExperiment!")
        
    .checkNumArgs("Association p-value",psig,"numeric",c(0,1),"both")
    usepc <- usepc[1]    
    .checkTextArgs("Use PCA",usepc,c("auto","estim","fixed","none"),
        multiarg=FALSE)
    if (!is.null(npcs)) {
        npcs <- npcs[1]
        .checkNumArgs("Number of PCs",npcs,"numeric",0,"gte")
    }
    else {
        if (usepc == "fixed")
            warning("The number of PCs to use must be provided when ",
                "usepc = 'fixed'! Switching to usepc = 'estim'...",
                immediate.=TRUE)
        usepc <- "estim"
    }
    
    # Preprocess phenotypes, similarly to gwaGlm
    p <- phenotypes(obj)
    chResCov <- .validateResponseAndCovariates(p,response,covariates)
    response <- chResCov$res
    covariates <- chResCov$cvs
    p <- p[,c(response,covariates)]
    
    later <- FALSE
    if (usepc == "auto") {
        # If pcaCovariates exist, use as many PCs as in there (e.g. from TW)
        # otherwise, let rrBLUP estimate
        if (.hasPcaCovariates(obj)) {
            pcov <- pcaCovariates(obj)
            npcs <- ncol(pcov)
        }
        else
            usepc <- "estim"
    }
    if (usepc == "estim")
        later <- TRUE
    if (usepc == "fixed")
        # Silence, will use provided npcs which have been checked before
    if (usepc == "none")
        npcs <- 0
    
    # Add rownames at first column so as to be compatible with rrBLUP::GWAS
    pheno <- cbind(rownames(p),p)
    names(pheno)[1] <- "gid"
    
    # Prepare the genotypes
    disp("Preparing genotypes for rrBLUP...")
    geno <- .prepareGenotypesForBlup(obj)
    
    # Time to decide on number of PCs
    if (later) {
        disp("Estimating number of PCs to include in rrBLUP model...")
        npcs <- .estimateNPCinBlup(pheno,geno,covariates)
        disp("Finished! ",npcs," PCs will be included in rrBLUP model...")
    }

    # Already filtered MAF, so min.MAF=.Machine$double.eps
    # Also, GWAS does not like robust PCA... Either we chose to estimate the
    # best n.PC based on some kind of "elbow" in discovered SNPs, or based on
    # the number of PCs returned by Tracy-Widom test (ncol(pcov))
    disp("Performing GWAS with rrBLUP model")
    disp("Trait(s)    : ",response)
    disp("Covariate(s): ",paste0(covariates,collapse=", "))
    disp("Use PCs     : ",ifelse(npcs==0,"No","Yes"))
    log <- capture.output({
        rrb <- GWAS(pheno,geno,fixed=covariates,K=NULL,
            min.MAF=.Machine$double.eps,n.PC=npcs,n.core=.coresFrac(rc),
            P3D=TRUE,plot=FALSE)
    })
    
    rrb$pvalue <- 10^-rrb[,ncol(rrb)]
    disp("Found ",length(which(rrb$pvalue<psig))," associations ",
        "(uncorrected p-values).")
    
    return(rrb) # To be changed after statgenGWAS
}

gwaStatgen <- function(obj,response,covariates,pcs=TRUE,psig=0.05,rc=NULL) {
    if (!is(obj,"GWASExperiment"))
        stop("obj must be an object of class GWASExperiment!")
        
    .checkNumArgs("Association p-value",psig,"numeric",c(0,1),"both")
    
    # Preprocess phenotypes, similarly to gwaGlm
    p <- phenotypes(obj)
    chResCov <- .validateResponseAndCovariates(p,response,covariates)
    response <- chResCov$res
    covariates <- chResCov$cvs
    p <- p[,c(response,covariates)]
    
    # Convert to gData object
    disp("Preparing data for statgenGWAS...")
    sg <- GWASExperiment2gData(obj,covariates,pcs)
    
    disp("Performing GWAS with statgenGWAS model")
    disp("Trait(s)    : ",response)
    disp("Covariate(s): ",paste0(covariates,collapse=", "))
    disp("Use PCs     : ",ifelse(pcs,"Yes","No"))
    
    ssres <- runSingleTraitGwas(
        gData=sg,
        traits=response,
        covar=seq_len(ncol(sg$covar)),
        MAF=.Machine$double.eps,
        thrType="fixed",
        LODThr=-log10(psig),
        nCores=.coresFrac(rc)
    )

    ssfinal <- as.data.frame(ssres$GWAResult[[1]])
    rownames(ssfinal) <- ssfinal$snp
    ssfinal <- ssfinal[,c("effect","effectSe","pValue")]
    
    disp("Found ",length(which(ssfinal$pValue<psig))," associations ",
        "(uncorrected p-values).")
    
    return(ssfinal) # To be changed after statgenGWAS
}

.prepareGenotypesForBlup <- function(obj) {
    tmp <- assay(obj,1)
    geno <- as(tmp,"numeric")
    geno[geno==0] <- -1
    geno[geno==1] <- 0
    geno[geno==2] <- 1
    geno <- as.data.frame(geno)
    phold <- data.frame(snp=rownames(geno),chrom=rep(NA,nrow(geno)),
        bp=rep(NA,nrow(geno)))
    return(cbind(phold,geno))
}

.estimateNPCinBlup <- function(pheno,geno,fixed,n=20) {
    nsig <- numeric(n)
    for (i in 1:n) {
        disp("  Testing rrBLUP with ",i," PCs")
        log <- capture.output({
            tmp <- GWAS(pheno,geno,fixed=fixed,n.core=.coresFrac(rc),P3D=TRUE,
                plot=FALSE,min.MAF=.Machine$double.eps,n.PC=i)
        })
        nsig[i] <- length(which(tmp[,ncol(tmp)]> -log10(0.05)))
    }
    # Select the first closest to the median
    m <- abs(nsig - ceiling(median(nsig)))
    return(which(m==min(m))[1])
}

.validateResponseAndCovariates <- function(p,res,cvs) {
    if (res %in% cvs)
        stop("Response variable ",res," cannot be also a covariate!")
            
    if (is.numeric(res))
        res <- names(p)[res] # Will be NA if out of bounds
    if (!(res %in% names(p)))
        stop("Response variable ",res," could not be found in the input ",
            "GWASExperiment phenotypes!")
    
    if (is.numeric(cvs))
        cvs <- names(p)[cvs]
    if (!all(cvs %in% names(p))) {
        nf <- cvs[!(cvs %in% names(p))]
        stop("Covariates ",paste0(nf,collapse=", ")," could not be found in ",
            "the input GWASExperiment phenotypes!")
    }
    
    return(list(res=res,cvs=cvs))
}
