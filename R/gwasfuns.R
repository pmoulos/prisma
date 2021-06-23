gwa <- function(obj,response,covariates=NULL,pcs=FALSE,psig=0.05,
    methods=c("glm","rrblup","statgen","snptest"),
    combine=c("fisher","simes","max","min","harmonic","whitlock","pandora"),
    #args=getDefaults("gwargs"),
    family=NULL,usepcblup=c("auto","estim","fixed","none"),npcsblup=NULL,
    rc=NULL,...) {
    if (!is(obj,"GWASExperiment"))
        stop("obj must be an object of class GWASExperiment!")
    
    combine <- combine[1]
    usepcblup <- usepcblup[1]
    .checkNumArgs("Association p-value",psig,"numeric",c(0,1),"both")
    if (!is.null(npcsblup))
        .checkNumArgs("rrBLUP number of PCs",npcsblup,"numeric",0,"gt")
    .checkTextArgs("GWAS methods",methods,c("glm","rrblup","statgen","snptest"),
        multiarg=TRUE)
    .checkTextArgs("p-value combination",combine,c("fisher","simes","max","min",
        "harmonic","whitlock","pandora"),multiarg=FALSE)
    .checkTextArgs("PCA in rrBLUP",usepcblup,c("auto","estim","fixed","none"),
        multiarg=FALSE)
    
    # Initialize pvalue matrix
    pMatrix <- matrix(1,nrow=nrow(obj),ncol=length(methods))
    rownames(pMatrix) <- rownames(theTrain)
    colnames(pMatrix) <- methods
    
    # Run GWAS
    for (m in methods) {
        switch(m,
            glm = {
                glmRes <- gwaGlm(obj,response,covariates,pcs,family,psig,rc)
                pMatrix[,m] <- glmRes[,4]
            },
            rrblup = {
                rrbRes <- gwaBlup(obj,response,covariates,usepc=usepcblup,
                    npcs=npcsblup,psig=psig,rc=rc)
                pMatrix[,m] <- rrbRes$pvalue
            },
            statgen = {
                sgRes <- gwaStatgen(obj,response,covariates,pcs,psig,rc)
                pMatrix[,m] <- sgRes$pValue
            },
            snptest = {
                sgSnp <- gwaSnptest(obj,response,covariates,pcs,psig)
                pMatrix[,m] <- sgSnp$pvalue
            }
        )
    }
    
    # Combine p-values
    if (length(methods) > 1) {
        disp("Combining p-values from multi-GWAS using ",combine," method")
        switch(combine,
            fisher = {
                tmp <- fisherMethod(pMatrix,p.corr="none",
                    zeroSub=.Machine$double.xmin)
                pComb <- tmp$p.value
            },
            simes = {
                pComb <- apply(pMatrix,1,combineSimes)
            },
            max = {
                pComb <- apply(pMatrix,1,combineMaxp)
            },
            min = {
                pComb <- apply(pMatrix,1,combineMinp)
            },
            harmonic = {
                # NYI until weights
            },
            whitlock = {
                # NYI until weights
            },
            pandora = {
                # NYI until weights
            }
        )
    }
    else
        pComb <- pMatrix
   
    disp("Overall, found ",length(which(pComb<psig))," associations with ",
        "(uncorrected) combined p-values).")
   return(pComb)
}

gwaGlm <- function(obj,response,covariates=NULL,pcs=FALSE,family=NULL,psig=0.05,
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
    # We must ensure that if binomial requested, res is 0-1 and binary
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
            warning("PC covariates reqeusted in the model, but not calculated ",
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

gwaBlup <- function(obj,response,covariates=NULL,usepc=c("auto","estim","fixed",
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
        npcs <- .estimateNPCinBlup(pheno,geno,covariates,rc=rc)
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

gwaStatgen <- function(obj,response,covariates=NULL,pcs=TRUE,psig=0.05,
    rc=NULL) {
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

gwaSnptest <- function(obj,response,covariates=NULL,pcs=TRUE,psig=0.05,
    test=c("frequentist","bayesian"),model=c("additive","dominant","recessive",
    "general","heterozygote")) {
    # Tool availability before all else
    if (!.toolAvailable("snptest"))
        stop("SNPTEST program not found in the system!")
    
    test <- test[1]
    .checkTextArgs("SNPTEST testing",test,c("frequentist","bayesian"),
        multiarg=FALSE)
    model <- model[1]
    .checkTextArgs("SNPTEST type",model,c("additive","dominant","recessive",
        "general","heterozygote"),multiarg=FALSE)
    
    modelCode <- c(additive=1,dominant=2,recessive=3,general=4,heterozygote=5)
    
    # 1. Prepare SNPTEST files (from PLINK)
    # 1a. Create a SNPTEST workspace in Rtmp
    stwork <- file.path(tempdir(),paste0("snptest_workspace_",.randomString()))
    if (!dir.exists(stwork))
        dir.create(stwork,recursive=TRUE)
    # 1b. obj to PLINK through snpStats::write.plink
    # 1c. Necessary phenotypes and covariates to SNPTEST sample file
    tmplink <- .preparePlinkInputForSnptest(obj,response,covariates,pcs,stwork)
    
    # 2. Run the SNPTEST command
    disp("Performing GWAS with SNPTEST")
    disp("Trait(s)    : ",response)
    disp("Covariate(s): ",paste0(covariates,collapse=", "))
    disp("Use PCs     : ",ifelse(pcs,"No","Yes"))
    
    snptest <- .findSnptest()
    command <- paste0(snptest," -data ",tmplink$plink," ",tmplink$sample,
        " -frequentist ",modelCode[model]," -method score -pheno ",response,
        " -cov_all -o ",file.path(stwork,"snptest.out")
    )
    
    message("Executing: ",command)
    out <- tryCatch(system(command),error=function(e) {
        message("Caught error: ",e$message)
        return(1L)
    },finally="")
    
    if (out == 1L) {
        message("SNPTEST run failed! Will return unit p-values...")
        return(rep(1,nrow(obj)))
    }
    
    # 3. Read and gather SNPTEST output in a table similar to others
    tmp <- read.table(file.path(stwork,"snptest.out"),header=TRUE,sep=" ")
    res <- tmp[,c(23,24,21)]
    colnames(res) <- c("effect","se","pvalue")
    rownames(res) <- rownames(obj)
    
    disp("Found ",length(which(res[,3]<psig))," associations (uncorrected ",
        "p-values).")
    return(res)
}

.preparePlinkInputForSnptest <- function(obj,response,covariates,pcs,wspace) {
    # Retrieve PLINK related data
    geno <- t(assay(obj,1))
    fam <- gsamples(obj)
    map <- gfeatures(obj)
    
    # Subset phenotypes and add PCs if requested
    p <- phenotypes(obj)
    chResCov <- .validateResponseAndCovariates(p,response,covariates)
    response <- chResCov$res
    covariates <- chResCov$cvs
    p <- p[,c(response,covariates)]
    ct <- .initSnptestSampleFirstRow(p,response,covariates)
    if (pcs) { # Include robust PCs in the model
        if (.hasPcaCovariates(obj)) {
            pcov <- pcaCovariates(obj)
            p <- cbind(p,pcov)
            ct <- c(ct,rep("C",ncol(pcov)))
        }
        else
            warning("PC covariates reqeusted in the model, but not calculated ",
                "PC covariates found! Ignoring...",immediate.=TRUE)
    }
    
    # The plink files
    disp("Writing PLINK files in ",wspace)
    write.plink(
        file.base=file.path(wspace,"forsnptest"),
        snps=geno,
        pedigree=fam$pedigree,
        id=rownames(fam),
        father=fam$father,
        mother=fam$mother,
        sex=fam$sex,
        phenotype=fam$affected,
        chromosome=map$chromosome,
        position=map$position,
        allele.1=map$allele.1,
        allele.2=map$allele.2
    )
    
    # Modify the phenotypes data frame
    disp("Writing SNPTEST sample file in ",wspace)
    p <- data.frame(sample_id=rownames(p),p)
    ct <- c("0",ct)
    sfile <- file.path(wspace,"pheno.sample")
    r1 <- paste0(colnames(p),collapse=" ")
    r2 <- paste0(ct,collapse=" ")
    writeLines(c(r1,r2),sfile)
    write.table(p,file=sfile,append=TRUE,quote=FALSE,sep=" ",row.names=FALSE,
        col.names=FALSE)
    
    return(list(
        plink=file.path(wspace,"forsnptest.bed"),
        sample=file.path(wspace,"pheno.sample")
    ))
}

.initSnptestSampleFirstRow <- function(p,response,covariates) {
    # Initialize the column type vector
    r <- p[,response,drop=FALSE]
    ctr <- character(ncol(r))
    for (i in seq_len(ncol(r))) {
        tmp <- r[,i]
        if (is.numeric(tmp))
            ctr[i] <- "P"
        else if (is.integer(tmp) || is.character(tmp) || is.factor(tmp)) {
            if (length(unique(tmp)) == 2)
                ctr[i] <- "B"
            else
                ctr[i] <- "D"
        }
    }
    r <- p[,covariates,drop=FALSE]
    ctc <- character(ncol(r))
    for (i in seq_len(ncol(r))) {
        tmp <- r[,i]
        if (is.numeric(tmp)) {
            if (is.integer(tmp)) {
                if (length(unique(tmp)) <= 10)
                    ctc[i] <- "D"
                else
                    ctc[i] <- "C"
            }
            else
                ctc[i] <- "C"
        }
        else if (is.character(tmp) || is.factor(tmp))
            ctc[i] <- "D"
    }
    return(c(ctr,ctc))
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

.estimateNPCinBlup <- function(pheno,geno,fixed,n=20,rc=NULL) {
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

.validateResponseAndCovariates <- function(p,res,cvs=NULL) {
    if (res %in% cvs)
        stop("Response variable ",res," cannot be also a covariate!")
            
    if (is.numeric(res))
        res <- names(p)[res] # Will be NA if out of bounds
    if (!(res %in% names(p)))
        stop("Response variable ",res," could not be found in the input ",
            "GWASExperiment phenotypes!")
    
    if (is.numeric(cvs))
        cvs <- names(p)[cvs]
    if (is.character(cvs) && !all(cvs %in% names(p))) {
        nf <- cvs[!(cvs %in% names(p))]
        stop("Covariates ",paste0(nf,collapse=", ")," could not be found in ",
            "the input GWASExperiment phenotypes!")
    }
    
    return(list(res=res,cvs=cvs))
}
