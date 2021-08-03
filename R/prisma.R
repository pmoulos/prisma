prsPipeline <- function(
    gwe,
    phenotype,
    covariates,
    pcs=FALSE,
    npcs=0,
    snpSelection=NULL,
    trainSize=0.8,
    niter=10,
    filters=getDefaults("filters"),
    pcaMethod=c("auto","snprel","grid","hubert"),
    imputeMissing=FALSE,
    imputeMethod=c("single","split"),
    gwaMethods=c("glm","rrblup","statgen","snptest"), # lasso later
    gwaCombine=c("fisher","simes","max","min","harmonic","whitlock","pandora"),
    family=NULL,
    glmOpts=getDefaults("glm"),
    rrblupOpts=getDefaults("rrblup"),
    statgenOpts=getDefaults("statgen"),
    snptestOpts=getDefaults("snptest"),
    prsMethods=c("lassosum","prsice"),
    lassosumOpts=getDefaults("lassosum"),
    prsiceOpts=getDefaults("prsice"),
    prsWorkspace=NULL,
    cleanup=c("none","intermediate","all"),
    logging=c("screen","sink"),
    rc=NULL
) {
    # Can we run a GWA?
    .canRunGwa(gwe)
    
    # Validate rest
    pcaMethod <- pcaMethod[1]
    imputeMethod <- imputeMethod[1]
    gwaCombine <- gwaCombine[1]
    cleanup <- cleanup[1]
    
    .checkNumArgs("Number of iterations (niter)",niter,"numeric",0,"gt")
    .checkNumArgs("Training dataset size fraction (trainSize)",trainSize,
        "numeric",0,"gt")
    .checkTextArgs("PCA method (pcaMethod)",pcaMethod,c("auto","snprel",
        "grid","hubert"),multiarg=FALSE)
    .checkNumArgs("Number of PCs (npcs)",npcs,"numeric",0,"gte")
    .checkTextArgs("Imputation mode (mode)",imputeMethod,c("single","split"),
        multiarg=FALSE)
    .checkTextArgs("GWA methods (gwaMethods)",gwaMethods,
        c("glm","rrblup","statgen","snptest"),multiarg=TRUE)
    .checkTextArgs("p-value combination (gwaCombine)",gwaCombine,
        c("fisher","simes","max","min","harmonic","whitlock","pandora"),
        multiarg=FALSE)
    .checkTextArgs("Logging option",logging,c("screen","sink"),multiarg=FALSE)
    .checkTextArgs("Cleanup option",cleanup,c("none","intermediate","all"),
        multiarg=FALSE)
    
    glmOpts <- .checkGwaArgs(glmOpts,"glm")
    rrblupOpts <- .checkGwaArgs(rrblupOpts,"rrblup")
    #statgenOpts <- .checkGwaArgs(glmOpts,"statgen")
    snptesOpts <- .checkGwaArgs(snptestOpts,"snptest")
    lassosumOpts <- .checkPrsArgs(lassosumOpts,"lassosum")
    prsiceOpts <- .checkPrsArgs(prsiceOpts,"prsice")
    
    # Workspace valid?
    prsWorkspace <- .validateWorkspacePath(prsWorkspace,"prisma")
    
    # Validate response and covariates
    p <- phenotypes(gwe)
    chResCov <- .validateResponseAndCovariates(p,phenotype,covariates)
    phenotype <- chResCov$res
    covariates <- chResCov$cvs
    
    # Missing value imputation before all else
    if (imputeMissing) {
        gwe <- imputeGWAS(gwe,mode=imputeMethod)
        imputeMissing <- FALSE # No need to pass downstream!
    }
    
    # Check if the filters are given properly
    filters <- .checkFilters(filters)
    
    # Main iteration
    if (is.null(snpSelection))
        theResult <- .prsPipelineDenovo(gwe,phenotype,covariates,pcs,npcs,
            trainSize,niter,filters,pcaMethod,imputeMissing,imputeMethod,
            gwaMethods,gwaCombine,family,glmOpts,rrblupOpts,
            statgenOpts,snptestOpts,prsMethods,lassosumOpts,prsiceOpts,
            prsWorkspace,logging,rc)
    else
        # This time is a PGS Catalog data frame
        theResult <- .prsPipelineExternal(gwe,phenotype,covariates,pcs,npcs,
            snpSelection,trainSize,niter,filters,pcaMethod,imputeMissing,
            imputeMethod,gwaMethods[1],family,gwaOpts,prsiceOpts,prsWorkspace,
            logging,rc)

    switch(cleanup,
        none = {
            disp("All pipeline output can be found at ",prsWorkspace)
        },
        intermediate = {
            disp("Cleaning up temporary program-specific intermediate files ",
                "from ",prsWorkspace)
            # Some routine to delete PLINK and summary stats etc.
        },
        all = {
            disp("Cleaning up all pipeline files!")
            unlink(prsWorkspace,recursive=TRUE,force=TRUE)
            theResult <- lapply(theResult,function(x) {
                attr(x,"workspace") <- NULL
            })
        }
    )
    
    return(theResult)
}

aggregatePrsMarkers <- function(gwaList,mode=c("intersect","union"),qcut=0.9) {
    # Check if prsbetas non-empty everywhere
    check <- vapply(gwaList,function(x) {
        return(is.null(prsbetas(x)))
    },logical(1))
    if (any(check)) {
        w <- seq_along(check)[check]
        warning(length(w)," objects do not have an associated PRS analysis! ",
            "These are ",paste(w,collapse=", ")," and will be discarded")
        gwaList <- gwaList[!check]
    }
    
    # Argument validation
    mode <- mode[1]
    .checkTextArgs("Aggregation mode (mode)",mode,c("intersect","union"),
        multiarg=FALSE)
    .checkNumArgs("Quantile cutoff",qcut,"numeric",c(0,1),"both")
    
    # TODO: Argument to control which association to choose from?
    preCandidates <- lapply(gwaList,function(x,m) {
        b <- prsbetas(x)
        g <- lapply(colnames(b),function(n,b) {
            s <- rownames(b)[which(b[,n] != 0)]
        },b)
        return(Reduce(mode,g))
    },mode)
    prsCandidates <- unlist(preCandidates)
    freq <- table(prsCandidates)
    goods <- names(freq)[freq >= floor(quantile(freq,qcut))]
    
    out <- as.data.frame(freq[goods])
    names(out) <- c("snp","freq")
    rownames(out) <- out$snp
    return(out[order(out$freq,decreasing=TRUE),,drop=FALSE])
}

.prsPipelineDenovo <- function(gwe,phenotype,covariates,pcs,npcs,trainSize,
    niter,filters,pcaMethod,imputeMissing,imputeMethod,gwaMethods,gwaCombine,
    family,glmOpts,rrblupOpts,statgenOpts,snptestOpts,prsMethods,lassosumOpts,
    prsiceOpts,prsWorkspace,logging,rc) {
    # Initialize the list of GWASExperiment s
    theResult <- vector("list",niter)
    
    # For output directory (workspace) format
    dig <- nchar(as.character(niter))
    pcno <- 0
    
    # The worker
    for (i in seq_len(niter)) {
        pad <- paste0(rep("0",dig - nchar(as.character(i))),collapse="")
        iterWspace <- file.path(prsWorkspace,paste0(format(Sys.time(),
            "%Y%m%d%H%M%S"),"_denovo_",pad,i))
        dir.create(iterWspace,recursive=TRUE,showWarnings=FALSE)
        
        if (logging == "sink") {
            fh <- file(file.path(prsWorkspace,paste0("iter_",pad,i,".log")),
                open="wt")
            sink(fh,type="message")
        }
        
        disp("\n==============================================================")
        disp("-----> Denovo pipeline iteration ",i)
        disp("==============================================================\n")
        
        # Partition the object
        disp("\n----- Dataset partitioning -----\n")
        tmp <- partitionGWAS(gwe,by=phenotype,n=1,frac=trainSize,out="ttboth")
        theTrain <- tmp$train
        theTest <- tmp$test
        
        # Base QC
        disp("\n----- Training (base) QC -----\n")
        theTrain <- filterGWAS(theTrain,filters=filters,
            imputeMissing=imputeMissing)
        if (pcs) {
            disp("\n----- Training (base) PCA -----\n")
            theTrain <- calcPcaCovar(theTrain,method=pcaMethod,npc=npcs)
            pcno <- ncol(pcaCovariates(theTrain))
        }
                
        # Target QC
        disp("\n----- Test (target) QC -----\n")
        theTest <- filterGWAS(theTest,filters=filters,
            imputeMissing=imputeMissing)
        if (pcs) {
            disp("\n----- Test (target) PCA -----\n")
            theTest <- calcPcaCovar(theTest,method=pcaMethod,npc=pcno)
        }
        
        # Run GWAS on base
        disp("\n----- Training (base) GWAS -----\n")
        # TODO: Pass method options - defaults for time being
        theTrain <- gwa(theTrain,phenotype,covariates,pcs=pcs,
            methods=gwaMethods,combine=gwaCombine,usepcblup="rrint",
            npcsblup=pcno,rc=rc)
        
        # Run PRS
        disp("\n----- PRS analysis with base and target -----\n")
        theTest <- runPRS(base=theTrain,target=theTest,phenotype,covariates,
            pcs=TRUE,methods=prsMethods,prsiceOpts=prsiceOpts,wspace=iterWspace,
            rc=rc)
        
        disp("==============================================================\n")
        
        if (logging == "sink")
            sink(type="message")
        
        # A temporary attribute to the GWASExperiment output to maintain the
        # workspace directory for evaluation. Will be removed if super clean
        attr(theTest,"workspace") <- iterWspace
        # And also keep the PRSice R2 triplet
        if ("prsice" %in% prsMethods)
            attr(theTest,"PR2") <- .readPrsiceR2(iterWspace)
        
        theResult[[i]] <- theTest
    }
    
    return(theResult)
}

.prsPipelineExternal <- function(gwe,phenotype,covariates,pcs,npcs,snpSelection,
    trainSize,niter,filters,pcaMethod,imputeMissing,imputeMethod,gwaMethod,
    family,gwaOpts,prsiceOpts,prsWorkspace,logging,rc) {
    # Initialize the list of GWASExperiment s
    theResult <- vector("list",niter)
    
    # For output directory (workspace) format
    dig <- nchar(as.character(niter))
    pcno <- 0
    
    # The worker
    for (i in seq_len(niter)) {
        pad <- paste0(rep("0",dig - nchar(as.character(i))),collapse="")
        iterWspace <- file.path(prsWorkspace,paste0(format(Sys.time(),
            "%Y%m%d%H%M%S"),"_external_",pad,i))
        dir.create(iterWspace,recursive=TRUE,showWarnings=FALSE)
        
        if (logging == "sink") {
            fh <- file(file.path(prsWorkspace,paste0("iter_",pad,i,".log")),
                open="wt")
            sink(fh,type="message")
        }
        
        disp("\n==============================================================")
        disp("-----> External pipeline iteration ",i)
        disp("==============================================================\n")
        
        pad <- paste0(rep("0",dig - nchar(as.character(i))),collapse="")
        iterWspace <- file.path(prsWorkspace,paste0(format(Sys.time(),
            "%Y%m%d%H%M%S"),"_validation_",pad,i))
        dir.create(iterWspace,recursive=TRUE)
        
        # Final validation with PRSice
        disp("\n----- Dataset partitioning -----\n")
        tmp <- partitionGWAS(gwe,by=phenotype,n=1,frac=trainSize,out="ttboth")
        base <- tmp$train
        target <- tmp$test

        # Base QC
        disp("\n----- Base QC -----\n")
        base <- filterGWAS(base,filters=filters,imputeMissing=imputeMissing)
        if (pcs) {
            disp("\n----- Base PCA -----\n")
            base <- calcPcaCovar(base,method=pcaMethod,npc=npcs)
            pcno <- ncol(pcaCovariates(base))
        }

        # Target QC
        disp("\n----- Target QC -----\n")
        target <- filterGWAS(target,filters=filters,imputeMissing=imputeMissing)
        if (pcs) {
            disp("\n----- Target PCA -----\n")
            target <- calcPcaCovar(target,method=pcaMethod,npc=pcno)
        }

        # Run basic GWAS on base to get coefficients
        # TODO: Pass method options - defaults and glm for time being
        disp("\n----- Base GWAS -----\n")
        base <- gwa(base,phenotype,covariates,pcs=pcs,methods="glm",rc=rc)

        # Run PRS with a safeguard for possibly filtered SNPs
        safeguard <- Reduce("intersect",list(snpSelection,rownames(base),
            rownames(target)))
        if (length(safeguard) < length(snpSelection)) {
            disp("  ",length(snpSelection) - length(safeguard)," SNPs not ",
                "found in base data! Excluding...")
            subBase <- base[safeguard,,drop=FALSE]
            subTarget <- target[safeguard,,drop=FALSE]
            if (pcs) {
                disp("  ...and recalculating PCs...")
                subBase <- calcPcaCovar(subBase,method=pcaMethod,npc=pcno)
                subTarget <- calcPcaCovar(subBase,method=pcaMethod,npc=pcno)
            }
        }
        else {
            subBase <- base
            subTarget <- target
        }
        
        disp("\n----- PRS analysis with base and target -----\n")
        prsiceOut <- prsicePRS(base=subBase,target=subTarget,response=phenotype,
            covariates=covariates,pcs=pcs,mode="apply",wspace=iterWspace,rc=rc)
        
        # A temporary attribute to the GWASExperiment output to maintain the
        # workspace directory for evaluation. Will be removed if super clean
        attr(prsiceOut,"workspace") <- iterWspace
        # And also keep the PRSice R2 triplet
        attr(prsiceOut,"PR2") <- .readPrsiceR2(iterWspace)
        
        theResult[[i]] <- prsiceOut
    }
    
    return(theResult)
}

.readPrsiceR2 <- function(wspace) {
    sum <- dir(wspace,pattern=".summary$",full.names=TRUE)
    if (is.character(sum) && file.exists(sum)) {
        tmp <- read.delim(sum)
        return(c(r2m=tmp$Full.R2,r2n=tmp$Null.R2,r2p=tmp$PRS.R2))
    }
    else
        return(c(r2m=NA,r2n=NA,r2p=NA))
}

.getR2 <- function(obj) {
    tmp <- lapply(obj,function(x) {
        if (!is.null(attr(x,"PR2")))
            return(attr(x,"PR2"))
        else {
            if (!is.null(attr(x,"workspace")))
                return(.readPrsiceR2(attr(x,"workspace")))
            else
                return(c(r2m=NA,r2n=NA,r2p=NA))
        }
    })
    return(do.call("rbind",tmp))
}
