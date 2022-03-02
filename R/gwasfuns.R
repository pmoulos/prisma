gwa <- function(obj,response,covariates=NULL,pcs=FALSE,psig=0.05,
    methods=c("glm","rrblup","statgen","snptest","plink","lasso"),
    combine=c("fisher","simes","max","min","harmonic","whitlock","pandora"),
    glmOpts=getDefaults("glm"),rrblupOpts=getDefaults("rrblup"),
    statgenOpts=getDefaults("statgen"),snptestOpts=getDefaults("snptest"),
    plinkOpts=getDefaults("plink"),rc=NULL,...) {
    # Cannot run GWA in an empty object without genotypes/phenotypes or if the
    # object is not of class GWASExperiment
    .canRunGwa(obj)
    
    # We might need to deactivate parallel calculations for GLM if called from
    # the total PRISMA pipeline
    noGlmPrl <- grepl("doTryCatch|eval\\(parse",
        deparse(sys.calls()[sys.nframe()-1])[1])
    
    # Check input variables
    combine <- combine[1]
    .checkNumArgs("Association p-value",psig,"numeric",c(0,1),"both")
    .checkTextArgs("GWAS methods (methods)",methods,c("glm","rrblup","statgen",
        "snptest","plink","lasso"),multiarg=TRUE)
    .checkTextArgs("p-value combination",combine,c("fisher","simes","max","min",
        "harmonic","whitlock","pandora"),multiarg=FALSE)
    
    glmOpts <- .checkGwaArgs(glmOpts,"glm")
    rrblupOpts <- .checkGwaArgs(rrblupOpts,"rrblup")
    #statgenOpts <- .checkGwaArgs(statgenOpts,"statgen")
    snptesOpts <- .checkGwaArgs(snptestOpts,"snptest")
    plinkOpts <- .checkGwaArgs(plinkOpts,"plink")
    
    # rrBLUP does not like binary phenotypes
    if ("rrblup" %in% methods) {
        if (rrblupOpts$pcblup == "rrint" && is.null(rrblupOpts$npcs) && pcs) {
            if (.hasPcaCovariates(obj))
                rrblupOpts$npcs <- ncol(pcaCovariates(obj))
            else
                rrblupOpts$npcs <- 0
        }
        
        if (.maybeBinaryForBinomial(phenotypes(obj)[,response])) {
            warning("The trait ",response," is binary and not supported by ",
                "rrBLUP which will be removed.\nConsider vanilla ridge ",
                "regression (glmnet) if not already included.",immediate.=TRUE)
            methods <- methods[methods!="rrblup"]
        }
    }
    
    # Splits the association analysis per chromosome, ortherwise, some chunks
    map <- gfeatures(obj)
    if ("chromosome" %in% names(map))
        parts <- split(obj,map$chromosome)
    else {
        splitFactor <- .splitFactorForParallel(nrow(obj),rc)
        parts <- split(obj,splitFactor)
    }
    
    # Parallelization inside chunks where possible
    disp("\nStarting GWA analysis in ",length(parts)," chunks with ",
        length(methods)," methods: ",paste0(methods,collapse=", "))
    P <- lapply(names(parts),function(x,prt) {
        disp("\n========== Analyzing chromosome/part ",x)
        o <- prt[[x]]
        y <- .gwaWorker(o,response,covariates,pcs,psig,methods,combine,
            glmOpts$family,rrblupOpts$pcblup,rrblupOpts$npcs,
            snptestOpts$test,snptestOpts$model,snptestOpts$workspace,
            plinkOpts$effect,plinkOpts$seed,plinkOpts$workspace,
            plinkOpts$cliopts,glmOpts$size,noGlmPrl,rc)
        gc(verbose=FALSE)
        return(y)
    },parts)
    names(P) <- names(parts)
    
    disp("\nFinished, the following putative associations were found per part:")
    lapply(names(P),function(x,s,pop) {
        disp("  Part ",x,": ",length(which(pop[[x]][,ncol(pop[[x]])] < s)),
            " associations found out of ",nrow(pop[[x]])," markers")
    },psig,P)
    
    tmp <- do.call("rbind",P)
    eassoc <- tmp[,seq_along(methods),drop=FALSE]
    passoc <- tmp[,(length(methods)+1):ncol(tmp),drop=FALSE]
    disp("\nOverall, found ",length(which(passoc[,ncol(passoc)]<psig)),
        " associations with (uncorrected) combined p-values).")
    
    # Covariates will be sorted prior to storing, if previously exist, replace
    npcs <- ifelse(pcs && .hasPcaCovariates(obj),ncol(pcaCovariates(obj)),0)
    pvalues(obj,response,covariates,npcs) <- passoc
    effects(obj,response,covariates,npcs) <- eassoc
    
    return(obj)
}

# Returns only the p-value matrix. For assigning p-values to object, use GWA
.gwaWorker <- function(obj,response,covariates,pcs,psig,methods,combine,
    family,usepcblup,npcsblup,snptest,snpmodel,snpspace,pleff,plseed,plspace,
    plcliopts,size,noglmprl,rc=NULL,...) {
    # Initialize pvalue matrix
    pMatrix <- eMatrix <- matrix(1,nrow=nrow(obj),ncol=length(methods))
    rownames(pMatrix) <- rownames(eMatrix) <- rownames(obj)
    colnames(pMatrix) <- colnames(eMatrix) <- methods
    
    # Run GWAS
    rcLocal <- rc
    if (noglmprl)
        rcLocal <- NULL
        
    for (m in methods) {
        switch(m,
            glm = {
                glmRes <- gwaGlm(obj,response,covariates,pcs,family,psig,
                    penalized=FALSE,size,rc=rcLocal)
                pMatrix[,m] <- glmRes[,4]
                eMatrix[,m] <- glmRes[,1]
            },
            rrblup = {
                rrbRes <- gwaBlup(obj,response,covariates,usepc=usepcblup,
                    npcs=npcsblup,psig=psig,rc=rcLocal)
                pMatrix[,m] <- rrbRes$pvalue
                eMatrix[,m] <- rrbRes$effect
            },
            statgen = {
                sgRes <- gwaStatgen(obj,response,covariates,pcs,psig,rcLocal)
                pMatrix[,m] <- sgRes$pValue
                eMatrix[,m] <- sgRes$effect
            },
            snptest = {
                stRes <- gwaSnptest(obj,response,covariates,pcs,psig,snptest,
                    snpmodel,snpspace)
                pMatrix[,m] <- stRes$pvalue
                eMatrix[,m] <- stRes$effect
            },
            plink = {
                plRes <- gwaPlink(obj,response,covariates,pcs,psig,plseed,
                    plcliopts,plspace,rcLocal)
                pMatrix[,m] <- plRes$pvalue
                eMatrix[,m] <- plRes$effect
            },
            lasso = {
                lsRes <- gwaGlm(obj,response,covariates,pcs,family,psig,size,
                    penalized=TRUE,rc=rcLocal)
                pMatrix[,m] <- lsRes[,4]
                eMatrix[,m] <- lsRes[,1]
            }
        )
    }
    
    # NA fix, pvalue 1, effect 0
    if (any(is.na(pMatrix)))
        pMatrix[is.na(pMatrix)] <- 1
    if (any(is.na(eMatrix)))
        eMatrix[is.na(eMatrix)] <- 0
    
    # Combine p-values
    if (length(methods) > 1) {
        disp("\nCombining p-values from multi-GWAS using ",combine," method")
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
   
    disp("\nOverall, found ",length(which(pComb<psig))," associations with ",
        "(uncorrected) combined p-values).")
   
    # Assign to object pvalues slot
    if (length(methods) > 1) {
        pMatrix <- cbind(pMatrix,pComb)
        colnames(pMatrix)[ncol(pMatrix)] <- combine
    }
    
    # Since this function is not used directly, easier to cbind effects and
    # p-values for later result display
    return(cbind(eMatrix,pMatrix))
}

gwaGlm <- function(obj,response,covariates=NULL,pcs=FALSE,family=NULL,psig=0.05,
    penalized=FALSE,size=1,rc=NULL,...) {
    .canRunGwa(obj)
    .checkNumArgs("Association p-value (psig)",psig,"numeric",c(0,1),"both")
    .checkNumArgs("SNPs to bucket test (size)",size, "numeric",0,"gt")
    
    # Construct the glm data frame, check if response and covariate names are
    # valid
    p <- phenotypes(obj)
    chResCov <- .validateResponseAndCovariates(p,response,covariates)
    response <- chResCov$res
    covariates <- chResCov$cvs
    
    # Reduce the phenotypes for testing and check GLM call integrity
    p <- p[,c(response,covariates),drop=FALSE]
    # We must ensure that if binomial requested, res is 0-1 and binary
    if (!is.null(family)) {
        family <- family[1]
        .checkTextArgs("Regression family",family,
            c("gaussian","binomial","poisson"),multiarg=FALSE)
        fpredHelper <- "provided"
    }
    else {
        family <- ifelse(.maybeBinaryForBinomial(p[,response]),"binomial",
            "gaussian")
        fpredHelper <- "predicted"
    }
    
    # If binomial decided, check response and convert if necessary
    if (family == "binomial")
        p[,response] <- .validateBinaryForBinomial(p[,response])
    
    if (pcs) { # Include robust PCs in the model
        if (.hasPcaCovariates(obj)) {
            pcov <- pcaCovariates(obj)
            p <- cbind(p,pcov)
        }
        else
            warning("PC covariates requested in the model, but no calculated ",
                "PC covariates found! Ignoring...",immediate.=TRUE)
    }
    
    # We now need to split the SnpMatrix, as glm is quite fast, parallelizng
    # per snp will most probably cause more overhead than chunking and iterating
    # We split according to the number of available cores for analysis.
    snps <- genotypes(obj)
    ##!!!
    #snps <- switch.alleles(snps,seq_len(ncol(snps)))
    ##!!!
    splitFactor <- .splitFactorForParallel(nrow(snps),rc)
    # We now need to split the SnpMatrix, as glm is quite fast, parallelizng
    # per snp will most probably cause more overhead than chunking and iterating
    # We split according to the number of available cores for analysis.
    #splits <- split(seq_len(nrow(snps)),splitFactor)
    splits <- split(rownames(snps),splitFactor)
    
    wh <- ifelse(penalized,"lasso","glm")
    whEx <- ifelse(penalized,"Penalized Regression (LASSO)","GLM")
    disp("\nPerforming GWAS with ",whEx," over ",length(splits)," chunks")
    disp("Trait(s)    : ",response)
    disp("Covariate(s): ",paste0(covariates,collapse=", "))
    disp("Regression  : ",family," (",fpredHelper,")")
    disp("Use PCs     : ",ifelse(pcs,"Yes","No"))
        
    if (size == 1)
        assocList <- cmclapply(splits,function(n,p,r,f,w,snps,...) {
            disp("  testing SNPs from ",n[1]," to ",n[length(n)])
            batch <- lapply(n,.glmTestSingleSNP,p,r,f,w,snps,...)
            batch <- do.call("rbind",batch)
            rownames(batch) <- n
            return(batch)
        },p,response,family,wh,snps,...,rc=rc)
    else
        assocList <- cmclapply(splits,function(n,p,r,f,w,m,snps,...) {
            disp("  testing SNPs from ",n[1]," to ",n[length(n)])
            # Further split per chunks of m - we should enforce ordering
            sf <- .splitFactorForParallel(length(n),floor(length(n)/m))
            inSplits <- split(n,sf)
            batch <- lapply(inSplits,.glmTestMultiSNP,p,r,f,w,snps,...)
            batch <- do.call("rbind",batch)
            rownames(batch) <- n
            return(batch)
        },p,response,family,wh,size,snps,...,rc=rc)
     
    assocResult <- do.call("rbind",assocList)
    scol <- ifelse(penalized,5,4)
    disp("Found ",length(which(assocResult[,scol]<psig))," associations ",
        "(uncorrected p-values).")
    
    return(assocResult)
}

.glmTestSingleSNP <- function(s,p,r,f,w,snps,...) {
    disp("    testing ",s,level="full")
    dat <- cbind(p,t(as(snps[s,],"numeric")))
    tmp <- coefficients(.gwaGlmWorker(w,dat,r,f,...))
    rownames(tmp)[nrow(tmp)] <- s
    return(tmp[s,,drop=FALSE])
}

.glmTestMultiSNP <- function(s,p,r,f,w,snps,...) {
    disp("    ",paste(s,collapse="\n    "),level="full")
    dat <- cbind(p,t(as(snps[s,],"numeric")))
    tmp <- coefficients(.gwaGlmWorker(w,dat,r,f,...))
    if (w=="glm") { # Deal with non-pruning collinearity
        ret <- matrix(NA,length(s),4)
        colnames(ret) <- colnames(tmp)
        rownames(ret) <- s
        com <- intersect(rownames(tmp),s)
        ret[com,] <- tmp[com,]
        return(ret)
    }
    else {
        rownames(tmp)[(nrow(tmp)-length(s)+1):nrow(tmp)] <- s
        return(tmp[s,,drop=FALSE])
    }
}

# ... is other options passed to GLM
.gwaGlmWorker <- function(wh,dat,res,fam,...) {
    # We assume that dat contains only the deserved covariates
    ii <- which(colnames(dat)==res)
    # Before constructing formula, all the covariates must have proper names
    colnames(dat) <- make.names(colnames(dat))
    covs <- colnames(dat)[-ii]
    cres <- colnames(dat)[ii]
    f <- as.formula(paste(cres,paste0(covs,collapse="+"),sep="~"))
    disp("Formula is: ",level="full")
    disp(show(f),level="full")
    if (wh=="glm")
        return(summary(glm(f,data=dat,family=fam,...)))
    else if (wh=="lasso") {
        ret <- tryCatch({
            summary(islasso(f,data=dat,family=fam,...))
        },error=function(e) {
            disp("Caught error: ",e$message)
            disp("Removing penalizing factor (lambda=0)")
            return(summary(islasso(f,data=dat,family=fam,unpenalized=covs,
                lambda=1,...)))
        },finally="")
        return(ret)
    }
}

gwaBlup <- function(obj,response,covariates=NULL,usepc=c("auto","estim",
    "rrint","fixed","none"),npcs=NULL,psig=0.05,rc=NULL) {
    .canRunGwa(obj)
    .checkNumArgs("Association p-value (psig)",psig,"numeric",c(0,1),"both")
    usepc <- usepc[1]    
    .checkTextArgs("Use PCA (usepc)",usepc,c("auto","estim","rrint","fixed",
        "none"),multiarg=FALSE)
    if (!is.null(npcs)) {
        npcs <- npcs[1]
        .checkNumArgs("Number of PCs (npcs)",npcs,"numeric",0,"gte")
    }
    else {
        if (usepc %in% c("fixed","rrint")) {
            warning("The number of PCs to use must be provided when ",
                "usepc is 'fixed' or 'rrint'! Switching to usepc = 'estim'...",
                immediate.=TRUE)
            usepc <- "estim"
        }
    }
    
    # Preprocess phenotypes, similarly to gwaGlm
    p <- phenotypes(obj)
    chResCov <- .validateResponseAndCovariates(p,response,covariates)
    response <- chResCov$res
    covariates <- chResCov$cvs
    p <- p[,c(response,covariates),drop=FALSE]
    
    later <- FALSE
    switch(usepc,
        auto = {
            # If pcaCovariates exist, use as many PCs as in there (e.g. from 
            # TW) otherwise, let rrBLUP estimate
            if (.hasPcaCovariates(obj)) {
                pcov <- pcaCovariates(obj)
                #npcs <- ncol(pcov)
                p <- cbind(p,pcov)
                covariates <- c(covariates,colnames(pcov))
                npcs <- 0
            }
            else
                usepc <- "estim"
        },
        estim = {
            later <- TRUE
        },
        rrint = {
            # Nothing, uses npcs without auto estimation
        },
        fixed = {
            if (.hasPcaCovariates(obj)) {
                pcov <- pcaCovariates(obj)
                p <- cbind(p,pcov)
                covariates <- c(covariates,colnames(pcov))
                npcs <- 0
            }
            else
                stop("When usepc = 'fixed', PCA must be performed prior to ",
                    "rrBLUP GWAS, e.g. with calcPcaCovar().")
        },
        none = {
            npcs <- 0
        }
    )
    
    # Add rownames at first column so as to be compatible with rrBLUP::GWAS
    pheno <- cbind(rownames(p),p)
    names(pheno)[1] <- "gid"
    
    # Prepare the genotypes
    disp("\nPreparing genotypes for rrBLUP...")
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
    disp("Use PCs     : ",ifelse(npcs==0 || is.null(npcs),"No","Yes"))
    log <- capture.output({
        rrb <- GWAS(pheno,geno,fixed=covariates,K=NULL,
            min.MAF=.Machine$double.eps,n.PC=npcs,n.core=.coresFrac(rc),
            P3D=TRUE,plot=FALSE)
        est <- mixed.solve(pheno[,response],Z=t(geno[,4:ncol(geno)]))
    })
    
    rrb$effect <- est$u
    rrb$pvalue <- 10^-rrb[,ncol(rrb)-1]
    disp("Found ",length(which(rrb$pvalue<psig))," associations ",
        "(uncorrected p-values).")
    
    return(rrb) # To be changed after statgenGWAS
}

gwaStatgen <- function(obj,response,covariates=NULL,pcs=TRUE,psig=0.05,
    rc=NULL) {
    .canRunGwa(obj)
    .checkNumArgs("Association p-value (psig)",psig,"numeric",c(0,1),"both")
    
    # Preprocess phenotypes, similarly to gwaGlm
    p <- phenotypes(obj)
    chResCov <- .validateResponseAndCovariates(p,response,covariates)
    response <- chResCov$res
    covariates <- chResCov$cvs
    p <- p[,c(response,covariates),drop=FALSE]
    # statgenGWAS works only with numeric phenotypes...
    p <- .checkAndCorrectFactorsFormat(p)
    phenotypes(obj) <- p
    
    # Convert to gData object
    disp("\nPreparing data for statgenGWAS...")
    #sg <- GWASExperiment2gData(obj,covariates,pcs,reverse=TRUE)
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
    "general","heterozygote"),wspace=NULL) {
    # Tool availability before all else
    if (!.toolAvailable("snptest"))
        stop("SNPTEST program not found in the system!")
    .canRunGwa(obj)
    
    test <- test[1]
    .checkTextArgs("SNPTEST testing (test)",test,c("frequentist","bayesian"),
        multiarg=FALSE)
    model <- model[1]
    .checkTextArgs("SNPTEST type (model)",model,c("additive","dominant",
        "recessive","general","heterozygote"),multiarg=FALSE)
    
    modelCode <- c(additive=1,dominant=2,recessive=3,general=4,heterozygote=5)
    
    # 1. Prepare SNPTEST files (from PLINK)
    # 1a. Create a SNPTEST workspace in Rtmp
    #stwork <- file.path(tempdir(),paste0("snptest_workspace_",.randomString()))
    #if (!dir.exists(stwork))
    #    dir.create(stwork,recursive=TRUE)
    wspace <- .validateWorkspacePath(wspace,"snptest")
    # 1b. obj to PLINK through snpStats::write.plink
    # 1c. Necessary phenotypes and covariates to SNPTEST sample file
    prepList <- .preparePlinkInputForSnptest(obj,response,covariates,pcs,wspace)
    
    # 2. Run the SNPTEST command
    disp("\nPerforming GWAS with SNPTEST")
    disp("Trait(s)    : ",response)
    disp("Covariate(s): ",paste0(covariates,collapse=", "))
    disp("Use PCs     : ",ifelse(pcs,"Yes","No"))
    
    snptest <- .getToolPath("snptest")
    command <- paste(
       paste0(snptest," \\"),
       paste0("  -data ",prepList$plink," \\"),
       paste0("    ",prepList$sample," \\"),
       paste0("  -frequentist ",modelCode[model]," \\"),
       "  -method score \\",
       paste0("  -pheno ",response," \\"),
       "  -cov_all \\",
       paste0("  -o ",file.path(wspace,"snptest.out")),
       sep="\n"
    )
    
    message("\nExecuting:\n",command)
    #out <- tryCatch(system(command,ignore.stdout=TRUE,ignore.stderr=TRUE),
    #    error=function(e) {
    #    message("Caught error: ",e$message)
    #    return(1L)
    #},finally="")
    
    out <- tryCatch({
        log <- .formatSystemOutputForDisp(capture.output({
            system(command,intern=TRUE)
        }))
        disp("\nSNPTEST output is:\n",level="full")
        disp(paste(log,collapse="\n"),"\n",level="full")
    },error=function(e) {
        message("Caught error: ",e$message)
        return(1L)
    },finally="")
    
    if (!.isEmpty(out) && out == 1L) {
        message("SNPTEST run failed! Will return unit p-values...")
        return(rep(1,nrow(obj)))
    }
    
    # 3. Read and gather SNPTEST output in a table similar to others
    tmp <- read.table(file.path(wspace,"snptest.out"),header=TRUE,sep=" ")
    res <- tmp[,c(23,24,21)]
    colnames(res) <- c("effect","se","pvalue")
    rownames(res) <- rownames(obj)
    
    disp("Found ",length(which(res[,3]<psig))," associations (uncorrected ",
        "p-values).")
    return(res)
}

gwaPlink <- function(obj,response,covariates=NULL,pcs=TRUE,psig=0.05,
    seed=42,cliopts=NULL,wspace=NULL,rc=NULL) {
    # Tool availability before all else
    if (!.toolAvailable("plink"))
        stop("PLINK program not found in the system!")
    .canRunGwa(obj)
    
    # Type of trait
    useLogistic <- .maybeBinaryForBinomial(phenotypes(obj)[,response])
    
    wspace <- .validateWorkspacePath(wspace,"plink")
    
    prepList <- .preparePlinkInputForPlink(obj,response,covariates,pcs,wspace)
    
    disp("\nPerforming GWAS with PLINK")
    disp("Trait(s)    : ",response)
    disp("Covariate(s): ",paste0(covariates,collapse=", "))
    disp("Use PCs     : ",ifelse(pcs,"Yes","No"))
    
    # Do we need to run permutations?
    perm <- ifelse(ncol(obj)<1000,TRUE,FALSE)
    # Threads?
    threads <- .coresFrac(rc)
    outBase <- file.path(wspace,paste0("out_",prepList$runid))
    plink <- .getToolPath("plink")
    command <- paste(
       paste0(plink," \\"),
       paste0("  --bfile ",prepList$bfile," \\"),
       paste0("  --pheno ",prepList$pheno," \\"),
       "  --allow-no-sex \\",
       #"  --keep-allele-order \\",
       "  --ci 0.95 \\",
       paste0("  --seed ",seed," \\"),
       paste0("  --out ",outBase),
       sep="\n"
    )
    
    # Rest arguments
    if (!is.null(prepList$covar)) {
        cvr <- paste0("  --covar ",prepList$covar)
        command <- paste0(command," \\\n",cvr)
    }
    
    model <- ifelse(is.null(prepList$covar)," --assoc",
        ifelse(useLogistic,"  --logistic beta","  --linear"))
    if (perm)
        model <- paste0(model," perm")
    command <- paste0(command," \\\n",model)
    
    if (threads > 1) {
        addThreads <- paste0("  --threads ",threads)
        command <- paste0(command," \\\n",addThreads)
    }
    
    if (!is.null(cliopts)) {
        addOpts <- paste0("  ",cliopts)
        command <- paste0(command," \\\n",addOpts)
    }
    
    message("\nExecuting:\n",command)
    out <- tryCatch({
        system(command,ignore.stdout=TRUE,ignore.stderr=TRUE)
        ## Does not look very good in R cli...
        #log <- .formatSystemOutputForDisp(capture.output({
        #    system(command,intern=TRUE)
        #}))
        if (prismaVerbosity() == "full") {
            logfile <- paste0(outBase,".log")
            log <- .formatSystemOutputForDisp(readLines(logfile))
            disp("\nPLINK output is:\n")
            disp(paste(log,collapse="\n"),"\n")
        }
    },error=function(e) {
        message("Caught error: ",e$message)
        return(1L)
    },finally="")
    
    if (!.isEmpty(out) && out == 1L) {
        message("PLINK run failed! Will return unit p-values...")
        return(rep(1,nrow(obj)))
    }
    
    # Read the results
    suffix <- ifelse(useLogistic,"logistic","linear")
    tmp <- read.table(paste0(outBase,".assoc.",suffix),header=TRUE,sep="")
    tmp <- tmp[!(tmp$TEST %in% prepList$covs),,drop=FALSE]
    res <- tmp[,c("BETA","SE","STAT","P")]
    if (perm) {
        ptmp <- read.table(paste0(outBase,".assoc.",suffix,".perm"),header=TRUE,
            sep="")
        res$P <- ptmp$EMP1
    }
    colnames(res) <- c("effect","se","stat","pvalue")
    rownames(res) <- rownames(obj)
    
    disp("Found ",length(which(res[,4]<psig))," associations (uncorrected ",
        "p-values).")
    return(res)
}

.preparePlinkInputForSnptest <- function(obj,response,covariates,pcs,wspace) {
    disp("Preparing SNPTEST association run at workspace directory: ",wspace)
    
    # Subset phenotypes and add PCs if requested
    p <- phenotypes(obj)
    chResCov <- .validateResponseAndCovariates(p,response,covariates)
    response <- chResCov$res
    covariates <- chResCov$cvs
    p <- p[,c(response,covariates),drop=FALSE]
    ct <- .initSnptestSampleFirstRow(p,response,covariates)
    if (pcs) { # Include robust PCs in the model
        if (.hasPcaCovariates(obj)) {
            pcov <- pcaCovariates(obj)
            p <- cbind(p,pcov)
            ct <- c(ct,rep("C",ncol(pcov)))
        }
        else
            warning("PC covariates requested in the model, but no calculated ",
                "PC covariates found! Ignoring...",immediate.=TRUE)
    }
    
    # Random run id
    runId <- .randomString()
    
    # Modify the phenotypes data frame
    disp(" writing SNPTEST sample file in ",wspace)
    p <- .checkAndCorrectFactorsFormat(p)
    p <- data.frame(sample_id=rownames(p),p)
    ct <- c("0",ct)
    sfile <- file.path(wspace,"pheno.sample")
    r1 <- paste0(colnames(p),collapse=" ")
    r2 <- paste0(ct,collapse=" ")
    writeLines(c(r1,r2),sfile)
    write.table(p,file=sfile,append=TRUE,quote=FALSE,sep=" ",row.names=FALSE,
        col.names=FALSE)
    
    # The plink files
    disp("  writing PLINK files in ",wspace)
    base <- file.path(wspace,paste0("plink_snptest_",runId))
    writePlink(obj,response,outBase=base)
    # FIXME: We reverse until find out wtf is going on...
    # But PLINK always returns specific effects no matter what
    #writePlink(obj,response,outBase=base,reverse=TRUE)
    
    disp("Preparation done!")
    return(list(plink=paste0(base,".bed"),sample=sfile,runid=runId))
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

.preparePlinkInputForPlink <- function(obj,response,covariates,pcs,wspace) {
    disp("Preparing PLINK association run at workspace directory: ",wspace)
    
    # Subset phenotypes and add PCs if requested
    p <- phenotypes(obj)
    chResCov <- .validateResponseAndCovariates(p,response,covariates)
    response <- chResCov$res
    covariates <- chResCov$cvs
    
    # Do we need a covariates file?
    needCov <- TRUE
    if (is.null(covariates) && !pcs)
        needCov <- FALSE
    
    # What PLINK expects
    pheno <- p[,c("FID","IID",response)]
    if (needCov)
        covar <- p[,c("FID","IID",covariates)]
        
    # Make sure FID and IID matches, otherwise cannot be properly handled by
    # snpStats::write.plink
    if (!identical(pheno[,"FID"],pheno[,"IID"])) {
        pheno[,"IID"] <- pheno[,"FID"]
        if (needCov)
            covar[,"IID"] <- covar[,"FID"]
    }
    
    # PCs needed? Will not run if needCov is FALSE anyway as pcs is FALSE
    if (pcs) {
        if (.hasPcaCovariates(obj)) {
            pcov <- pcaCovariates(obj)
            covar <- cbind(covar,pcov)
        }
        else
             warning("PC covariates reqeusted in the model, but no calculated ",
                "PC covariates found! Ignoring...",immediate.=TRUE)
    }
    
    # Random run id
    runId <- .randomString()
    
    # Prepare the phenotype - a file must be written with space delimited
    disp("  writing the PLINK phenotype file...")
    phenoFile <- file.path(wspace,paste0("pheno_",runId))
    pheno <- .checkAndCorrectFactorsFormat(pheno)
    write.table(pheno,file=phenoFile,quote=FALSE,row.names=FALSE)
    
    # Prepare the covariates - a file must be written with space delimited
    if (needCov) {
        disp("  writing the covariates file...")
        covarFile <- file.path(wspace,paste0("covar_",runId))
        covar <- .checkAndCorrectFactorsFormat(covar)
        write.table(covar,file=covarFile,quote=FALSE,row.names=FALSE)
    }
    else
        covarFile <- NULL
    
    # The plink files
    disp("  writing PLINK genotype files in ",wspace)
    base <- file.path(wspace,paste0("plink_plink_",runId))
    writePlink(obj,response,outBase=base)
    
    disp("Preparation done!")
    return(list(bfile=base,pheno=phenoFile,covar=covarFile,runid=runId,
        covs=names(covar)[3:ncol(covar)]))
}

.prepareGenotypesForBlup <- function(obj) {
    tmp <- genotypes(obj)
    geno <- as(tmp,"numeric")
    geno[geno==0] <- -1
    geno[geno==1] <- 0
    geno[geno==2] <- 1
    #geno[geno==0] <- 1
    #geno[geno==1] <- 0
    #geno[geno==2] <- -1
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
    
    if (res %in% cvs)
        stop("Response variable ",res," cannot be also a covariate!")
    
    return(list(res=res,cvs=cvs))
}

.validateBinaryForBinomial <- function(x) {
    if (length(unique(x)) > 2)
        stop("The response variable cannot have more than two values when ",
            "GLM family is 'binomial'!")
    if (!is.factor(x)) {
        if (!all(unique(x) %in% c(0,1))) {
            warning("When GLM family is 'binomial', the response variable ",
                "must be either 0, 1 or a 2-level factor! Converting to ",
                "factor...")
            x <- as.factor(x)
        }
    }
    return(x)
}

# 0-1 or 2-level factor
.isBinaryForBinomial <- function(x) {
    if (is.factor(x) && length(levels(x)) == 2)
        return(TRUE)
    if (all(unique(x) %in% c(0,1)))
        return(TRUE)
    return(FALSE)
}

.maybeBinaryForBinomial <- function(x) {
    return(tryCatch({
        suppressWarnings(.validateBinaryForBinomial(x))
        TRUE
    },error=function(e) { 
        return(FALSE)
    },finally=""))
}

.canRunGwa <- function(o) {
    # Object class?
    if (!is(o,"GWASExperiment"))
        stop("GWAS input must be an object of class GWASExperiment!")
        
    gv <- all(dim(genotypes(o)) > 0)
    pv <- !is.null(phenotypes(o)) && all(dim(phenotypes(o)) > 0)
    if (!gv)
        stop("Cannot run GWAS with a GWASExperiment without genotypes!")
    if (!pv)
        stop("Cannot run GWAS with a GWASExperiment without phenotypes!")
}


.checkGwaArgs <- function(args,test=c("glm","rrblup","statgen","snptest",
    "plink")) {
    test <- test[1]
    
    # Allowed and given values
    defaults <- getDefaults(test)
    allowed <- names(defaults)
    given <- names(args)
    
    # Check if illegal filter names have been provided
    check <- given %in% allowed
    if (!any(check))
        stop("No valid ",test," parameter name found!")
    if (!all(check)) {
        warning("The following ",test," parameters names are invalid and will ",
            "be ignored:\n",paste(given[!check],collapse=", "))
        args <- args[given[check]]
    }
    
    # Proceed with per algorithm check, will stop here if something wrong
    switch(test,
        glm = {
            .checkGlmArgs(args)
        },
        rrblup = {
            .checkRrblupArgs(args)
        },
        statgen = {
            .checkStatgenArgs(args)
        },
        snptest = {
            .checkSnptestArgs(args)
        },
        plink = {
            .checkPlinkArgs(args)
        }
    )
    
    # If all OK
    return(.setArg(defaults,args))
}

.checkGlmArgs <- function(a) {
    if (!.isEmpty(a$family))
        .checkTextArgs("Regression family (family)",a$family,c("gaussian",
            "binomial","poisson"),multiarg=FALSE)
    if (!.isEmpty(a$size))
        .checkNumArgs("SNPs to bucket test (size)",a$size, "numeric",0,"gt")
}

.checkRrblupArgs <- function(a) {
    if (!.isEmpty(a$pcblup))
        .checkTextArgs("PCs for rrBLUP (pcblup)",a$pcblup,c("auto","estim",
            "rrint","fixed","none"),multiarg=FALSE)
    if (!.isEmpty(a$npcs))
        .checkNumArgs("Number of PCs for rrBLUP (npcs)",a$npcs,"numeric",0,
            "gt")
}

.checkStatgenArgs <- function(a) {
}

.checkPlinkArgs <- function(a) {
    if (!.isEmpty(a$effect))
        .checkTextArgs("Genotype effect model for PLINK (effect)",a$effect,
            c("genotypic","hethom","dominant","recessive"),multiarg=FALSE)
    if (!.isEmpty(a$seed))
        .checkNumArgs("Seed for permutations for PLINK (seed)",a$seed,
            "numeric",0,"gt")
    if (!.isEmpty(a$workspace)) {
        if (!is.character(a$workspace))
            stop("The PLINK workspace directory must be a character!")
    }
    if (!.isEmpty(a$cliopts)) {
        if (!is.character(a$cliopts))
            stop("PLINK additional parameters must be a string!")
    }
}

.checkSnptestArgs <- function(a) {
    if (!.isEmpty(a$snptest))
        .checkTextArgs("Test class for SNPTEST (snptest)",a$snptest,
            c("frequentist","bayesian"),multiarg=FALSE)
    if (!.isEmpty(a$snpmodel))
        .checkTextArgs("Genotype model for SNPTEST (snpmodel)",a$snpmodel,
            c("additive","dominant","recessive","general","heterozygote"),
            multiarg=FALSE)
    if (!.isEmpty(a$workspace)) {
        if (!is.character(a$workspace))
            stop("The SNPTEST workspace directory must be a character!")
    }
}
