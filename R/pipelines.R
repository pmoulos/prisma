.prsPipelineDenovo <- function(gwe,phenotype,covariates,pcs,npcs,trainSize,
    niter,filters,pcaMethod,imputeMissing,imputeMethod,gwaMethods,gwaCombine,
    glmOpts,rrblupOpts,statgenOpts,snptestOpts,plinkOpts,prsMethods,
    lassosumOpts,prsiceOpts,prsWorkspace,logging,output,runId,dig,rc) {
    # Do we continue a crashed/interrupted run?
    if (length(niter) > 1)
        iters <- niter
    else
        iters <- seq_len(niter)
    
    # Get the GDS file for copying to subfolders if needed
    gorig <- gdsfile(gwe)
    
    # The worker
    theResult <- cmclapply(iters,function(i) {
        pad <- paste0(rep("0",dig - nchar(as.character(i))),collapse="")
        #iterWspace <- file.path(prsWorkspace,paste0(format(Sys.time(),
        #    "%Y%m%d%H%M%S"),"_denovo_",pad,i))
        iterWspace <- file.path(prsWorkspace,paste0(runId,"_denovo_",pad,i))
        dir.create(iterWspace,recursive=TRUE,showWarnings=FALSE)
        
        # Record progress
        .recordProgress(i,iterWspace)
        
        # If not provided by the user, PLINK and SNPTEST workspaces should live
        # within the main PRS analysis workspace
        if ("snptest" %in% gwaMethods && is.null(snptestOpts$workspace)) {
            snptestOpts$workspace <- file.path(iterWspace,
                paste0("snptest_",.randomString()))
            dir.create(snptestOpts$workspace,recursive=TRUE,showWarnings=FALSE)
        }
        if ("plink" %in% gwaMethods && is.null(plinkOpts$workspace)) {
            plinkOpts$workspace <- file.path(iterWspace,
                paste0("plink_",.randomString()))
            dir.create(plinkOpts$workspace,recursive=TRUE,showWarnings=FALSE)
        }
        
        # Local result saving for restoring after possible crash
        saveFile <- file.path(iterWspace,"denovo_result.RData")
        
        if (logging == "file") {
            # Print something basic before sinking
            ilf <- file.path(prsWorkspace,paste0(runId,"_iter_",pad,i,".log"))
            disp("Executing denovo pipeline iteration ",i,", log is at ",ilf)
            fh <- file(ilf,open="wt")
            sink(fh,type="output")
            sink(fh,type="message")
        }
        else
            disp("\n")
        
        disp("==============================================================")
        disp("-----> Denovo pipeline iteration ",i)
        disp("==============================================================\n")
        
        # Partition the object
        splitResult <- .denovoDatasetPartition(gwe,phenotype,trainSize,output,
            logging)
        theTrain <- splitResult$train
        theTest <- splitResult$test
        trainIndex <- splitResult$itrain
        testIndex <- splitResult$itest
        
        # Copy the gds file in the iteration workspace and change training and
        # test objects metadata, if running in parallel
        if (!is.null(rc)) {
            gdest <- file.path(iterWspace,basename(gorig))
            file.copy(from=gorig,to=gdest,overwrite=TRUE)
            gdsfile(theTrain) <- gdsfile(theTest) <- gdest
        }
                
        # Base and target QC
        qcResult <- .localQc(theTrain,theTest,filters,imputeMissing,pcs,
            pcaMethod,npcs,logging)
        theTrain <- qcResult$base
        theTest <- qcResult$target
        
        # Run GWAS on base
        theTrain <- .gwaForPrsWorker(theTrain,phenotype,covariates,pcs,
            gwaMethods,gwaCombine,glmOpts,rrblupOpts,statgenOpts,snptestOpts,
            plinkOpts,rc,logging)
        
        # Run PRS
        theTest <- .runPRSWorkerCreate(theTrain,theTest,phenotype,covariates,
            pcs,prsMethods,prsiceOpts,iterWspace,rc,logging)
        
        # The theTest variable will be returned if output is "gwaslist", we
        # need to supply the original effects and p-values
        eff <- pv <- matrix(0,nrow(theTest),length(gwaMethods))
        rownames(eff) <- rownames(pv) <- rownames(theTest)
        colnames(eff) <- colnames(pv) <- colnames(pvalues(theTrain))
        sf <- intersect(rownames(theTrain),rownames(theTest))
        eff[sf,] <- effects(theTrain)[sf,,drop=FALSE]
        pv[sf,] <- pvalues(theTrain)[sf,,drop=FALSE]
        pcno <- ifelse(pcs,ifelse(.hasPcaCovariates(theTest),
            ncol(pcaCovariates(theTest)),0),0)
        effects(theTest,phenotype,covariates,pcno) <- eff
        pvalues(theTest,phenotype,covariates,pcno) <- pv
        
        # A temporary attribute to the GWASExperiment output to maintain the
        # workspace directory for evaluation. Will be removed if super clean
        attr(theTest,"workspace") <- iterWspace
        # And also keep the PRSice R2 triplet
        if ("prsice" %in% prsMethods)
            attr(theTest,"PR2") <- .readPrsiceR2(iterWspace)
        
        if (output == "gwaslist") {
            # Restore original GDS
            if (!is.null(rc))
                gdsfile(theTest) <- gorig
            result <- theTest
        }
        else if (output == "summaries")
            result <- list(
                baseIndex=trainIndex,
                targetIndex=testIndex,
                effects=effects(theTest),
                betas=prsbetas(theTest),
                npcs=pcno,
                pr2=attr(theTest,"PR2")
            )
        
        # Remove copied GDS
        if (!is.null(rc))
            unlink(gdest,force=TRUE)
        
        disp("\n----- Saving iteration ",i," result to ",saveFile,"-----\n")
        save(result,file=saveFile)
        
        disp("==============================================================")
        
        if (logging == "file") {
            sink(type="message")
            sink(type="output")
        }
        else
            disp("\n")
        
        # Record progress
        .recordProgress(i,iterWspace,TRUE)
        
        return(result)
    },rc=rc,setseed=TRUE)
    
    return(theResult)
}

.prsPipelineExternal <- function(gwe,phenotype,covariates,pcs,npcs,snpSelection,
    trainSize,niter,filters,pcaMethod,imputeMissing,imputeMethod,gwaMethod,
    gwaOpts,prsiceOpts,prsWorkspace,logging,runId,dig,rc) {
    # Do we continue a crashed/interrupted run?
    if (length(niter) > 1)
        iters <- niter
    else
        iters <- seq_len(niter)
        
    # Get the GDS file for copying to subfolders if needed
    gorig <- gdsfile(gwe)
        
    # The actual worker
    theResult <- cmclapply(iters,function(i) {
        pad <- paste0(rep("0",dig - nchar(as.character(i))),collapse="")
        #iterWspace <- file.path(prsWorkspace,paste0(format(Sys.time(),
        #    "%Y%m%d%H%M%S"),"_external_",pad,i))
        iterWspace <- file.path(prsWorkspace,paste0(runId,"_external_",pad,i))
        dir.create(iterWspace,recursive=TRUE,showWarnings=FALSE)
        
        # Record progress
        .recordProgress(i,iterWspace)
        
        # If not provided by the user, PLINK and SNPTEST workspaces should live
        # within the main PRS analysis workspace
        if (gwaMethod %in% c("snptest","plink") && is.null(gwaOpts$workspace)) {
            gwaOpts$workspace <- file.path(iterWspace,
                paste0(gwaMethod,"_",.randomString()))
            dir.create(gwaOpts$workspace,recursive=TRUE,showWarnings=FALSE)
        }
        
        # Local result saving for restoring after possible crash
        saveFile <- file.path(iterWspace,"external_result.RData")
        
        if (logging == "file") {
            # Print something basic before sinking
            ilf <- file.path(prsWorkspace,paste0(runId,"_iter_",pad,i,".log"))
            disp("Executing external pipeline iteration ",i,", log is at ",ilf)
            fh <- file(ilf,open="wt")
            sink(fh,type="output")
            sink(fh,type="message")
        }
        else
            disp("\n")
        
        disp("==============================================================")
        disp("-----> External pipeline iteration ",i)
        disp("==============================================================\n")
        
        # Partition the object
        splitRes <- .externalDatasetPartition(gwe,phenotype,trainSize,logging)
        base <- splitRes$base
        target <- splitRes$target
        
        # Copy the gds file in the iteration workspace and change training and
        # test objects metadata, if running in parallel
        if (!is.null(rc)) {
            gdest <- file.path(iterWspace,basename(gorig))
            file.copy(from=gorig,to=gdest,overwrite=TRUE)
            gdsfile(base) <- gdsfile(target) <- gdest
        }
        
        # Base and target QC
        qcResult <- .localQc(base,target,filters,imputeMissing,pcs,
            pcaMethod,npcs,logging)
        base <- qcResult$base
        target <- qcResult$target
        
        # Run basic GWAS on base to get coefficients (from one method for now)
        base <- .gwaForPrsReducedWorker(base,phenotype,covariates,pcs,gwaMethod,
            gwaOpts,rc,logging)
            
        # If snpSelection is a data.frame, use its rownames
        if (is.data.frame(snpSelection)) {}
            snpSelection <- rownames(snpSelection)
        
        nsnp <- length(snpSelection)
        
        # Run PRS with a safeguard for possibly filtered SNPs
        safeguard <- Reduce("intersect",list(snpSelection,rownames(base),
            rownames(target)))
        if (length(safeguard) < length(snpSelection))
            disp(length(snpSelection) - length(safeguard)," SNPs not found in ",
                "base data! Excluding...")
        else
            safeguard <- snpSelection
        subBase <- base[safeguard,,drop=FALSE]
        #subBase <- base
        subTarget <- target[safeguard,,drop=FALSE]

        if (pcs) {
            disp("Recalculating PCs...")
            pcno <- ncol(pcaCovariates(base))
            subBase <- calcPcaCovar(subBase,method=pcaMethod,npc=pcno)
            subTarget <- calcPcaCovar(subTarget,method=pcaMethod,npc=pcno)
        }
        
        # This approach may not work well, as some algorithms are potentially
        # using SNP relationships internally (e.g. SNPTEST, or statgenGWAS)
        ## Run basic GWAS on base to get coefficients (from one method for now)
        #subBase <- .gwaForPrsReducedWorker(subBase,phenotype,covariates,pcs,
        #   gwaMethod,gwaOpts,rc,logging)
        
        result <- .runPRSWorkerApply(subBase,subTarget,phenotype,covariates,
            pcs,iterWspace,rc,logging)
        
        # A temporary attribute to the GWASExperiment output to maintain the
        # workspace directory for evaluation. Will be removed if super clean
        attr(result,"workspace") <- iterWspace
        # And also keep the PRSice R2 triplet
        attr(result,"PR2") <- .readPrsiceR2(iterWspace)
        # And the number of SNPs used
        attr(result,"nsnp") <- nsnp
        
        # Remove copied GDS
        if (!is.null(rc))
            unlink(gdest,force=TRUE)
        
        disp("\n----- Saving iteration ",i," result to ",saveFile,"-----\n")
        save(result,file=saveFile)
        
        disp("==============================================================")
        
        if (logging == "file") {
            sink(type="message")
            sink(type="output")
        }
        else
            disp("\n")
        
        # Record progress
        .recordProgress(i,iterWspace,TRUE)
        
        return(result)
    },rc=rc,setseed=TRUE)
    
    return(theResult)
}

.prsPipelineValidate <- function(dnList,gwe,phenotype,covariates,pcs,
    snpSelection,niter,pcaMethod,gwaMethod,prsWorkspace,logging,runId,dig,rc) {
    # Do we continue a crashed/interrupted run?
    if (length(niter) > 1)
        iters <- niter
    else
        iters <- seq_len(niter)
        
    # Get the GDS file for copying to subfolders if needed
    gorig <- gdsfile(gwe)
        
    # The actual worker
    theResult <- cmclapply(iters,function(i,D) {
        pad <- paste0(rep("0",dig - nchar(as.character(i))),collapse="")
        iterWspace <- file.path(prsWorkspace,paste0(runId,"_evaluate_",pad,i))
        dir.create(iterWspace,recursive=TRUE,showWarnings=FALSE)
        
        # Record progress
        .recordProgress(i,iterWspace)
        
        # Local result saving for restoring after possible crash
        saveFile <- file.path(iterWspace,"evaluate_result.RData")
        
        if (logging == "file") {
            # Print something basic before sinking
            ilf <- file.path(prsWorkspace,paste0(runId,"_iter_",pad,i,".log"))
            disp("Executing evaluation pipeline iteration ",i,", log is at ",
                ilf)
            fh <- file(ilf,open="wt")
            sink(fh,type="output")
            sink(fh,type="message")
        }
        else
            disp("\n")
        
        disp("==============================================================")
        disp("-----> Evaluation pipeline iteration ",i)
        disp("==============================================================\n")
        
        # Run PRS with a safeguard for possibly filtered SNPs
        dnTmp <- D[[i]]
        
        # If dnList is a list of GWASExperiments/targets - OK, else we must
        # retrieve the genotypes from the total object.
        if (!is(dnTmp,"GWASExperiment")) { # Should be a summary list
            dnObj <- gwe[,dnTmp$targetIndex,drop=FALSE]
            effs <- dnTmp$effects
            betas <- dnTmp$betas
            pcno <- dnTmp$npcs
            if (is.null(pcno))
                pcno <- ifelse(.hasPcaCovariates(gwe),
                    ncol(pcaCovariates(gwe)),0)
            #colnames(betas)[colnames(betas)=="prsice"] <- gwaMethod
        }
        else { # Is already target (may have PCs)
            dnObj <- dnTmp
            effs <- effects(dnObj)
            betas <- prsbetas(dnObj)
            pcno <- ifelse(.hasPcaCovariates(dnObj),
                    ncol(pcaCovariates(dnObj)),0)
            #colnames(betas)[colnames(betas)=="prsice"] <- gwaMethod
        }
        
        # Copy the gds file in the iteration workspace and change training and
        # test objects metadata, if running in parallel
        if (!is.null(rc)) {
            gdest <- file.path(iterWspace,basename(gorig))
            file.copy(from=gorig,to=gdest,overwrite=TRUE)
            gdsfile(dnObj) <- gdest
        }
        
        # Check of snpSelection has been performed upstream
        hasAvgEffs <- FALSE
        if (is.data.frame(snpSelection)) {
            avgEffs <- snpSelection[,"effect",drop=FALSE]
            snpSelection <- rownames(snpSelection)
            hasAvgEffs <- TRUE
        }
        
        nsnp <- length(snpSelection)
        
        safeguard <- Reduce("intersect",list(snpSelection,rownames(dnObj),
            rownames(betas)))
        if (length(safeguard) < length(snpSelection))
            disp(length(snpSelection) - length(safeguard)," SNPs not found in ",
                "target data! Excluding...")
        else
            safeguard <- snpSelection
        
        # Construct the effects of the reduced object
        if (hasAvgEffs) {
            feff <- as.matrix(avgEffs[safeguard,,drop=FALSE])
            colnames(feff) <- gwaMethod
        }
        else
            feff <- effs[safeguard,gwaMethod,drop=FALSE]
            
        if ("lassosum" %in% colnames(betas)) {
            # Lassosum always has > SNPs than PRSice
            notInPrsice <- rownames(betas)[which(betas[,"prsice"] == 0)]
            inLasso <- rownames(betas)[which(betas[,"lassosum"] != 0)]
            toReplace <- intersect(safeguard,intersect(inLasso,notInPrsice))
            feff[toReplace,1] <- betas[toReplace,"lassosum"]
        }
        # Some zero effects may remain due to aggregation
        feff <- feff[feff[,1] != 0,,drop=FALSE]
        # The final safeguard
        subObj <- dnObj[rownames(feff),,drop=FALSE]
        
        if (pcs) {
            disp("Recalculating PCs...")
            subObj <- calcPcaCovar(subObj,method=pcaMethod,npc=pcno)
            pcno <- ncol(pcaCovariates(subObj))
        }
        
        effects(subObj,phenotype,covariates,pcno) <- feff
        
        result <- .runPRSWorkerApply(subObj,subObj,phenotype,covariates,
            pcs,iterWspace,rc,logging)
        
        # A temporary attribute to the GWASExperiment output to maintain the
        # workspace directory for evaluation. Will be removed if super clean
        attr(result,"workspace") <- iterWspace
        # And also keep the PRSice R2 triplet
        attr(result,"PR2") <- .readPrsiceR2(iterWspace)
        # And the number of SNPs used
        attr(result,"nsnp") <- nsnp
        
        # Remove copied GDS
        if (!is.null(rc))
            unlink(gdest,force=TRUE)
        
        disp("\n----- Saving iteration ",i," result to ",saveFile,"-----\n")
        save(result,file=saveFile)
        
        disp("==============================================================")
        
        if (logging == "file") {
            sink(type="message")
            sink(type="output")
        }
        else
            disp("\n")
        
        # Record progress
        .recordProgress(i,iterWspace,TRUE)
        
        return(result)
    },dnList,rc=rc,setseed=TRUE)
    
    return(theResult)
}

.prsPipelineEvaluate <- function(dnList,gwe,response,covariates,pcs,
    snpSelection,niter,pcaMethod,prsWorkspace,logging,runId,dig,rc=NULL) {
    # Do we continue a crashed/interrupted run?
    if (missing(niter))
        niter <- seq_along(dnList)
    
    if (length(niter) > 1)
        iters <- niter
    else
        iters <- seq_len(niter)
        
    # Get the GDS file for copying to subfolders if needed
    gorig <- gdsfile(gwe)
        
    # The actual worker
    theResult <- cmclapply(iters,function(i,D) {
        pad <- paste0(rep("0",dig - nchar(as.character(i))),collapse="")
        iterWspace <- file.path(prsWorkspace,paste0(runId,"_evaluate_",pad,i))
        dir.create(iterWspace,recursive=TRUE,showWarnings=FALSE)
        
        # Record progress
        .recordProgress(i,iterWspace)
        
        # Local result saving for restoring after possible crash
        saveFile <- file.path(iterWspace,"evaluate_result.RData")
        
        ## We only need a temp workspace in this case, just to copy the GDS
        #iterWspace <- file.path(prsWorkspace,paste0("_evaluate_",i))
        #dir.create(iterWspace,recursive=TRUE,showWarnings=FALSE)
        
        if (logging == "file") {
            # Print something basic before sinking
            ilf <- file.path(prsWorkspace,paste0(runId,"_iter_",pad,i,".log"))
            disp("Executing evaluation pipeline iteration ",i,", log is at ",
                ilf)
            fh <- file(ilf,open="wt")
            sink(fh,type="output")
            sink(fh,type="message")
        }
        else
            disp("\n")
        
        disp("==============================================================")
        disp("-----> Evaluation pipeline iteration ",i)
        disp("==============================================================\n")
        
        # Run PRS with a safeguard for possibly filtered SNPs
        dnTmp <- D[[i]]
        
        # If dnList is a list of GWASExperiments/targets - OK, else we must
        # retrieve the genotypes from the total object.
        if (!is(dnTmp,"GWASExperiment")) { # Should be a summary list
            dnObj <- gwe[,dnTmp$targetIndex,drop=FALSE]
            pcno <- dnTmp$npcs
            if (is.null(pcno))
                pcno <- ifelse(.hasPcaCovariates(gwe),
                    ncol(pcaCovariates(gwe)),0)
        }
        else { # Is already target (may have PCs)
            dnObj <- dnTmp
            pcno <- ifelse(.hasPcaCovariates(dnObj),
                    ncol(pcaCovariates(dnObj)),0)
        }
        
        # Copy the gds file in the iteration workspace and change training and
        # test objects metadata, if running in parallel
        if (!is.null(rc)) {
            gdest <- file.path(iterWspace,basename(gorig))
            file.copy(from=gorig,to=gdest,overwrite=TRUE)
            gdsfile(dnObj) <- gdest
        }
        
        safeguard <- intersect(rownames(snpSelection),rownames(dnObj))
        if (length(safeguard) < nrow(snpSelection))
            disp(nrow(snpSelection) - length(safeguard)," SNPs not found in ",
                "target data! Excluding...")
        else
            safeguard <- rownames(snpSelection)
        
        subObj <- dnObj[safeguard,,drop=FALSE]
        
        if (pcs) {
            disp("Recalculating PCs...")
            subObj <- calcPcaCovar(subObj,method=pcaMethod,npc=pcno)
            pcno <- ncol(pcaCovariates(subObj))
        }
        
        # The expected call is
        result <- tryCatch({
            prsRegressionMetrics(snpSelection,subObj,response,covariates,
                pcs,step=nrow(snpSelection))
        },error=function(e) {
            .exitFromSink(logging)
            stop("Caught error during target PRS evaluation: ",e$message,
                call.=FALSE)
        },interrupt=function(i) {
            .exitFromSink(logging)
            stop("Caught keyboard interruption during target PRS evaluation!",
                call.=FALSE)
        },finally="")
        
        # Remove copied GDS and its folder
        if (!is.null(rc))
            unlink(gdest,force=TRUE)
        
        disp("\n----- Saving iteration ",i," result to ",saveFile,"-----\n")
        save(result,file=saveFile)
        
        disp("==============================================================")
        
        if (logging == "file") {
            sink(type="message")
            sink(type="output")
        }
        else
            disp("\n")
        
        # Record progress
        .recordProgress(i,iterWspace,TRUE)
        
        return(result)
    },dnList,rc=rc,setseed=TRUE)
    
    ## Collapse the logs
    #logList <- lapply(seq_along(dnList),function(i) {
    #    lof <- file.path(prsWorkspace,paste0("_iter_",i,".log"))
    #    tmp <- readLines(lof)
    #    unlink(lof,force=TRUE)
    #    return(c(tmp,"\n"))
    #})
    #theLog <- do.call("c",logList)
    #writeLines(theLog,file.path(prsWorkspace,paste0("evaluation.log")))
    
    return(theResult)
}

.denovoDatasetPartition <- function(gwe,phenotype,trainSize,output,logging) {
    splitResult <- tryCatch({
        .denovoDatasetPartitionWorker(gwe,phenotype,trainSize,output)
    },error=function(e) {
        .exitFromSink(logging)
        stop("Caught error during dataset partitioning: ",e$message,
            call.=FALSE)
    },interrupt=function(i) {
        .exitFromSink(logging)
        stop("Caught keyboard interruption during dataset partitioning!",
            call.=FALSE)
    },finally="")
    return(splitResult)
}

.denovoDatasetPartitionWorker <- function(gwe,phenotype,trainSize,output) {
    disp("----- Dataset partitioning -----\n")
    if (output == "gwaslist") {
        tmp <- partitionGWAS(gwe,by=phenotype,n=1,frac=trainSize,
            out="ttboth")
        theTrain <- tmp$train
        theTest <- tmp$test
        trainIndex <- testIndex <- NULL
    }
    else if (output == "summaries") {
        tmp <- partitionGWAS(gwe,by=phenotype,n=1,frac=trainSize,
            out="index")
        trainIndex <- tmp[[1]]
        testIndex <- setdiff(seq_len(ncol(gwe)),trainIndex)
        theTrain <- gwe[,trainIndex,drop=FALSE]
        theTest <- gwe[,testIndex,drop=FALSE]
    }
    return(list(train=theTrain,test=theTest,itrain=trainIndex,itest=testIndex))
}

.localQc <- function(base,target,filters,imputeMissing,pcs,pcaMethod,
    npcs,logging) {
    qcResult <- tryCatch({
        .localQcWorker(base,target,filters,imputeMissing,pcs,pcaMethod,npcs)
    },error=function(e) {
        .exitFromSink(logging)
        stop("Caught error during dataset QC: ",e$message,call.=FALSE)
    },interrupt=function(i) {
        .exitFromSink(logging)
        stop("Caught keyboard interruption during dataset QC!",call.=FALSE)
    },finally="")
    return(qcResult)
}

.localQcWorker <- function(base,target,filters,imputeMissing,pcs,
    pcaMethod,npcs) {
    disp("\n----- Training (base) QC -----")
    base <- filterGWAS(base,filters=filters,imputeMissing=imputeMissing)
    if (pcs) {
        disp("\n----- Training (base) PCA -----\n")
        base <- calcPcaCovar(base,method=pcaMethod,npc=npcs)
        pcno <- ncol(pcaCovariates(base))
    }
            
    # Target QC
    disp("\n----- Test (target) QC -----")
    target <- filterGWAS(target,filters=filters,imputeMissing=imputeMissing)
    if (pcs) {
        disp("\n----- Test (target) PCA -----\n")
        target <- calcPcaCovar(target,method=pcaMethod,npc=pcno)
    }
    return(list(base=base,target=target))
}

.gwaForPrsWorker <- function(theTrain,phenotype,covariates,pcs,gwaMethods,
    gwaCombine,glmOpts,rrblupOpts,statgenOpts,snptestOpts,plinkOpts,rc,
    logging) {
    disp("\n----- Training (base) GWAS -----")
    theTrain <- tryCatch({
        gwa(theTrain,phenotype,covariates,pcs=pcs,methods=gwaMethods,
            combine=gwaCombine,glmOpts=glmOpts,rrblupOpts=rrblupOpts,
            statgenOpts=statgenOpts,snptestOpts=snptestOpts,
            plinkOpts=plinkOpts,rc=rc)
    },error=function(e) {
        .exitFromSink(logging)
        stop("Caught error during base GWA: ",e$message,call.=FALSE)
    },interrupt=function(i) {
        .exitFromSink(logging)
        stop("Caught keyboard interruption during base GWA!",call.=FALSE)
    },finally="")
    return(theTrain)
}

.runPRSWorkerCreate <- function(theTrain,theTest,phenotype,covariates,pcs,
    prsMethods,prsiceOpts,iterWspace,rc,logging) {
    disp("\n----- PRS analysis with base and target -----")
    theTest <- tryCatch({
        runPRS(base=theTrain,target=theTest,phenotype,covariates,pcs=pcs,
            methods=prsMethods,prsiceOpts=prsiceOpts,wspace=iterWspace,
            rc=rc)
    },error=function(e) {
        .exitFromSink(logging)
        stop("Caught error during target PRS: ",e$message,call.=FALSE)
    },interrupt=function(i) {
        .exitFromSink(logging)
        stop("Caught keyboard interruption during target PRS!",call.=FALSE)
    },finally="")
    return(theTest)
}

.externalDatasetPartition <- function(gwe,phenotype,trainSize,logging) {
    splitResult <- tryCatch({
        .externalDatasetPartitionWorker(gwe,phenotype,trainSize)
    },error=function(e) {
        .exitFromSink(logging)
        stop("Caught error during dataset partitioning: ",e$message,
            call.=FALSE)
    },interrupt=function(i) {
        .exitFromSink(logging)
        stop("Caught keyboard interruption during dataset partitioning!",
            call.=FALSE)
    },finally="")
    return(splitResult)
}

.externalDatasetPartitionWorker <- function(gwe,phenotype,trainSize) {
    disp("----- Dataset partitioning -----\n")
    tmp <- partitionGWAS(gwe,by=phenotype,n=1,frac=trainSize,
        out="ttboth")
    base <- tmp$train
    target <- tmp$test
    return(list(base=base,target=target))
}

.gwaForPrsReducedWorker <- function(base,phenotype,covariates,pcs,gwaMethod,
    gwaOpts,rc,logging) {
    disp("\n----- Base GWAS -----")
    base <- tryCatch({
        rcm <- paste0('gwa(base,phenotype,covariates,pcs=pcs,',
            'methods=gwaMethod,',paste0(gwaMethod,"Opts"),'=gwaOpts,rc=rc)')
        eval(parse(text=rcm))
    },error=function(e) {
        .exitFromSink(logging)
        stop("Caught error during reduced GWAS: ",e$message,call.=FALSE)
    },interrupt=function(i) {
        .exitFromSink(logging)
        stop("Caught keyboard interruption during reduced GWAS!",call.=FALSE)
    },finally="")
    return(base)
}

.runPRSWorkerApply <- function(base,target,phenotype,covariates,pcs,
    iterWspace,rc,logging) {
    disp("\n----- PRS analysis with base and target -----")
    prsiceOut <- tryCatch({
        prsicePRS(base=base,target=target,
            response=phenotype,covariates=covariates,pcs=pcs,
            mode="apply",wspace=iterWspace,rc=rc)
    },error=function(e) {
        .exitFromSink(logging)
        stop("Caught error during target PRS application: ",e$message,
            call.=FALSE)
    },interrupt=function(i) {
        .exitFromSink(logging)
        stop("Caught keyboard interruption during target PRS application!",
            call.=FALSE)
    },finally="")
    return(prsiceOut)
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

.recordProgress <- function(i,wspace,done=FALSE) {
    progFile <- file.path(wspace,".prog.json")
    if (!file.exists(progFile)) { # Initial write
        prog <- list(iter=i,done=done)
        write_json(prog,path=progFile,auto_unbox=TRUE,pretty=TRUE)
    }
    else { # Exists, update
        curr <- fromJSON(progFile)
        curr$iter <- i
        curr$done <- done
        write_json(curr,path=progFile,auto_unbox=TRUE,pretty=TRUE)
    }
}

.exitFromSink <- function(s) {
    if (s == "file") {
        sink(type="message")
        sink(type="output")
    }
}
