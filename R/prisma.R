# The prsPipeline, if fed with a filtered object (e.g. according to GWAS
# p-value cutoffs and the SNPs therein) can be used for PGS evaluation
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
    gwaMethods=c("glm","rrblup","statgen","snptest","plink"), # lasso later
    gwaCombine=c("fisher","simes","max","min","harmonic","whitlock","pandora"),
    glmOpts=getDefaults("glm"),
    rrblupOpts=getDefaults("rrblup"),
    statgenOpts=getDefaults("statgen"),
    snptestOpts=getDefaults("snptest"),
    plinkOpts=getDefaults("plink"),
    prsMethods=c("lassosum","prsice"),
    lassosumOpts=getDefaults("lassosum"),
    prsiceOpts=getDefaults("prsice"),
    prsWorkspace=NULL,
    cleanup=c("none","intermediate","all"),
    logging=c("screen","file"),
    output=c("gwaslist","summaries"),
    rc=NULL
) {
    #TODO: Log options in effect (or all), like in metaseqR
    #TODO: Unique run identifier to append to RData files
    
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
        c("glm","rrblup","statgen","snptest","plink"),multiarg=TRUE)
    .checkTextArgs("p-value combination (gwaCombine)",gwaCombine,
        c("fisher","simes","max","min","harmonic","whitlock","pandora"),
        multiarg=FALSE)
    .checkTextArgs("Logging option",logging,c("screen","file"),multiarg=FALSE)
    .checkTextArgs("Cleanup option",cleanup,c("none","intermediate","all"),
        multiarg=FALSE)
    .checkTextArgs("Output option",output,c("gwaslist","summaries"),
        multiarg=FALSE)
    
    glmOpts <- .checkGwaArgs(glmOpts,"glm")
    rrblupOpts <- .checkGwaArgs(rrblupOpts,"rrblup")
    #statgenOpts <- .checkGwaArgs(statgenOpts,"statgen")
    snptesOpts <- .checkGwaArgs(snptestOpts,"snptest")
    plinkOpts <- .checkGwaArgs(plinkOpts,"plink")
    lassosumOpts <- .checkPrsArgs(lassosumOpts,"lassosum")
    prsiceOpts <- .checkPrsArgs(prsiceOpts,"prsice")
    
    # .prettyLogOptions(...)
    
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
    
    ## Copy the gds file in the master workspace. Keep in mind that in the
    ## future, if we parallelize the iterations, the GDS file should be copied
    ## to each workspace subdir
    #m <- metadata(gwe)
    #gorig <- m$gdsfile
    #gdest <- file.path(prsWorkspace,basename(gorig))
    #file.copy(from=gorig,to=gdest,overwrite=TRUE)
    #m$gdsfile <- gdest
    #metadata(gwe) <- m
    
    # Main iteration
    if (is.null(snpSelection))
        theResult <- .prsPipelineDenovo(gwe,phenotype,covariates,pcs,npcs,
            trainSize,niter,filters,pcaMethod,imputeMissing,imputeMethod,
            gwaMethods,gwaCombine,glmOpts,rrblupOpts,statgenOpts,snptestOpts,
            plinkOpts,prsMethods,lassosumOpts,prsiceOpts,prsWorkspace,logging,
            output,rc)
    else
        # This time is a PGS Catalog data frame
        theResult <- .prsPipelineExternal(gwe,phenotype,covariates,pcs,npcs,
            snpSelection,trainSize,niter,filters,pcaMethod,imputeMissing,
            imputeMethod,gwaMethods[1],gwaOpts,prsiceOpts,prsWorkspace,logging,
            output,rc)

    switch(cleanup,
        none = {
            disp("All pipeline output can be found at ",prsWorkspace)
        },
        intermediate = {
            disp("Cleaning up temporary program-specific intermediate files ",
                "from ",prsWorkspace)
            .partialCleanup(prsWorkspace)
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
    
    ## If output is gwaslist, restore the original GDS file location
    #if (output == "gwalist")
    #    theResult <- lapply(theResult,function(x,g) {
    #        m <- metadata(x)
    #        m$gdsfile <- g
    #        metadata(x) <- m
    #        return(x)
    #    },gorig)
    ## Delete the copy from the workspace
    #unlink(gdest,force=TRUE)
    
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
        return(Reduce(m,g))
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
    glmOpts,rrblupOpts,statgenOpts,snptestOpts,plinkOpts,prsMethods,
    lassosumOpts,prsiceOpts,prsWorkspace,logging,output,rc) {
    # Initialize the list of GWASExperiment s
    theResult <- vector("list",niter)
    
    # For output directory (workspace) format
    dig <- nchar(as.character(niter))
    
    # Do we continue a crashed/interrupted run?
    if (length(niter) > 1)
        iters <- niter
    else
        iters <- seq_len(niter)
    
    # The RData file to save iteratively the growing result oject
    saveFile <- file.path(prsWorkspace,"denovo_iter_data.RData")
    
    # Get the GDS file for copying to subfolders
    gorig <- metadata(gwe)$gdsfile
    
    # The worker
    #for (i in iters) {
    theResult <- cmclapply(iters,function(i) {
        pad <- paste0(rep("0",dig - nchar(as.character(i))),collapse="")
        iterWspace <- file.path(prsWorkspace,paste0(format(Sys.time(),
            "%Y%m%d%H%M%S"),"_denovo_",pad,i))
        dir.create(iterWspace,recursive=TRUE,showWarnings=FALSE)
        
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
            disp("Executing denovo pipeline iteration ",i)
            fh <- file(file.path(prsWorkspace,paste0("iter_",pad,i,".log")),
                open="wt")
            sink(fh,type="output")
            sink(fh,type="message")
        }
        
        disp("\n==============================================================")
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
        # test objects metadata
        gdest <- file.path(iterWspace,basename(gorig))
        file.copy(from=gorig,to=gdest,overwrite=TRUE)
        gdsfile(theTrain) <- gdsfile(theTest) <- gdest
        
        # Base QC
        qcResult <- .denovoQc(theTrain,theTest,filters,imputeMissing,pcs,
            pcaMethod,npcs,logging)
        theTrain <- qcResult$train
        theTest <- qcResult$test
        
        # Run GWAS on base
        theTrain <- .gwaForPrsWorker(theTrain,phenotype,covariates,pcs,
            gwaMethods,gwaCombine,glmOpts,rrblupOpts,statgenOpts,snptestOpts,
            plinkOpts,rc,logging)
        
        # Run PRS
        theTest <- .runPRSWorker(theTrain,theTest,phenotype,covariates,pcs,
            prsMethods,prsiceOpts,iterWspace,rc,logging)
        
        # A temporary attribute to the GWASExperiment output to maintain the
        # workspace directory for evaluation. Will be removed if super clean
        attr(theTest,"workspace") <- iterWspace
        # And also keep the PRSice R2 triplet
        if ("prsice" %in% prsMethods)
            attr(theTest,"PR2") <- .readPrsiceR2(iterWspace)
        
        if (output == "gwaslist") { # Restore original GDS
            gdsfile(theTest) <- gorig
            #theResult[[i]] <- theTest
            result <- theTest
        }
        else if (output == "summaries")
            #theResult[[i]] <- list(
            result <- list(
                baseIndex=trainIndex,
                targetIndex=testIndex,
                betas=prsbetas(theTest),
                pr2=attr(theTest,"PR2")
            )
        # Remove copied GDS
        unlink(gdest,force=TRUE)
        
        # Iteratively save the result on the root workspace to be able to
        # continue later in case of crash
        #disp("\n----- Saving results up to iteration ",i," to ",saveFile,
        #    "-----\n")
        #save(theResult,file=saveFile)
        
        disp("\n----- Saving iteration ",i," result to ",saveFile,"-----\n")
        #result <- theResult[[i]]
        save(result,file=saveFile)
        
        disp("==============================================================\n")
        
        if (logging == "file") {
            sink(type="message")
            sink(type="output")
        }
        
        return(result)
        #theResult[[i]] <- result
    #}
    },rc=rc,setseed=TRUE)
    
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

.denovoQc <- function(theTrain,theTest,filters,imputeMissing,pcs,pcaMethod,
    npcs,logging) {
    qcResult <- tryCatch({
        .denovoQcWorker(theTrain,theTest,filters,imputeMissing,pcs,
            pcaMethod,npcs)
    },error=function(e) {
        .exitFromSink(logging)
        stop("Caught error during dataset QC: ",e$message,call.=FALSE)
    },interrupt=function(i) {
        .exitFromSink(logging)
        stop("Caught keyboard interruption during dataset QC!",call.=FALSE)
    },finally="")
    return(qcResult)
}

.denovoQcWorker <- function(theTrain,theTest,filters,imputeMissing,pcs,
    pcaMethod,npcs) {
    disp("\n----- Training (base) QC -----")
    theTrain <- filterGWAS(theTrain,filters=filters,
        imputeMissing=imputeMissing)
    if (pcs) {
        disp("\n----- Training (base) PCA -----\n")
        theTrain <- calcPcaCovar(theTrain,method=pcaMethod,npc=npcs)
        pcno <- ncol(pcaCovariates(theTrain))
    }
            
    # Target QC
    disp("\n----- Test (target) QC -----")
    theTest <- filterGWAS(theTest,filters=filters,
        imputeMissing=imputeMissing)
    if (pcs) {
        disp("\n----- Test (target) PCA -----\n")
        theTest <- calcPcaCovar(theTest,method=pcaMethod,npc=pcno)
    }
    return(list(train=theTrain,test=theTest))
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

.runPRSWorker <- function(theTrain,theTest,phenotype,covariates,pcs,
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
}


.prsPipelineExternal <- function(gwe,phenotype,covariates,pcs,npcs,snpSelection,
    trainSize,niter,filters,pcaMethod,imputeMissing,imputeMethod,gwaMethod,
    gwaOpts,prsiceOpts,prsWorkspace,logging,rc) {
    # Initialize the list of GWASExperiment s
    theResult <- vector("list",niter)
    
    # For output directory (workspace) format
    dig <- nchar(as.character(niter))
    pcno <- 0
    
    # Do we continue a crashed/interrupted run?
    if (length(niter) > 1)
        iters <- niter
    else
        iters <- seq_len(niter)
        
    # The RData file to save iteratively the growing result oject
    saveFile <- file.path(prsWorkspace,"external_iter_data.RData")
        
    # The actual worker
    for (i in iters) {
        pad <- paste0(rep("0",dig - nchar(as.character(i))),collapse="")
        iterWspace <- file.path(prsWorkspace,paste0(format(Sys.time(),
            "%Y%m%d%H%M%S"),"_external_",pad,i))
        dir.create(iterWspace,recursive=TRUE,showWarnings=FALSE)
        
        # If not provided by the user, PLINK and SNPTEST workspaces should live
        # within the main PRS analysis workspace
        if (gwaMethod %in% c("snptest","plink") && is.null(gwaOpts$workspace)) {
            gwaOpts$workspace <- file.path(iterWspace,
                paste0(gwaMethod,"_",.randomString()))
            dir.create(gwaOpts$workspace,recursive=TRUE,showWarnings=FALSE)
        }
        
        if (logging == "file") {
            # Print something basic before sinking
            disp("Executing external pipeline iteration ",i)
            fh <- file(file.path(prsWorkspace,paste0("iter_",pad,i,".log")),
                open="wt")
            sink(fh,type="output")
            sink(fh,type="message")
        }
        
        disp("\n==============================================================")
        disp("-----> External pipeline iteration ",i)
        disp("==============================================================\n")
        
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

        # Run basic GWAS on base to get coefficients (from one method for now)
        disp("\n----- Base GWAS -----\n")
        #base <- gwa(base,phenotype,covariates,pcs=pcs,methods="glm",rc=rc)
        rcm <- paste0('gwa(base,phenotype,covariates,pcs=pcs,',
            'methods=gwaMethod,',paste0(gwaMethod,"Opts"),'=gwaOpts,rc=rc)')
        base <- eval(parse(text=rcm))
        
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
        
        # Iteratively save the result on the root workspace to be able to
        # continue later in case of crash
        disp("\n----- Saving results up to iteration ",i," to ",saveFile,
            "-----\n")
        save(theResult,file=saveFile)
        
        disp("==============================================================\n")
        
        if (logging == "file") {
            sink(type="message")
            sink(type="output")
        }
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

.partialCleanup <- function(wspace) {
    subdirs <- dir(wspace,full.names=TRUE)
    lapply(subdirs,function(x) {
        plinks <- dir(x,pattern=".bed$|.bim$|.fam$",full.names=TRUE)
        sums <- dir(x,pattern="^base|^covar|^pheno",full.names=TRUE)
        gwas <- dir(x,pattern="^snptest_|^plink_",full.names=TRUE)
        if (length(plinks) > 0)
            unlink(plinks,recursive=TRUE,force=TRUE)
        if (length(sums) > 0)
            unlink(sums,recursive=TRUE,force=TRUE)
        if (length(gwas) > 0)
            unlink(gwas,recursive=TRUE,force=TRUE)
    })
}

.exitFromSink <- function(s) {
    if (s == "file") {
        sink(type="message")
        sink(type="output")
    }
}

.prettyLogOptions <- function() {
    # Display initialization report
    disp(strftime(Sys.time()),": Data processing started...\n")
    ############################################################################
    disp("Read counts file: ",countsName)
    disp("Conditions: ",paste(names(sampleList),collapse=", "))
    disp("Samples to include: ",paste(unlist(sampleList),collapse=", "))
    if (!is.null(excludeList) && !is.na(excludeList))
        disp("Samples to exclude: ",paste(unlist(excludeList),collapse=", "))
    else
        disp("Samples to exclude: none")
    if (!is.null(contrast))
        disp("Requested contrasts: ",paste(contrast,collapse=", "))
    else
        disp("Requested contrasts: none")
    if (!is.null(libsizeList)) {
        disp("Library sizes: ")
        for (n in names(libsizeList))
            disp("  ",paste(n,libsizeList[[n]],sep=": "))
    }
    if (!is.null(annotation) && !is.list(annotation))
        disp("Annotation: ",annotation)
    if (is.list(annotation))
        disp("Annotation: user provided GTF file")
    disp("Organism: ",org)
    disp("Reference source: ",refdb)
    disp("Count type: ",countType)
    if (countType == "utr") {
        disp("3' UTR fraction: ",utrOpts$frac)
        disp("3' UTR minimum length: ",utrOpts$minLength,"bps")
        disp("3' UTR downstream: ",utrOpts$downstream,"bps")
    }
    if (!is.null(preset))
        disp("Analysis preset: ",preset)
    disp("Transcriptional level: ",transLevel)
    if (!is.null(exonFilters)) {
        disp("Exon filters: ",paste(names(exonFilters),collapse=", "))
        for (ef in names(exonFilters)) {
            disp("  ",ef,": ")
            for (efp in names(exonFilters[[ef]])) {
                if (length(exonFilters[[ef]][[efp]])==1 && 
                    is.function(exonFilters[[ef]][[efp]]))
                    print(exonFilters[[ef]][[efp]])
                else if (length(exonFilters[[ef]][[efp]])==1)
                    disp("    ",paste(efp,exonFilters[[ef]][[efp]],sep=": "))
                else if (length(exonFilters[[ef]][[efp]])>1)
                    disp("    ",paste(efp,paste(exonFilters[[ef]][[efp]],
                        collapse=", "),sep=": "))
            }
        }
    }
    else
        disp("Exon filters: none applied")
    if (!is.null(geneFilters)) {
        disp("Gene filters: ",paste(names(geneFilters),collapse=", "))
        for (gf in names(geneFilters)) {
            disp("  ",gf,": ")
            for (gfp in names(geneFilters[[gf]])) {
                if (length(geneFilters[[gf]][[gfp]])==1 && 
                    is.function(geneFilters[[gf]][[gfp]]))
                    print(geneFilters[[gf]][[gfp]])
                else if (length(geneFilters[[gf]][[gfp]])==1)
                    disp("    ",paste(gfp,geneFilters[[gf]][[gfp]],sep=": "))
                else if (length(geneFilters[[gf]][[gfp]])>1)
                    disp("    ",paste(gfp,paste(geneFilters[[gf]][[gfp]],
                        collapse=", "),sep=": "))
            }
        }
    }
    else
        disp("Gene filters: none applied")
    disp("Filter application: ",whenApplyFilter)
    disp("Normalization algorithm: ",normalization)
    if (!is.null(normArgs)) {
        disp("Normalization arguments: ")
        for (na in names(normArgs)) {
            if (length(normArgs[[na]])==1 && is.function(normArgs[[na]])) {
                disp("  ",na,": ")
                disp(as.character(substitute(normArgs[[na]])))
            }
            else if (length(normArgs[[na]])==1)
                disp("  ",paste(na,normArgs[[na]],sep=": "))
            else if (length(normArgs[[na]])>1)
                disp("  ",paste(na,paste(normArgs[[na]],collapse=", "),
                    sep=": "))
        }
    }
    if (!any(is.na(statistics)))
        disp("Statistical algorithm(s): ",paste(statistics,collapse=", "))
    else
        disp("Statistical algorithm(s): no testing selected")
    if (!is.null(statArgs)) {
        if (!any(is.na(statistics))) {
            disp("Statistical arguments: ")
            for (sa in names(statArgs)) {
                if (length(statArgs[[sa]])==1 && is.function(statArgs[[sa]])) {
                    disp("  ",sa,": ")
                    disp(as.character(substitute(statArgs[[na]])))
                }
                else if (length(statArgs[[sa]])==1)
                    disp("  ",paste(sa,statArgs[[sa]],sep=": "))
                else if (length(statArgs[[sa]])>1)
                    disp("  ",paste(sa,paste(statArgs[[sa]],collapse=", "),
                        sep=": "))
            }
        }
        else
            disp("Statistical arguments: no testing selected")
    }
    disp("Meta-analysis method: ",metaP)
    disp("Multiple testing correction: ",adjustMethod)
    if (!is.na(pcut)) 
        disp("p-value threshold: ",pcut)
    disp("Logarithmic transformation offset: ",logOffset)
    if (!is.null(preset)) 
        disp("Analysis preset: ",preset)
    disp("Quality control plots: ",paste(qcPlots,collapse=", "))
    disp("Figure format: ",paste(figFormat,collapse=", "))
    if (!is.na(exportWhere)) 
        disp("Output directory: ",exportWhere)
    disp("Output data: ",paste(exportWhat,collapse=", "))
    disp("Output scale(s): ",paste(exportScale,collapse=", "))
    disp("Output values: ",paste(exportValues,collapse=", "))
    if ("stats" %in% exportWhat)
        disp("Output statistics: ",paste(exportStats,collapse=", "),"\n")
    ############################################################################
}
