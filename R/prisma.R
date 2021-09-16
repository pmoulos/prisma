prisma <- function(
    gwe,
    phenotype,
    covariates,
    pcs=FALSE,
    npcs=0,
    trainSize=0.8,
    niter=10,
    quantiles=c(0.1,0.2,0.3,0.4,0.5,0.75,0.8,0.9,0.95,0.99),
    dropSameQuantiles=TRUE,
    aggregation=c("intersection","union"),
    evalR2=c("adjusted","full","null"),
    selectionTest=c("empirical","wilcoxon","ttest"),
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
    continue=FALSE,
    useDenovoWorkspace=NULL,
    runId=NULL,
    rc=NULL
) {
    # Arguments are checked downstream
    
    # Only one used for now, until we run linear optimization in effects
    gwaMethods <- gwaMethods[1]
    #
    
    # If a run id is not provided, construct one   
    if (is.null(runId))
        runId <- .randomString(1,11)
    
    # If a workspace is not provided, we should pick one, but not in /tmp like
    # the other cases, as this process is more sensitive
    if (is.null(prsWorkspace)) {
        # Construct a workspace based on tools and date
        wspaceName <- paste(format(Sys.time(),"%Y%m%d%H%M%S"),
            paste0(gwaMethods,collapse="-"),aggregation,sep="_")
        prsWorkspace <- file.path(getwd(),wspaceName)
        prsWorkspace <- .validateWorkspacePath(prsWorkspace,"prisma")
    }
    
    # First part - PRS candidate SNPs. prsPipeline takes care for the collection
    # of thw whole dnResult list, whether that run has failed or not. When
    # continue=TRUE, prsPipeline will continue from this point and the dnResult
    # goinf to the next part will be intact. Also, the work is split to two
    # directories in the workspace, baseline and selection
    dnResult <- prsPipeline(
        gwe=gwe,
        phenotype=phenotype,
        covariates=covariates,
        pcs=pcs,
        npcs=npcs,
        trainSize=trainSize,
        niter=niter,
        filters=filters,
        pcaMethod=pcaMethod,
        imputeMissing=imputeMissing,
        imputeMethod=imputeMethod,
        gwaMethods=gwaMethods,
        gwaCombine=gwaCombine,
        glmOpts=glmOpts,
        rrblupOpts=rrblupOpts,
        statgenOpts=statgenOpts,
        snptestOpts=snptestOpts,
        plinkOpts=plinkOpts,
        prsMethods=prsMethods,
        lassosumOpts=lassosumOpts,
        prsiceOpts=prsiceOpts,
        prsWorkspace=file.path(prsWorkspace,"baseline"),
        cleanup=cleanup,
        logging=logging,
        output=output,
        continue=continue,
        runId=runId,
        rc=rc
    )
    
    # Second part - if continue=TRUE, dnResult will be collected from the 
    # previous step and passed here.
    prsSelection(
        dnList=dnResult,
        gwe=gwe,
        phenotype=phenotype,
        covariates=covariates,
        pcs=pcs,
        npcs=npcs,
        trainSize=trainSize,
        niter=niter,
        quantiles=quantiles,
        dropSameQuantile=dropSameQuantile,
        aggregation=aggregation,
        filters=filters,
        pcaMethod=pcaMethod,
        imputeMissing=imputeMissing,
        imputeMethod=imputeMethod,
        gwaMethods=gwaMethods,
        gwaCombine=gwaCombine,
        glmOpts=glmOpts,
        rrblupOpts=rrblupOpts,
        statgenOpts=statgenOpts,
        snptestOpts=snptestOpts,
        plinkOpts=plinkOpts,
        prsMethods=prsMethods,
        lassosumOpts=lassosumOpts,
        prsiceOpts=prsiceOpts,
        prsWorkspace=file.path(prsWorkspace,"selection"),
        cleanup=cleanup,
        logging=logging,
        output=output,
        continue=continue,
        runId=runId,
        rc=rc
    )

    # WIP
}

prsSelection <- function(
    dnList,
    gwe,
    phenotype,
    covariates,
    pcs=FALSE,
    npcs=0,
    trainSize=0.8,
    niter=10,
    quantiles=c(0.1,0.2,0.3,0.4,0.5,0.75,0.8,0.9,0.95,0.99),
    dropSameQuantiles=TRUE,
    aggregation=c("intersection","union"),
    evalR2=c("adjusted","full","null"),
    selectionTest=c("empirical","wilcoxon","ttest"),
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
    continue=FALSE,
    useDenovoWorkspace=NULL,
    runId=NULL,
    rc=NULL
) {
    # The rest are checked downstream
    aggregation <- aggregation[1]
    evalR2 <- evalR2[1]
    selectionTest <- selectionTest[1]
    .checkTextArgs("SNP aggregation method (aggregation)",aggregation,
        c("intersection","union"),multiarg=FALSE)
    .checkTextArgs("R2 type for evaluation (evalR2)",evalR2,
        c("adjusted","full","null"),multiarg=FALSE)
    .checkTextArgs("PRS selection test (selectionTest)",selectionTest,
        c("empirical","wilcoxon","ttest"),multiarg=FALSE)
    
    disp("Creating aggregation marker list, in ",aggregation," mode.")
    snpList <- cmclapply(quantiles,function(x,a,D) {
        disp("  Quantile ",paste0(100*x,"%"))
        return(aggregatePrsMarkers(D,mode=a,qcut=x))
    },aggregation,dnList)
    names(snpList) <- paste0("Q",as.character(quantiles))

    # Remove quantiles with the same number of SNPs - keep the largest quantile
    if (dropSameQuantiles) {
        rows <- unlist(lapply(snpList,nrow))
        snpList <- snpList[!duplicated(rows,fromLast=TRUE)]
    }

    # Validation - selection runs
    # If continue=TRUE, prsPipeline will detect if the run is complete in the
    # workspace and return a full exResult, otherwise it will continue.
    r2Df <- list()
    for (qu in quantiles) {
        qun <- paste0("Q",as.character(qu))
        if (!is.null(snpList[[qun]])) {
            disp("============================================================")
            disp("----------> Quantile ",paste0(100*qu,"%")," <----------")
            disp("============================================================")
            exResult <- prsPipeline(
                gwe=gwe,
                phenotype=phenotype,
                covariates=covariates,
                pcs=pcs,
                npcs=npcs,
                trainSize=trainSize,
                niter=niter,
                snpSelection=rownames(snpList[[qun]]),
                filters=filters,
                pcaMethod=pcaMethod,
                imputeMissing=imputeMissing,
                imputeMethod=imputeMethod,
                gwaMethods=gwaMethods,
                gwaCombine=gwaCombine,
                glmOpts=glmOpts,
                rrblupOpts=rrblupOpts,
                statgenOpts=statgenOpts,
                snptestOpts=snptestOpts,
                plinkOpts=plinkOpts,
                prsMethods=prsMethods,
                lassosumOpts=lassosumOpts,
                prsiceOpts=prsiceOpts,
                prsWorkspace=prsWorkspace,
                cleanup=cleanup,
                logging=logging,
                output="summaries",
                continue=continue,
                useDenovoWorkspace=useDenovoWorkspace,
                runId=runId,
                rc=rc
            )
            
            r2Df[[qun]] <- .getR2(exResult)
        }
    }
    
    # A list of niter x 3 matrices with PRSice2 metrics
    #return(r2Df)
    
    # Or, since the function is called 'selection', continue with graphs?
    # A loop for how many tests we want
    .evalPrismaParts(dnList,r2Df,evalR2,selectionTest)
    # This function should probably become internal and the main prisma function
    # should return some summary statistics, plot/report object which is input
    # later to some report generation function.
}

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
    continue=FALSE,
    useDenovoWorkspace=NULL,
    runId=NULL,
    rc=NULL
) {
    #TODO: Log options in effect (or all), like in metaseqR
    #TODO: Check mutual exclusive or warn
    #      e.g. trainSize and niter is ignored when useDenovoWorkspace
    # snpSelection may be:
    # i) a vector of SNP names, in this case individual run effects are used
    # ii) a data frame output from aggregatePrsMarkers, in this case, the
    # effects present in this data frame are used
    
    # Can we run a GWA?
    .canRunGwa(gwe)
    # Variable to collect any previous results
    prevResult <- NULL
    
    if (continue) {
    # Strategy for pipeline continuation
    # 1. If continue=TRUE, read the parameters JSON file from workspace
    #    (params.json), with a check of course, if workspace does not exist, a
    #    warning is issued and pipeline starts from the beginning
    # 2. The completed iterations are determined by the presence of a .prog.json
    #    file inside each subdirectory ({ iter: i, done: true }). This also
    #    solves a parallelization crash. We read them and gather the i's as well
    #    as done is true or false. So the process of working with .prog.json is
    # 2.1. Write it with {iter: i, done: false} on iteration start
    # 2.2. Read it and set done: true on iteration finish
    # 3. Finally, run the pipeline workers with niter the gathered i's
        if (is.null(prsWorkspace))
            stop("Please provide a valid PRISMA workspace to continue an ",
                "interrupted pipeline!")
        if (is.character(prsWorkspace) && !dir.exists(prsWorkspace))
            stop("The PRISMA workspace ",prsWorkspace," does not exist! ",
                "Please start a new run.")
        
        if (is.null(runId))
            callParams <- fromJSON(file.path(prsWorkspace,"params.json"))
        else {
            pfile <- file.path(prsWorkspace,paste0(runId,"params.json"))
            if (file.exists(pfile))
                callParams <- fromJSON(pfile)
            else {
                warning("Parameters file for run ",runId," was not found! ",
                    "Reading the latest one...",call.=FALSE,immediate.=TRUE)
                callParams <- fromJSON(file.path(prsWorkspace,"params.json"))
            }
        }
        
        phenotype <- callParams$phenotype
        covariates <- callParams$covariates
        pcs <- callParams$pcs
        npcs <- callParams$npcs
        snpSelection <- callParams$snpSelection
        trainSize <- callParams$trainSize
        niter <- callParams$niter
        filters <- lapply(callParams$filters,function(x) {
            return(ifelse(is.null(x),NA,x))
        })
        pcaMethod <- callParams$pcaMethod
        imputeMissing <- callParams$imputeMissing
        imputeMethod <- callParams$imputeMethod
        gwaMethods <- callParams$gwaMethods
        gwaCombine <- callParams$gwaCombine
        glmOpts <- callParams$glmOpts
        rrblupOpts <- callParams$rrblupOpts
        statgenOpts <- callParams$statgenOpts
        snptestOpts <- callParams$snptestOpts
        plinkOpts <- callParams$plinkOpts
        prsMethods <- callParams$prsMethods
        lassosumOpts <- callParams$lassosumOpts
        prsiceOpts <- callParams$prsiceOpts
        prsWorkspace <- callParams$prsWorkspace
        cleanup <- callParams$cleanup
        logging <- callParams$logging
        output <- callParams$output
        useDenovoWorkspace <- callParams$useDenovoWorkspace
        runId <- callParams$runId
        # For output directory (workspace) format
        dig <- nchar(as.character(max(niter)))
        
        # Now, we must try to guess how many iterations are complete so as to
        # remove them from the total
        denovo <- fast <- FALSE
        if (is.null(snpSelection)) {
            denovo <- TRUE
            sdirs <- dir(prsWorkspace,pattern=paste0(runId,"_denovo_"),
                full.names=TRUE)
        }
        else {
            if (.isDenovoWorkspace(useDenovoWorkspace)) {
                fast <- TRUE
                sdirs <- dir(prsWorkspace,pattern=paste0(runId,"_validate_"),
                    full.names=TRUE)
            }
            else
                sdirs <- dir(prsWorkspace,pattern=paste0(runId,"_external_"),
                    full.names=TRUE)
        }
        idone <- unlist(lapply(seq_along(sdirs),function(i,s) {
            f <- file.path(s[i],".prog.json")
            if (!file.exists(f)) # Very unlikely as it's created 1st thing
                return(NULL)
            js <- fromJSON(f)
            if (js$done)
                return(js$iter)
        },sdirs))
        niter <- setdiff(seq_len(niter),idone)
        if (length(niter) == 0) {
            disp("All requsted PRISMA iterations seem to have been performed!",
                "Nothing to do but re-gathering the results...")
            return(harvestWorkspace(prsWorkspace,runId,denovo,fast))
        }
        else {
            disp("Harvesting previous iterations...")
            prevResult <- harvestWorkspace(prsWorkspace,runId,denovo,fast)
            disp("Resuming PRISMA iterations ",paste0(niter,collapse=", "))
        }
    }
    
    # Validate in normal run or resume, users may have altered params file
    pcaMethod <- pcaMethod[1]
    imputeMethod <- imputeMethod[1]
    gwaCombine <- gwaCombine[1]
    cleanup <- cleanup[1]
    output <- output[1]
    
    # niter can be a scalar or a vector of missing iterations
    if (length(niter) == 1)
        .checkNumArgs("Number of iterations (niter)",niter,"numeric",0,"gt")
    # trainSize can be a scalar (fraction) or the index of the training split
    if (length(trainSize) == 1)
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
    .checkTextArgs("Logging option",logging,c("screen","file"),
        multiarg=FALSE)
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
    
    # If parallel runs, log on screen will be a mess
    if (!is.null(rc) && logging == "screen") {
        warning("Screen log requested with pipeline parallelization! ",
            "Setting logging to files as messages will be mixed from ",
            "parallel runs.",call.=FALSE,immediate.=TRUE)
        logging <- "file"
    }
    # Check the format of SNP selection and also compatible output format
    if (!is.null(snpSelection)) {
        # Verify what was provided
        # i) a character vector or a data.frame
        if (!is.character(snpSelection) && !is.data.frame(snpSelection))
            stop("When provided, snpSelection must be either a character ",
                "vector (SNP names to use) or a data.frame from the ",
                "aggregatePrsMarkers function!")
        # ii) if data.frame, must be output from aggregatePrsMarkers or an
        # acceptable subset (containing snps and effects)
        if (is.data.frame(snpSelection) && !("effect" %in% names(snpSelection)))
            stop("When snpSelection is a data.frame, it must have at least ",
                "a column named effects. Please provide the output of ",
                "the aggregatePrsMarkers function directly.")
        
        # Output cannot be a GWASExperiment list
        if (output == "gwaslist") {
            warning("A list of GWASExperiment-s requested as output but a set ",
                "of markers was provided! Switching to 'summaries'...",
                call.=FALSE,immediate.=TRUE)
            output <- "summaries"
        }
    }
        
    # Check if the filters are given properly
    filters <- .checkFilters(filters)
    
    # Workspace valid?
    prsWorkspace <- .validateWorkspacePath(prsWorkspace,"prisma")
    
    # Write parameters to JSON file
    if (!continue) {
        # For output directory (workspace) format
        dig <- nchar(as.character(max(niter)))
        if (!is.null(runId)) # A run id may be provided by the user...
            runId <- as.character(runId)
        else # ...or a 16-byte unique run id for continuation if required
            runId <- .randomString(1,11)
        
        # In this case only, we may provide a list of a denovo run
        dnwsHelper <- NULL
        if (!.isPreviousDenovoList(useDenovoWorkspace))
            dnwsHelper <- useDenovoWorkspace

        callParams <- list(
            phenotype=phenotype,
            covariates=covariates,
            pcs=pcs,
            npcs=npcs,
            snpSelection=snpSelection,
            trainSize=trainSize,
            niter=niter,
            filters=filters,
            pcaMethod=pcaMethod,
            imputeMissing=imputeMissing,
            imputeMethod=imputeMethod,
            gwaMethods=gwaMethods,
            gwaCombine=gwaCombine,
            glmOpts=glmOpts,
            rrblupOpts=rrblupOpts,
            statgenOpts=statgenOpts,
            snptestOpts=snptestOpts,
            plinkOpts=plinkOpts,
            prsMethods=prsMethods,
            lassosumOpts=lassosumOpts,
            prsiceOpts=prsiceOpts,
            prsWorkspace=prsWorkspace,
            cleanup=cleanup,
            logging=logging,
            output=output,
            useDenovoWorkspace=dnwsHelper,
            runId=runId
        )
        
        # Check if there is a previous run in the same workspace and save the
        # parameters file
        newpfile <- file.path(prsWorkspace,"params.json")
        if (file.exists(newpfile)) {
            warning("A parameters file from a previous PRISMA run was ",
                "detected in the same workspace (",prsWorkspace,")!\nIt is ",
                "advised to use individual workspaces for each new run.\n",
                "I will keep the previous file with its run id attached to ",
                "its filename and write the new one in params.json.",
                call.=FALSE,immediate.=TRUE)
            oldPars <- fromJSON(newpfile)
            file.copy(from=newpfile,to=file.path(prsWorkspace,
                paste0(oldPars$runId,"_params.json")),overwrite=TRUE)
        }
        # Then, write the new one
        write_json(callParams,path=newpfile,auto_unbox=TRUE,null="null",
            pretty=TRUE)
    }
    
    # .prettyLogOptions(...)
    
    currResult <- .prsPipeline(gwe,phenotype,covariates,pcs,npcs,snpSelection,
        trainSize,niter,filters,pcaMethod,imputeMissing,imputeMethod,
        gwaMethods,gwaCombine,glmOpts,rrblupOpts,statgenOpts,snptestOpts,
        plinkOpts,prsMethods,lassosumOpts,prsiceOpts,prsWorkspace,cleanup,
        logging,output,useDenovoWorkspace,runId,dig,rc)
    
    return(c(prevResult,currResult))
}

# TODO: Average effects! - Add one more column to out with averaged effect for
# each SNP... M4BU
aggregatePrsMarkers <- function(gwaList,mode=c("intersect","union"),qcut=0.9,
    assoc=c("auto","glm","rrblup","statgen","snptest","plink","lasso"),
    avgfun=c("mean","median","weight")) {
    # Check if prsbetas non-empty everywhere
    isGwaExp <- FALSE
    if (is(gwaList[[1]],"GWASExperiment")) {
        check <- vapply(gwaList,function(x) {
            return(is.null(prsbetas(x)))
        },logical(1))
        isGwaExp <- TRUE
    }
    else
        check <- vapply(gwaList,function(x) {
            return(is.null(x$betas))
        },logical(1))
    if (any(check)) {
        w <- seq_along(check)[check]
        warning(length(w)," objects do not have an associated PRS analysis! ",
            "These are ",paste(w,collapse=", ")," and will be discarded",
            immediate.=TRUE)
        gwaList <- gwaList[!check]
    }
    
    # Argument validation
    mode <- mode[1]
    assoc <- assoc[1]
    avgfun <- avgfun[1]
    .checkTextArgs("Aggregation mode (mode)",mode,c("intersect","union"),
        multiarg=FALSE)
    .checkTextArgs("Association test to choose from (assoc)",assoc,
        c("auto","glm","rrblup","statgen","snptest","plink","lasso"),
        multiarg=FALSE)
    .checkTextArgs("Effect averaging method (avgfun)",avgfun,
        c("mean","median","weight"),multiarg=FALSE)
    .checkNumArgs("Quantile cutoff",qcut,"numeric",c(0,1),"both")
    
    # Select which effects to use (until all)
    if (isGwaExp)
        tmpe <- effects(gwaList[[1]])
    else
        tmpe <- gwaList[[1]]$effects
    if (assoc != "auto" && !(assoc %in% names(tmpe))) {
        warning("The requested association method (",assoc,") cannot be ",
            "found in the imput object! Switching to auto...",
            immediate.=TRUE)
        assoc <- "auto"     
    }
    if (assoc == "auto") {
        pri <- .getGwaLinArgPrior()
        if (any(pri %in% colnames(tmpe)))
            assoc <- colnames(tmpe)[which(pri %in% colnames(tmpe))[1]]
    }
    
    preCandidates <- lapply(gwaList,function(x,m,ge) {
        if (isGwaExp)
            b <- prsbetas(x)
        else
            b <- x$betas
        g <- lapply(colnames(b),function(n,b) {
            s <- rownames(b)[which(b[,n] != 0)]
        },b)
        return(Reduce(m,g))
    },mode,isGwaExp)
    prsCandidates <- unlist(preCandidates)
    freq <- table(prsCandidates)
    goods <- names(freq)[freq >= floor(quantile(freq,qcut))]
    
    # Average effects for all?
    
    allEffs <- do.call("cbind",lapply(gwaList,function(x) {
        if (isGwaExp)
            e <- effects(x)[,assoc]
        else
            e <- x$effects[,assoc]
        eff <- numeric(length(freq))
        names(eff) <- names(freq)
        curr <- intersect(names(e),names(freq))
        eff[curr] <- e[curr]
        return(eff)
    }))
    if (avgfun %in% c("mean","median"))
        avgEffs <- apply(allEffs,1,avgfun)
    else { # Weighting
        r1 <- .getR2(gwaList)[,"r2p"]
        #w <- sqrt(1/r1)
        r2 <- 1-r1
        w <- (sum(r2)/r2)/sum(sum(r2)/r2)
        #w <- r2/sum(r2)
        #hlp <- sign(diff(rev(r1)))
        #s <- c(sign(diff(r1)),hlp[1])
        #s[s==-1] <- 0.5
        #s[s==1] <- 1.5
        #w <- w + s
        avgEffs <- apply(allEffs,1,function(x,w) {
            return(sum(x*w))
        },w)
    }
    
    out <- as.data.frame(freq[goods])
    names(out) <- c("snp","freq")
    out$effect <- avgEffs[goods]
    rownames(out) <- out$snp
    return(out[order(out$freq,decreasing=TRUE),,drop=FALSE])
}

harvestWorkspace <- function(wspace,rid,denovo=TRUE,fast=FALSE) {
    pat <- ifelse(denovo,"denovo_",ifelse(fast,"validate_","external_"))
    sdirs <- dir(wspace,pattern=paste0(rid,"_",pat),full.names=TRUE)
    disp("Harvesting PRISMA workspace ",wspace)
    results <- lapply(sdirs,function(x) {
        disp("  Loading result from ",x)
        resFile <- file.path(x,paste0(pat,"result.RData"))
        if (file.exists(resFile)) {
            load(file.path(x,paste0(pat,"result.RData")))
            return(result)
        }
        #else
        #   disp("  Warning: File ",resFile," does not exist! Have all ",
        #       "iterations completed successfully?")
    })
    disp("Done!")
    return(results)
}

# The prsPipeline, if fed with a filtered object (e.g. according to GWAS
# p-value cutoffs and the SNPs therein) can be used for PGS evaluation
.prsPipeline <- function(gwe,phenotype,covariates,pcs,npcs,snpSelection,
    trainSize,niter,filters,pcaMethod,imputeMissing,imputeMethod,gwaMethods,
    gwaCombine,glmOpts,rrblupOpts,statgenOpts,snptestOpts,plinkOpts,prsMethods,
    lassosumOpts,prsiceOpts,prsWorkspace,cleanup,logging,output,
    useDenovoWorkspace,runId,dig,rc) {
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
    
    # Copy the gds file in the master workspace only in the case of no 
    # parallelization. Multi-opening is not possible...
    if (is.null(rc)) {
        gorig <- gdsfile(gwe)
        gdest <- file.path(prsWorkspace,basename(gorig))
        file.copy(from=gorig,to=gdest,overwrite=TRUE)
        gdsfile(gwe) <- gdest
    }
    
    # Main iteration
    if (is.null(snpSelection))
        theResult <- .prsPipelineDenovo(gwe,phenotype,covariates,pcs,npcs,
            trainSize,niter,filters,pcaMethod,imputeMissing,imputeMethod,
            gwaMethods,gwaCombine,glmOpts,rrblupOpts,statgenOpts,snptestOpts,
            plinkOpts,prsMethods,lassosumOpts,prsiceOpts,prsWorkspace,logging,
            output,runId,dig,rc)
    else {
        switch(gwaMethods[1],
            glm = { gwaOpts <- glmOpts },
            rrblup = { gwaOpts <- rrblupOpts },
            statgen = { gwaOpts <- statgenOpts },
            snptest = { gwaOpts <- snptestOpts },
            plink = { gwaOpts <- plinkOpts }
        )
        # This time the output is a PGS Catalog data frame
        if (is.null(useDenovoWorkspace))
            theResult <- .prsPipelineExternal(gwe,phenotype,covariates,pcs,npcs,
                snpSelection,trainSize,niter,filters,pcaMethod,imputeMissing,
                imputeMethod,gwaMethods[1],gwaOpts,prsiceOpts,prsWorkspace,
                logging,runId,dig,rc)
        else {
            # Has been validated upstream that is a valid denovo workspace or
            # a list from harvested workspace
            if (is.list(useDenovoWorkspace))
                dnList <- useDenovoWorkspace
            else {
                prm <- fromJSON(file.path(useDenovoWorkspace,"params.json"))
                dnList <- harvestWorkspace(useDenovoWorkspace,prm$runId)
            }
            theResult <- .prsPipelineValidate(dnList,gwe,phenotype,covariates,
                pcs,snpSelection,trainSize,niter,pcaMethod,gwaMethods[1],
                prsiceOpts,prsWorkspace,logging,runId,dig,rc)
        }
    }

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
    
    # If output is gwaslist and no parallelization, restore the original GDS 
    # file location (otherwise already restored)
    if (is.null(rc)) {
        if (output == "gwalist") {
            theResult <- lapply(theResult,function(x,g) {
                gdsfile(x) <- g
                return(x)
            },gorig)
        }
        # Delete the copy from the workspace
        unlink(gdest,force=TRUE)
    }
    
    return(theResult)
}

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
        if (is.data.frame(snpSelection))
            snpSelection <- rownames(snpSelection)
        
        # Run PRS with a safeguard for possibly filtered SNPs
        safeguard <- Reduce("intersect",list(snpSelection,rownames(base),
            rownames(target)))
        if (length(safeguard) < length(snpSelection)) {
            disp(length(snpSelection) - length(safeguard)," SNPs not found in ",
                "base data! Excluding...")
            subBase <- base[safeguard,,drop=FALSE]
            #subBase <- base
            subTarget <- target[safeguard,,drop=FALSE]
        }
        else {
            subBase <- base[snpSelection,,drop=FALSE]
            #subBase <- base
            subTarget <- target[snpSelection,,drop=FALSE]
        }
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
    snpSelection,trainSize,niter,pcaMethod,gwaMethod,prsiceOpts,prsWorkspace,
    logging,runId,dig,rc) {
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
        iterWspace <- file.path(prsWorkspace,paste0(runId,"_validate_",pad,i))
        dir.create(iterWspace,recursive=TRUE,showWarnings=FALSE)
        
        # Record progress
        .recordProgress(i,iterWspace)
        
        # Local result saving for restoring after possible crash
        saveFile <- file.path(iterWspace,"validate_result.RData")
        
        if (logging == "file") {
            # Print something basic before sinking
            ilf <- file.path(prsWorkspace,paste0(runId,"_iter_",pad,i,".log"))
            disp("Executing validation pipeline iteration ",i,", log is at ",
                ilf)
            fh <- file(ilf,open="wt")
            sink(fh,type="output")
            sink(fh,type="message")
        }
        else
            disp("\n")
        
        disp("==============================================================")
        disp("-----> Validation pipeline iteration ",i)
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
        
        #pcno <- 0
        if (pcs) {
            disp("Recalculating PCs...")
            #if (is(dnTmp,"GWASExperiment"))
            #    pcno <- ifelse(.hasPcaCovariates(dnObj),
            #        ncol(pcaCovariates(dnObj)),0)
            #else
            #    pcno <- ifelse(.hasPcaCovariates(gwe),
            #        ncol(pcaCovariates(gwe)),0)
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

harvestR2 <- function(obj) {
    return(.getR2(obj))
}

.getR2 <- function(obj) {
    if (is(obj[[1]],"GWASExperiment") || is.data.frame(obj[[1]])) {
        if (!is.null(attr(obj[[1]],"PR2"))) {
            tmp <- lapply(obj,function(x) {
                return(attr(x,"PR2"))
            })
        }
        else if (!is.null(attr(obj[[1]],"workspace"))) {
            tmp <- lapply(obj,function(x) {
                return(.readPrsiceR2(attr(x,"workspace")))
            })
        }
        else {
            tmp <- lapply(obj,function(x) {
                return(c(r2m=NA,r2n=NA,r2p=NA))
            })
        }
    }
    else if (is.list(obj[[1]])) {
        if (!is.null(obj[[1]]$pr2)) {
            tmp <- lapply(obj,function(x) {
                return(x$pr2)
            })
        }
    }
    
    return(do.call("rbind",tmp))
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

.isDenovoWorkspace <- function(d) {
    furtherCheck <- ifelse(!is.null(d) && is.character(d) && dir.exists(d),
        TRUE,FALSE)
    if (furtherCheck) {
        checks <- logical(3)
        # Does it contain individual directories with the _denovo_ pattern?
        sdirs <- dir(d,pattern="_denovo_",full.names=TRUE)
        checks[1] <- length(sdirs) > 0
        # Does it contain a params.json file?
        checks[2] <- file.exists(file.path(d,"params.json"))
        # Is a complete run?
        if (checks[2]) {
            prm <- fromJSON(file.path(d,"params.json"))
            nAsked <- prm$niter
            nDone <- length(unlist(lapply(seq_along(sdirs),function(i,s) {
                f <- file.path(s[i],".prog.json")
                if (!file.exists(f)) # Very unlikely as it's created 1st thing
                    return(NULL)
                js <- fromJSON(f)
                if (js$done)
                    return(js$iter)
            },sdirs)))
            checks[3] <- nAsked == nDone
        }
        return(all(checks))
    }
    else
        return(furtherCheck)
}

.isPreviousDenovoList <- function(x) {
    if (!is.list(x))
        return(FALSE)   
    if (is(x[[1]],"GWASExperiment") && !.isEmpty(effects(x[[1]]))
        && !.isEmpty(prsbetas(x[[1]])))
        return(TRUE)
    if (is(x[[1]],"list") && !.isEmpty(x[[1]]$effects)
        && !.isEmpty(x[[1]]$betas))
        return(TRUE)
    return(FALSE)
}

.exitFromSink <- function(s) {
    if (s == "file") {
        sink(type="message")
        sink(type="output")
    }
}

.prettyLogOptions <- function(callArgs,what=c("prisma","pipeline")) {
    what <- what[1]
    if (what == "prisma") {
        allArgs <- .getPrismaMainDefaults()
        allArgs[names(callArgs)] <- callArgs
        .prettyLogOptionsPrisma(allArgs)
    }
    else if (what == "pipeline") {
        allArgs <- .getPrsPipelineDefaults()
        allArgs[names(callArgs)] <- callArgs
        .prettyLogOptionsPrs(allArgs)
    }
    
    # For upstream
    #thisCall <- as.list(match.call())[-1]
    
    
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

.getPrismaMainDefaults <- function() {
    return(list(
        covariates=NULL,
        pcs=FALSE,
        npcs=0,
        trainSize=0.8,
        niter=10,
        quantiles=c(0.1,0.2,0.3,0.4,0.5,0.75,0.8,0.9,0.95,0.99),
        dropSameQuantiles=TRUE,
        aggregation="intersection",
        filters=getDefaults("filters"),
        pcaMethod="snprel",
        imputeMissing=FALSE,
        imputeMethod="single",
        gwaMethods="glm",
        gwaCombine="simes",
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
        continue=FALSE,
        useDenovoWorkspace=NULL,
        runId=NULL,
        rc=NULL
    ))
}

.getPrsPipelineDefaults <- function() {
    return(list(
        covariates=NULL,
        pcs=FALSE,
        npcs=0,
        snpSelection=NULL,
        trainSize=0.8,
        niter=10,
        filters=getDefaults("filters"),
        pcaMethod="snprel",
        imputeMissing=FALSE,
        imputeMethod="single",
        gwaMethods="glm",
        gwaCombine="simes",
        glmOpts=getDefaults("glm"),
        rrblupOpts=getDefaults("rrblup"),
        statgenOpts=getDefaults("statgen"),
        snptestOpts=getDefaults("snptest"),
        plinkOpts=getDefaults("plink"),
        prsMethods=c("lassosum","prsice"),
        lassosumOpts=getDefaults("lassosum"),
        prsiceOpts=getDefaults("prsice"),
        prsWorkspace=NULL,
        cleanup="intermediate",
        logging="screen",
        output="gwaslist",
        continue=FALSE,
        useDenovoWorkspace=NULL,
        runId=NULL,
        rc=NULL
    ))
}
