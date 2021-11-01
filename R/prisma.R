prisma <- function(
    gwe,
    phenotype,
    covariates,
    pcs=FALSE,
    npcs=0,
    trainSize=0.8,
    niter=10,
    resolution=c("frequency","quantile"),
    step=if (resolution=="frequency") 1 else 
        c(0.1,0.2,0.3,0.4,0.5,0.75,0.8,0.9,0.95,0.99),
    minFreq=2,
    minSnps=5,
    dropSameQuantiles=TRUE,
    aggregation=c("intersection","union"),
    effectWeight=c("mean","median","weight"),
    prsSelectMethod=c("maxima","elbow"),
    prsSelectCrit=c("prs_r2","prs_pvalue","prs_aic"),
    prsSelectStat=c("mean","median","none"),
    prsSelectR2=c("adjusted","raw"),
    sigTest=c("ttest","wilcoxon","empirical"),
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
    dnOutput=c("gwaslist","summaries"),
    output=c("normal","full"),
    continue=FALSE,
    useDenovoWorkspace=NULL,
    runId=NULL,
    evalOnSplit=c("original","new"),
    evalWith=c("vanilla","prscice"),
    rc=NULL
) {
    # To be attached to output object
    FUN_CALL <- deparse(sys.call())
    
    prsSelectMethod <- prsSelectMethod[1]
    prsSelectCrit <- prsSelectCrit[1]
    prsSelectStat <- prsSelectStat[1]
    prsSelectR2 <- prsSelectR2[1]
    sigTest <- sigTest[1]
    dnOutput <- dnOutput[1]
    evalOnSplit <- evalOnSplit[1]
    .checkTextArgs("PRS selection method (prsSelectMethod)",prsSelectMethod,
        c("maxima","elbow"),multiarg=FALSE)
    .checkTextArgs("PRS selection criterion (prsSelectCrit)",prsSelectCrit,
        c("prs_r2","prs_pvalue","prs_aic"),multiarg=FALSE)
    .checkTextArgs("PRS selection statistics (prsSelectStat)",prsSelectStat,
        c("mean","median","none"),multiarg=FALSE)
    .checkTextArgs("PRS R2 for selection (prsSelectStat)",prsSelectR2,
        c("adjusted","raw"),multiarg=FALSE)
    .checkTextArgs("Significance test for report/plots (sigTest)",sigTest,
        c("ttest","wilcoxon","empirical"),multiarg=FALSE)
    .checkTextArgs("Evaluation on initial or new split (evalOnSplit)",
        evalOnSplit,c("original","new"),multiarg=FALSE)
    .checkTextArgs("De novo PRS extraction output (dnOutput)",dnOutput,
        c("gwaslist","summaries"),multiarg=FALSE)
    .checkTextArgs("PRISMA pipeline output (output)",output,c("normal","full"),
        multiarg=FALSE)
    .checkTextArgs("GWA methods to test (gwaMethods)",gwaMethods,
        c("glm","rrblup","statgen","snptest","plink"),multiarg=TRUE)
    
    glmOpts <- .checkGwaArgs(glmOpts,"glm")
    rrblupOpts <- .checkGwaArgs(rrblupOpts,"rrblup")
    #statgenOpts <- .checkGwaArgs(statgenOpts,"statgen")
    snptesOpts <- .checkGwaArgs(snptestOpts,"snptest")
    plinkOpts <- .checkGwaArgs(plinkOpts,"plink")
    lassosumOpts <- .checkPrsArgs(lassosumOpts,"lassosum")
    prsiceOpts <- .checkPrsArgs(prsiceOpts,"prsice")
    # Rest arguments are checked downstream
    
    # For the time being, vanilla evaluation with new splits is not implemented!
    if (evalWith == "vanilla" && evalOnSplit == "new")
        stop("Evaluation with R internal functions with a new dataset split ",
            "is not yet implemented!\nPlease use evalWith = \"prsice\" or ",
            "evalOnSplit = \"original\".",call.=FALSE)
    
    # If a workspace is not provided, we should pick one, but not in /tmp like
    # the other cases, as this process is more sensitive
    if (is.null(prsWorkspace)) {
        # Construct a workspace based on tools and date
        wspaceName <- paste(format(Sys.time(),"%Y%m%d%H%M%S"),
            paste0(gwaMethods,collapse="-"),aggregation,sep="_")
        prsWorkspace <- file.path(getwd(),wspaceName)
        prsWorkspace <- .validateWorkspacePath(prsWorkspace,"prisma")
    }
    
    # Display initialization report
    TB <- Sys.time()
    disp("\n",strftime(Sys.time()),": Data processing started...")
    
    # Here, display options
    callArgs <- as.list(match.call())[-1]
    .prettyLogOptions(callArgs,"prisma")
    
    ## Only one used for now, until we run linear optimization in effects
    #gwaMethods <- gwaMethods[1]
    # Use them all in a loop, and later some combination
    outList <- vector("list",length(gwaMethods))
    names(outList) <- gwaMethods
    
    # If a run id is not provided, construct one - careful, it will be the same
    # TODO: Think of something about unique run ids...
    if (is.null(runId))
        runId <- .randomString(1,11)
    
    # Run discovery and evaluation loop for each selected method
    for (m in gwaMethods) {
        disp("\n",.symbolBar("*",64))
        disp("Executing PRISMA pipeline with ",m)
        disp(.symbolBar("*",64))
        
        # First part - PRS candidate SNPs. prsPipeline takes care for the 
        # collection of the whole dnResult list, whether that run has failed or 
        # not. When continue=TRUE, prsPipeline will continue from this point and
        # the dnResult going to the next part will be intact. Also, the work is 
        # split to two directories in the workspace, baseline and selection
        disp("\n",.symbolBar("#",64))
        disp("1. Discovery of de novo PRS candidates")
        disp(.symbolBar("#",64),"\n")
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
            gwaMethods=m,
            gwaCombine=gwaCombine,
            glmOpts=glmOpts,
            rrblupOpts=rrblupOpts,
            statgenOpts=statgenOpts,
            snptestOpts=snptestOpts,
            plinkOpts=plinkOpts,
            prsMethods=prsMethods,
            lassosumOpts=lassosumOpts,
            prsiceOpts=prsiceOpts,
            prsWorkspace=file.path(prsWorkspace,m,"baseline"),
            cleanup=cleanup,
            logging=logging,
            output=dnOutput,
            continue=continue,
            runId=runId,
            rc=rc
        )
        
        # Second part collect evaluation metrics - if continue=TRUE, dnResult 
        # will be collected from the previous step and passed here.
        disp("\n",.symbolBar("#",64))
        disp("2. Calculation of PRS evaluation metrics")
        disp(.symbolBar("#",64),"\n")
        uDnS <- NULL
        if (evalOnSplit == "original")
            uDnS <- dnResult
        evalMetrics <- prsSelection(
            dnList=dnResult,
            gwe=gwe,
            phenotype=phenotype,
            covariates=covariates,
            pcs=pcs,
            npcs=npcs,
            trainSize=trainSize,
            niter=niter,
            resolution=resolution,
            step=step,
            minFreq=minFreq,
            minSnps=minSnps,
            dropSameQuantile=dropSameQuantile,
            aggregation=aggregation,
            effectWeight=effectWeight,
            filters=filters,
            pcaMethod=pcaMethod,
            imputeMissing=imputeMissing,
            imputeMethod=imputeMethod,
            gwaMethods=m,
            gwaCombine=gwaCombine,
            glmOpts=glmOpts,
            rrblupOpts=rrblupOpts,
            statgenOpts=statgenOpts,
            snptestOpts=snptestOpts,
            plinkOpts=plinkOpts,
            prsMethods=prsMethods,
            lassosumOpts=lassosumOpts,
            prsiceOpts=prsiceOpts,
            prsWorkspace=file.path(prsWorkspace,m,"selection"),
            cleanup=cleanup,
            logging=logging,
            output="summary", # If full, run prsSelection directly
            continue=continue,
            runId=runId,
            evalWith=evalWith,
            useDenovoWorkspace=uDnS,
            rc=rc
        )
        
        # Third part, suggest PRS
        disp("\n",.symbolBar("#",64))
        disp("3. Selection of best PRS candidate markers")
        disp(.symbolBar("#",64),"\n")
        candidates <- selectPrs(
            metrics=evalMetrics$metrics,
            snpSelection=evalMetrics$pgs,
            gwe=gwe,
            method=prsSelectMethod,
            crit=prsSelectCrit,
            stat=prsSelectStat,
            r2type=prsSelectR2,
            base=evalMetrics$baseline
        )
        
        # Make the plots and pass them to the report object later
        disp("\n",.symbolBar("#",64))
        disp("4. Wrap-up and graphics")
        disp(.symbolBar("#",64),"\n")
        plStat <- prsSelectStat
        if (prsSelectStat == "none")
            plStat <- "mean"
        plval <- ifelse(sigTest=="ttest","p_ttest",ifelse(sigTest=="wilcoxon",
            "p_wilcox","p_emp"))
        plots <- .plotPrsEvaluation(evalMetrics$baseline,evalMetrics$metrics,
            by=prsSelectCrit,pval=plval,stat=plStat)
        if (!is.null(evalMetrics$full))
            plots$frden <- .plotFreqDensities(evalMetrics$baseline,
                evalMetrics$full,by="prs_r2")
        
        # Also PRS plots for each candidate
        plots$prsScatter <- plots$prsHist <- 
            vector("list",length(candidates$others)+1)
        names(plots$prsScatter) <- names(plots$prsHist) <- 
            c(as.character(nrow(candidates$main)),names(candidates$others))
        trait <- phenotypes(gwe)[,phenotype]
        # The first
        popPrs <- PRS(gwe,candidates$main,prsiceOpts$score)
        ppl <- .plotPrsTrait(popPrs,trait,phenotype)
        plots$prsScatter[[1]] <- ppl$scatter
        plots$prsHist[[1]] <- ppl$hist
        # The rest
        for (na in names(candidates$others)) {
            popPrs <- PRS(gwe,candidates$others[[na]],prsiceOpts$score)
            ppl <- .plotPrsTrait(popPrs,trait,phenotype)
            plots$prsScatter[[na]] <- ppl$scatter
            plots$prsHist[[na]] <- ppl$hist
        }
        
        ## Display a summary of CV metrics
        #summarizeCvMetrics(cvMetrics,nrow(candidates$main))
        
        # The final output should be an object that can be used to build a
        # report but also some kind of inspection
        # .report(dnList,evalMetrics)
        # or even better construct an environment/list with all inputs (for 
        # logging) and pass to a reporting function .report(env)
        
        mainPrs <- list(candidates$main)
        names(mainPrs) <- as.character(nrow(candidates$main))
        outList[[m]] <- list(
            candidates=c(mainPrs,candidates$others),
            iterations=NULL,
            reportData=list(
                evalMetrics=evalMetrics,
                plots=plots
            )
        )
        if (output == "full")
            outList[[m]]$iterations = dnResult
    }
    
    disp("\n",strftime(Sys.time()),": Data processing finished!\n")
    execTime <- .elap2human(TB)
    disp("Total processing time: ",execTime,"\n\n")
    
    # Record call args
    allArgs <- .getPrismaMainDefaults()
    allArgs[names(callArgs)] <- callArgs
    
    return(list(
        params=list(
            call=FUN_CALL,
            args=allArgs
        ),
        results=outList
    ))
    
    ## Extra code for readjusting final effects - not so useful
    #theMainPrs <- adjustPrsWeights(candidates$main,gwe,phenotype,covariates,
    #    pcs,rc=rc)
    #
    ## Redo CV to see if the final effects improve R2
    #newCvMetrics <- prsCrossValidate(
    #    snpSelection=theMainPrs,
    #    gwe=gwe,
    #    response=phenotype,
    #    covariates=covariates,
    #    pcs=pcs,
    #    leaveOut=cvOutSize,
    #    times=ncvs,
    #    rc=rc
    #)
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
    resolution=c("frequency","quantile"),
    step=if (resolution=="frequency") 1 else 
        c(0.1,0.2,0.3,0.4,0.5,0.75,0.8,0.9,0.95,0.99),
    minFreq=2,
    minSnps=5,
    dropSameQuantiles=TRUE,
    aggregation=c("intersection","union"),
    effectWeight=c("mean","median","weight"),
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
    evalWith=c("vanilla","prscice"),
    output=c("summary","full"),
    rc=NULL
) {
    # The rest are checked downstream
    aggregation <- aggregation[1]
    effectWeight <- effectWeight[1]
    resolution <- resolution[1]
    output <- output[1]
    evalWith <- evalWith[1]
    .checkTextArgs("SNP aggregation method (aggregation)",aggregation,
        c("intersection","union"),multiarg=FALSE)
    .checkTextArgs("Output data level (output)",output,c("summary","full"),
        multiarg=FALSE)
    .checkTextArgs("Evaluation framework (evalWith)",evalWith,
        c("vanilla","prsice"),multiarg=FALSE)
    .checkNumArgs("Minimum frequency (minFreq)",minFreq,"numeric",0,"gte")
    
    # Check if the selected GWA method is in dnList
    # Consider adding an auto function, will use .getPri to find or licomb later
    gwaMethods <- gwaMethods[1]
    if (is(dnList[[1]],"GWASExperiment"))
        tmpn <- colnames(effects(dnList[[1]]))
    else
        tmpn <- colnames(dnList[[1]]$effects)
    if (!(gwaMethods %in% tmpn))
        stop("The requested effects from ",gwaMethods," does not exist ",
            "in the input list object. Has it been performed?")
    
    # Check resolution steps
    if (resolution == "quantile") {
        if (length(step) == 1)
            stop("When resolution is \"quantile\", the step must be a ",
                "numeric vector of quantiles (0-1)!")
        if (any(unlist(lapply(step),function(x) {
            if (x<0 || x>1)
                return(TRUE)
            return(FALSE)
        })))
            stop("When resolution is \"quantile\", all the requested ",
                "quantiles must be between 0-1!")
    }
    else if (resolution == "frequency")
        .checkNumArgs("Frequency step (step)",step,"numeric",0,"gte")
    
    # Start doing the job
    if (aggregation == "intersection") aggregation <- "intersect"
    snpSummary <- aggregatePrsMarkers(dnList,mode=aggregation,qcut=0,
        assoc=gwaMethods,avgfun=effectWeight)
    # Apply minimum frequency threshold
    snpSummary <- snpSummary[snpSummary$freq>=minFreq,,drop=FALSE]
    # Define evaluation plot steps
    if (resolution == "quantiles") {
        fstep <- round(quantile(snpSummary$freq,step))
        if (dropSameQuantiles)
            # Remove quantiles with the same number of SNPs (same frequency) -
            # keep the largest quantile
            fstep <- fstep[!duplicated(fstep,fromLast=TRUE)]
    }
    else if (resolution == "frequency") {
        fstep <- sort(unique(snpSummary$freq))
        if (step > 1) {
            # We would need at least 20 data points in this case...
            if (step > round(max(fstep)/20)) {
                warning("The chosen frequency is too large and does not leave ",
                    "enough datapoints for reliable evaluation! ",
                    "Auto-adjusting...",immediate.=TRUE)
                step <- round(max(fstep)/20)
            }
            fstep <- seq(from=fstep[1],to=fstep[length(fstep)],by=step)
        }
    }
    
    # Finally, ensure that the final frequency step contains more than minSnps
    # SNPs - if yes, exclude them, otherwise crashes, besides PRS with very
    # few SNPs is not PRS
    fval <- snpSummary[snpSummary$freq>=fstep[length(fstep)],,drop=FALSE]
    while (nrow(fval) <= minSnps) {
        fstep <- fstep[-length(fstep)]
        fval <- snpSummary[snpSummary$freq>=fstep[length(fstep)],,drop=FALSE]
    }
    #if (nrow(fval) == 1)
    #    fstep <- fstep[-length(fstep)]
    
    # Here, display options if called independently 
    if (!(grepl("prisma\\(",deparse(sys.calls()[sys.nframe()-1])[1]))) {
        callArgs <- as.list(match.call())[-1]
        .prettyLogOptions(callArgs,"select")
    }
    
    # Validation - selection runs
    # If continue=TRUE, prsPipeline will detect if the run is complete in the
    # workspace and return a full exResult, otherwise it will continue.
    # In this way, one can simply collect the results of a previous runs
    # without a dedicated function for this and by simply re-running the
    # initial command with the initial arguments - the workspace must exist!
    # Some tweaking per selection threshold is required.
    metrics <- vector("list",length(fstep))
    counter <- 0
    for (n in fstep) {
        counter <- counter + 1
        snpSelection <- snpSummary[snpSummary$freq>=n,,drop=FALSE]
        
        currSpace <- file.path(prsWorkspace,n)
        currType <- ifelse(is.null(useDenovoWorkspace),"external","evaluate")
        currCont <- ifelse(.isPrismaWorkspace(currSpace,currType),TRUE,FALSE)
        
        message("\n",.symbolBar("=",64))
        message("-----> Frequency ",n," <-----")
        message(.symbolBar("=",64))
        exResult <- prsPipeline(
            gwe=gwe,
            phenotype=phenotype,
            covariates=covariates,
            pcs=pcs,
            npcs=npcs,
            trainSize=trainSize,
            niter=niter,
            snpSelection=snpSelection,
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
            prsWorkspace=currSpace,
            cleanup=cleanup,
            logging=logging,
            output="summaries",
            continue=currCont,
            useDenovoWorkspace=useDenovoWorkspace,
            runId=as.character(n),
            evalWith=evalWith,
            rc=rc
        )
        
        if (evalWith == "prsice") {
            metrics[[counter]] <- .getR2(exResult)
            nsnp <- attr(exResult[[1]],"nsnp")
            metrics[[counter]]$n_snp <- rep(nsnp,nrow(metrics[[counter]]))
        }
        else if (evalWith == "vanilla")
            metrics[[counter]] <- do.call("rbind",exResult)
    }

    # The PRS selection function never does denovo. So there are 3 cases of
    # pipeline execution and downstream calculations
    # - External validation: use the SNPs in the aggregated data frame and
    #   perform de novo training and test splits for each frequency cutoff and
    #   then use PRSice for validation (R2). Very time consuming, effects are
    #   recalculated each time, but maybe the most independent. Rarely used
    #   although quite tested. Maybe useful for PGS Catalog SNPs.
    #   Output: a list of niterx4 data frames with the three PRSice2 R2 metrics
    # TODO: Allow internal validation with this method
    # - Evaluation with PRSice: use the SNPs in the aggregated data frame and
    #   run PRSice in 'apply' mode with the given SNPs at specific frequencies.
    #   Much faster as it used the original train/test splits used in the de
    #   novo PRS construction.
    #   Output: a list of niterx4 data frames with the three PRSice2 R2 metrics
    # - Evaluation with internal regression: use the SNPs in the aggregated
    #   data frame with GLM regression in R. The fastest case which allows a
    #   more fine grained solution as it runs internally and uses the original
    #   train/test splits and returns more metrics.
    #   Output: a list of niterx10 data frames with the three PRSice2 R2 metrics
    
    # Here we collect
    # - Number of snps in fstep
    # - Mean/median of evaluation metrics in fstep
    # - Std/IQR of evaluation metrics in fstep
    # - Adjusted (sqrt(R2/log(N))) evaluation metrics
    # - p-value against the baseline (R2 only)
    baseline <- .getR2(dnList)[,"r2p"]
    if (evalWith == "prsice")
        freqMetrics <- lapply(metrics,function(x) {
            # for p-values
            b <- .subsampleBase(x$r2p,baseline)
            
            return(c(
                n_snp=x$n_snp[1],
                freq=x$freq[1],
                mean_prs_r2=mean(x$r2p,na.rm=TRUE),
                mean_full_r2=mean(x$r2m,na.rm=TRUE),
                mean_reduced_r2=mean(x$r2n,na.rm=TRUE),
                sd_prs_r2=sd(x$r2p,na.rm=TRUE),
                sd_full_r2=sd(x$r2m,na.rm=TRUE),
                sd_reduced_r2=sd(x$r2n,na.rm=TRUE),
                median_prs_r2=median(x$r2p,na.rm=TRUE),
                median_full_r2=median(x$r2m,na.rm=TRUE),
                median_reduced_r2=median(x$r2n,na.rm=TRUE),
                iqr_prs_r2=IQR(x$r2p,na.rm=TRUE),
                iqr_full_r2=IQR(x$r2m,na.rm=TRUE),
                iqr_reduced_r2=IQR(x$r2n,na.rm=TRUE),
                mean_prs_aic=NA,
                mean_full_aic=NA,
                mean_reduced_aic=NA,
                sd_prs_aic=NA,
                sd_full_aic=NA,
                sd_reduced_aic=NA,
                median_prs_aic=NA,
                median_full_aic=NA,
                median_reduced_aic=NA,
                iqr_prs_aic=NA,
                iqr_full_aic=NA,
                iqr_reduced_aic=NA,
                prs_p=NA,
                full_p=max(x$p),
                reduced_p=NA,
                p_emp=length(which(x$r2p < b))/length(b),
                p_ttest=t.test(x$r2p,b,alternative="greater")$p.value,
                p_wilcox=wilcox.test(x$r2p,b,alternative="greater")$p.value
            ))
        })
    else if (evalWith == "vanilla")
        freqMetrics <- lapply(metrics,function(x) {
            # for p-values
            b <- .subsampleBase(x$prs_r2,baseline)
            
            return(c(
                n_snp=x$n_snp[1],
                freq=x$freq[1],
                mean_prs_r2=mean(x$prs_r2,na.rm=TRUE),
                mean_full_r2=mean(x$full_r2,na.rm=TRUE),
                mean_reduced_r2=mean(x$reduced_r2,na.rm=TRUE),
                sd_prs_r2=sd(x$prs_r2,na.rm=TRUE),
                sd_full_r2=sd(x$full_r2,na.rm=TRUE),
                sd_reduced_r2=sd(x$reduced_r2,na.rm=TRUE),
                median_prs_r2=median(x$prs_r2,na.rm=TRUE),
                median_full_r2=median(x$full_r2,na.rm=TRUE),
                median_reduced_r2=median(x$reduced_r2,na.rm=TRUE),
                iqr_prs_r2=IQR(x$prs_r2,na.rm=TRUE),
                iqr_full_r2=IQR(x$full_r2,na.rm=TRUE),
                iqr_reduced_r2=IQR(x$reduced_r2,na.rm=TRUE),
                mean_prs_aic=mean(x$prs_aic,na.rm=TRUE),
                mean_full_aic=mean(x$full_aic,na.rm=TRUE),
                mean_reduced_aic=mean(x$reduced_aic,na.rm=TRUE),
                sd_prs_aic=sd(x$prs_aic,na.rm=TRUE),
                sd_full_aic=sd(x$full_aic,na.rm=TRUE),
                sd_reduced_aic=sd(x$reduced_aic,na.rm=TRUE),
                median_prs_aic=median(x$prs_aic,na.rm=TRUE),
                median_full_aic=median(x$full_aic,na.rm=TRUE),
                median_reduced_aic=median(x$reduced_aic,na.rm=TRUE),
                iqr_prs_aic=IQR(x$prs_aic,na.rm=TRUE),
                iqr_full_aic=IQR(x$full_aic,na.rm=TRUE),
                iqr_reduced_aic=IQR(x$reduced_aic,na.rm=TRUE),
                prs_p=max(x$prs_pvalue,na.rm=TRUE),
                full_p=max(x$full_pvalue,na.rm=TRUE),
                reduced_p=max(x$reduced_pvalue,na.rm=TRUE),
                p_emp=length(which(x$prs_r2 < b))/length(b),
                p_ttest=t.test(x$prs_r2,b,alternative="greater")$p.value,
                p_wilcox=wilcox.test(x$prs_r2,b,alternative="greater")$p.value
            ))
        })
    
    freqMetrics <- do.call("rbind",freqMetrics)
    rownames(freqMetrics) <- freqMetrics[,"n_snp"]
     
    # Consider changing the previous barchart to something like returned by
    # PRSice (heatmap in colorbars denoting significance) since it's hell to
    # properly create a second axis.
    #
    # Also, user should choose the level of output: metric summaries only or
    # metric summaries and metric details.
    
    if (output == "summary")
        return(list(
            pgs=snpSummary,
            baseline=baseline,
            metrics=freqMetrics
        ))
    else if (output == "full") {
        return(list(
            pgs=snpSummary,
            baseline=baseline,
            metrics=freqMetrics,
            full=metrics
        ))
    }
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
    evalWith=c("vanilla","prscice"),
    rc=NULL
) {
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
                    "This is normal if providing a runId for the first ", "time. Reading the latest one...",call.=FALSE,
                    immediate.=TRUE)
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
        # A previous denovo workspace may have been provided...
        if (!.isPreviousDenovoList(useDenovoWorkspace) 
            && !.isPreviousDenovoList(useDenovoWorkspace))
            useDenovoWorkspace <- callParams$useDenovoWorkspace
        runId <- callParams$runId
        evalWith <- callParams$evalWith
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
            if (.isPrismaWorkspace(useDenovoWorkspace,"denovo")
                || .isPreviousDenovoList(useDenovoWorkspace)) {
                fast <- TRUE
                sdirs <- dir(prsWorkspace,pattern=paste0(runId,"_evaluate_"),
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
    evalWith <- evalWith[1]
    
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
    .checkTextArgs("Logging option (logging)",logging,c("screen","file"),
        multiarg=FALSE)
    .checkTextArgs("Cleanup option (cleanup)",cleanup,
        c("none","intermediate","all"),multiarg=FALSE)
    .checkTextArgs("Output option (output)",output,c("gwaslist","summaries"),
        multiarg=FALSE)
    .checkTextArgs("Evaluation framework (evalWith)",evalWith,
        c("vanilla","prsice"),multiarg=FALSE)
    
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
            runId=runId,
            evalWith=evalWith
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
    
    # Here, display options, if called directly and not from prsSelection
    # (otherwise too much screen flooding)
    clr <- deparse(sys.calls()[sys.nframe()-1])[1]
    if (!(grepl("prisma\\(",clr) || grepl("prsSelection\\(",clr))) {
        callParams2 <- callParams
        callParams2$gwe <- gwe
        .prettyLogOptions(callParams2,"pipeline")
    }
    disp(" ")
    
    currResult <- .prsPipeline(gwe,phenotype,covariates,pcs,npcs,snpSelection,
        trainSize,niter,filters,pcaMethod,imputeMissing,imputeMethod,
        gwaMethods,gwaCombine,glmOpts,rrblupOpts,statgenOpts,snptestOpts,
        plinkOpts,prsMethods,lassosumOpts,prsiceOpts,prsWorkspace,cleanup,
        logging,output,useDenovoWorkspace,runId,evalWith,dig,rc)
    
    #if (evalWith == "vanilla")
    #    return(do.call("rbind",currResult))
    #else
    return(c(prevResult,currResult))
}

# TODO: Average effects! - Add one more column to out with averaged effect for
# each SNP... M4BU
aggregatePrsMarkers <- function(gwaList,mode=c("intersect","union"),qcut=0.9,
    assoc=c("auto","glm","rrblup","statgen","snptest","plink","lasso"),
    avgfun=c("mean","median","weight"),gwe=NULL) {
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
    if (assoc != "auto" && !(assoc %in% colnames(tmpe))) {
        warning("The requested association method (",assoc,") cannot be ",
            "found in the input object! Switching to auto...",
            immediate.=TRUE)
        assoc <- "auto"     
    }
    if (assoc == "auto") {
        pri <- .getGwaLinArgPrior()
        if (any(pri %in% colnames(tmpe)))
            assoc <- colnames(tmpe)[which(pri %in% colnames(tmpe))[1]]
    }
    
    # TODO: lassosum effects from lassosum
    
    disp("Aggregation mode: ",mode)
    preCandidates <- lapply(seq_along(gwaList),function(i,G,m,ge) {
        x <- G[[i]]
        if (ge)
            b <- prsbetas(x)
        else
            b <- x$betas
        g <- lapply(colnames(b),function(n,b) {
            s <- rownames(b)[which(b[,n] != 0)]
        },b)
        #return(rownames(b)[which(b[,"prsice"] != 0)])
        o <- Reduce(m,g)
        disp("  iteration ",i," yields ",length(o)," markers")
        return(o)
        #return(Reduce(m,g))
    },gwaList,mode,isGwaExp)
    
    
    ## Extract lassosum effects for lassosum-only variants
    #newEffs <- lapply(preCandidates,function(x,G,ge) {
    #    x <- G[[i]]
    #    if (ge)
    #        b <- prsbetas(x)[,"lassosum"]
    #    else
    #        b <- x$betas[,"lassosum"]
    #        
    #    eff <- numeric(length(freq))
    #    names(eff) <- names(freq)
    #    curr <- intersect(names(b),names(freq))
    #    eff[curr] <- b[curr]
    #    
    #    # setdiff(lassosum,prsice)
    #    # betas[lassosum,] for the setdiff
    #},gwaList,isGwaExp)
    
    disp("Calculating frequencies")
    prsCandidates <- unlist(preCandidates)
    freq <- table(prsCandidates)
    goods <- names(freq)[freq >= floor(quantile(freq,qcut))]
    
    # Average effects for all?
    disp("Averaging effects with ",avgfun)
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
    
    disp("Constructing output")
    out <- as.data.frame(freq[goods])
    names(out) <- c("snp","freq")
    out$effect <- avgEffs[goods]
    rownames(out) <- out$snp
    
    # Correct the strange thing with effects except PLINK by reversing...
    if (assoc != "plink")
        out$effect <- -avgEffs[goods]
    
    # If a GWAS object is provided, harmonize the output with PGS
    if (!is.null(gwe))
        out <- .constructOutput(rownames(out),out,gwe)
    
    disp("Done!")
    return(out[order(out$freq,decreasing=TRUE),,drop=FALSE])
}

harvestWorkspace <- function(wspace,rid,denovo=TRUE,fast=FALSE) {
    pat <- ifelse(denovo,"denovo_",ifelse(fast,"evaluate_","external_"))
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
    useDenovoWorkspace,runId,evalWith,dig,rc) {
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
            if (evalWith == "vanilla")
                theResult <- .prsPipelineEvaluate(dnList,gwe,phenotype,
                    covariates,pcs,snpSelection,niter,pcaMethod,prsWorkspace,
                    logging,runId,dig,rc)
            else
                theResult <- .prsPipelineValidate(dnList,gwe,phenotype,
                    covariates,pcs,snpSelection,niter,pcaMethod,gwaMethods[1],
                    prsWorkspace,logging,runId,dig,rc)
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
                return(c(r2m=NA,r2n=NA,r2p=NA,p=NA,emp=NA,nsnp=NA))
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

.isPrismaWorkspace <- function(d,type) {
    furtherCheck <- ifelse(!is.null(d) && is.character(d) && dir.exists(d),
        TRUE,FALSE)
    if (furtherCheck) {
        checks <- logical(3)
        # Does it contain individual directories with the _denovo_ pattern?
        pat <- ifelse(type=="denovo","_denovo_",ifelse(type=="external",
            "_external_","_evaluate_"))
        sdirs <- dir(d,pattern=pat,full.names=TRUE)
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

.prettyLogOptions <- function(callArgs,what=c("prisma","select","pipeline")) {
    what <- what[1]
    if (what == "prisma") {
        allArgs <- .getPrismaMainDefaults()
        allArgs[names(callArgs)] <- callArgs
        .prettyLogOptionsPrisma(allArgs)
    }
    else if (what == "select") {
        allArgs <- .getPrsSelectionDefaults()
        allArgs[names(callArgs)] <- callArgs
        .prettyLogOptionsSelect(allArgs)
    }
    else if (what == "pipeline") {
        allArgs <- .getPrsPipelineDefaults()
        allArgs[names(callArgs)] <- callArgs
        .prettyLogOptionsPrs(allArgs)
    }
}

.prettyLogOptionsPrisma <- function(args) {
    disp("\nInput")
    disp("--------------------------------------------------------------------")
    disp("Input object (gwe) : a GWASExperiment object with ",
        nrow(eval(args$gwe))," markers (SNPs) and ",ncol(eval(args$gwe)),
        " samples")
    
    disp("\nRegression parameters")
    disp("--------------------------------------------------------------------")
    disp("Phenotypic response (phenotype)          : ",eval(args$phenotype))
    disp("Regression covariates (covariates)       : ",
        if (is.null(eval(args$covariates))) "-" else 
            paste(eval(args$covariates),collapse=", "))
    disp("Principal Components in covariates (pcs) : ",
        ifelse(eval(args$pcs),"Yes","No"))
    disp("Number of Principal Components (npcs)    : ",
        ifelse(eval(args$pcs),ifelse(eval(args$npcs)==0,"Auto",
            eval(args$npcs)),"-"))
        
    disp("\nPreprocessing options")
    disp("--------------------------------------------------------------------")
    disp("Sample and genotype filtering options (filters)     :")
    .dispListOpts(eval(args$filters))
    disp("Principal Components calculation method (pcaMethod) : ",
        eval(args$pcaMethod))
    disp("Impute missing genotypes (imputeMissing)            : ",
        ifelse(eval(args$imputeMissing),"Yes","No"))
    disp("Missing genotype imputation method (imputeMethod)   : ",
        eval(args$imputeMethod))
    
    disp("\nGWA analysis options")
    disp("--------------------------------------------------------------------")
    disp("Genome-Wide Association tests (gwaMethods)        : ",
        paste(eval(args$gwaMethods),collapse=", "))
    disp("Association tests combination method (gwaCombine) : ",
        eval(args$gwaCombine))
    disp("Association tests options                         : ")
    for (g in eval(args$gwaMethods)) {
        disp("  ",g," : ")
        eval(parse(text=paste(".dispListOpts(eval(args$",g,"Opts),lo=2)",
            sep="")))
    }
    
    disp("\nPRS candidate extraction options")
    disp("--------------------------------------------------------------------")
    disp("PRS extraction algorithms (prsMethods) : ",
        paste(eval(args$prsMethods),collapse=", "))
    disp("PRS algorithms options (prsMethods)    : ")
    for (p in eval(args$prsMethods)) {
        disp("  ",p," : ")
        eval(parse(text=paste(".dispListOpts(eval(args$",p,"Opts),lo=2)",
            sep="")))
    }
    
    disp("\nPRS construction and evaluation options")
    disp("--------------------------------------------------------------------")
    disp("Training set size (trainSize)                               : ",
        paste0(100*eval(args$trainSize),"%"))
    disp("Number of PRS selection iterations (niter)                  : ",
        eval(args$niter))
    disp("PRS selection resolution method (resolution)                : ",
        eval(args$resolution))
    disp("PRS selection resolution step (step)                        : ",
        if (length(eval(args$step)) == 1) eval(args$step) else 
            paste(eval(args$step),collapse=", "))
    disp("Minimum SNP frequency in PRS (minFreq)                      : ",
        eval(args$minFreq))
    disp("Minimum number of SNPs in PRS (minSnps)                     : ",
        eval(args$minSnps))
    disp("Drop quantiles with same number of SNPs (dropSameQuantiles) : ",
        ifelse(eval(args$dropSameQuantiles),"Yes","No"))
    disp("PRS aggregation method (aggregation)                        : ",
        eval(args$aggregation))
    disp("SNP effect summarization method (effectWeight)              : ",
        eval(args$effectWeight))
    disp("PRS evaluation framework (evalWith)                         : ",
        eval(args$evalWith))
        
    disp("\nPRS selection options")
    disp("--------------------------------------------------------------------")
    disp("PRS selection method (prsSelectMethod)           : ",
        eval(args$prsSelectMethod))
    disp("PRS selection criterion (prsSelectCrit)          : ",
        eval(args$prsSelectCrit))
    disp("PRS criterion summarization (prsSelectStat)      : ",
        eval(args$prsSelectStat))
    disp("R^2 type if criterion includes R^2 (prsSelectR2) : ",
        eval(args$prsSelectR2))
    
    disp("\nMiscellaneous options")
    disp("--------------------------------------------------------------------")
    disp("Local workspace directory (prsWorkspace) : ",
        ifelse(is.null(eval(args$prsWorkspace)),"auto",eval(args$prsWorkspace)))
    disp("Cleanup level (cleanup)                  : ",eval(args$cleanup))
    disp("PRS pipeline output type (dnOutput)      : ",eval(args$output))
    disp("Output level (output)                    : ",eval(args$output))
    disp("Run ID (runId)                           : ",
        ifelse(is.null(eval(args$runId)),"auto",eval(args$runId)))
    disp("Available cores fraction (rc)            : ",
        ifelse(is.null(eval(args$rc)),"single core",
            paste0(100*eval(args$rc),"%")))
}

.prettyLogOptionsSelect <- function(args) {
    disp("\nInput")
    disp("--------------------------------------------------------------------")
    disp("Input object (gwe) : a GWASExperiment object with ",
        nrow(eval(args$gwe))," markers (SNPs) and ",ncol(eval(args$gwe)),
        " samples")
    
    disp("\nRegression parameters")
    disp("--------------------------------------------------------------------")
    disp("Phenotypic response (phenotype)          : ",eval(args$phenotype))
    disp("Regression covariates (covariates)       : ",
        if (is.null(eval(args$covariates))) "-" else 
            paste(eval(args$covariates),collapse=", "))
    disp("Principal Components in covariates (pcs) : ",
        ifelse(eval(args$pcs),"Yes","No"))
    disp("Number of Principal Components (npcs)    : ",
        ifelse(eval(args$pcs),ifelse(eval(args$npcs)==0,"Auto",
        eval(args$npcs)),"-"))
        
    disp("\nPreprocessing options")
    disp("--------------------------------------------------------------------")
    disp("Sample and genotype filtering options (filters)     :")
    .dispListOpts(eval(args$filters))
    disp("Principal Components calculation method (pcaMethod) : ",
        eval(args$pcaMethod))
    disp("Impute missing genotypes (imputeMissing)            : ",
        ifelse(eval(args$imputeMissing),"Yes","No"))
    disp("Missing genotype imputation method (imputeMethod)   : ",
        eval(args$imputeMethod))
    
    disp("\nGWA analysis options")
    disp("--------------------------------------------------------------------")
    disp("Genome-Wide Association tests (gwaMethods)        : ",
        paste(eval(args$gwaMethods),collapse=", "))
    disp("Association tests combination method (gwaCombine) : ",
        eval(args$gwaCombine))
    disp("Association tests options                         : ")
    for (g in eval(args$gwaMethods)) {
        disp("  ",g," : ")
        eval(parse(text=paste(".dispListOpts(eval(args$",g,"Opts),lo=2)",
            sep="")))
    }
    
    disp("\nPRS candidate extraction options")
    disp("--------------------------------------------------------------------")
    disp("PRS extraction algorithms (prsMethods) : ",
        paste(eval(args$prsMethods),collapse=", "))
    disp("PRS algorithms options (prsMethods)    : ")
    for (p in eval(args$prsMethods)) {
        disp("  ",p," : ")
        eval(parse(text=paste(".dispListOpts(eval(args$",p,"Opts),lo=2)",
            sep="")))
    }
    
    disp("\nPRS construction and evaluation options")
    disp("--------------------------------------------------------------------")
    disp("Training set size (trainSize)                               : ",
        paste0(100*eval(args$trainSize),"%"))
    disp("Number of PRS selection iterations (niter)                  : ",
        eval(args$niter))
    disp("PRS selection resolution method (resolution)                : ",
        eval(args$resolution))
    disp("PRS selection resolution step (step)                        : ",
        if (length(eval(args$step)) == 1) eval(args$step) else 
            paste(eval(args$step),collapse=", "))
    disp("Minimum SNP frequency in PRS (minFreq)                      : ",
        eval(args$minFreq))
    disp("Minimum number of SNPs in PRS (minSnps)                     : ",
        eval(args$minSnps))
    disp("Drop quantiles with same number of SNPs (dropSameQuantiles) : ",
        ifelse(eval(args$dropSameQuantiles),"Yes","No"))
    disp("PRS aggregation method (aggregation)                        : ",
        eval(args$aggregation))
    disp("SNP effect summarization method (effectWeight)              : ",
        eval(args$effectWeight))
    disp("PRS evaluation framework (evalWith)                         : ",
        eval(args$evalWith))
    
    disp("\nMiscellaneous options")
    disp("--------------------------------------------------------------------")
    disp("Local workspace directory (prsWorkspace) : ",
        ifelse(is.null(eval(args$prsWorkspace)),"auto",eval(args$prsWorkspace)))
    disp("Cleanup level (cleanup)                  : ",eval(args$cleanup))
    disp("Output type (output)                     : ",eval(args$output))
    disp("Run ID (runId)                           : ",
        ifelse(is.null(eval(args$runId)),"auto",eval(args$runId)))
    disp("Available cores fraction (rc)            : ",
        ifelse(is.null(eval(args$rc)),"single core",
            paste0(100*eval(args$rc),"%")))
}

.prettyLogOptionsPrs <- function(args) {
    disp("\nInput")
    disp("--------------------------------------------------------------------")
    disp("Input pipeline list (dnList) : ",
        if (is(eval(args$dnList[[1]]),"GWASExperiment")) 
            "a list of GWASExperiment objects" else
            "a list of PRS pipeline summaries")
    disp("Input dataset object (gwe)   : a GWASExperiment object with ",
        nrow(eval(args$gwe))," markers (SNPs) and ",
        ncol(eval(args$gwe))," samples")
    
    disp("\nRegression parameters")
    disp("--------------------------------------------------------------------")
    disp("Phenotypic response (phenotype)          : ",eval(args$phenotype))
    disp("Regression covariates (covariates)       : ",
        if (is.null(eval(args$covariates))) "-" else 
            paste(eval(args$covariates),collapse=", "))
    disp("Principal Components in covariates (pcs) : ",
        ifelse(eval(args$pcs),"Yes","No"))
    disp("Number of Principal Components (npcs)    : ",
        ifelse(eval(args$pcs),ifelse(eval(args$npcs)==0,"Auto",
            eval(args$npcs)),"-"))
        
    disp("\nPreprocessing options")
    disp("--------------------------------------------------------------------")
    disp("Sample and genotype filtering options (filters)     :")
    .dispListOpts(eval(args$filters))
    disp("Principal Components calculation method (pcaMethod) : ",
        eval(args$pcaMethod))
    disp("Impute missing genotypes (imputeMissing)            : ",
        ifelse(eval(args$imputeMissing),"Yes","No"))
    disp("Missing genotype imputation method (imputeMethod)   : ",
        eval(args$imputeMethod))
    
    disp("\nGWA analysis options")
    disp("--------------------------------------------------------------------")
    disp("Genome-Wide Association tests (gwaMethods)        : ",
        paste(eval(args$gwaMethods),collapse=", "))
    disp("Association tests combination method (gwaCombine) : ",
        eval(args$gwaCombine))
    disp("Association tests options                         : ")
    for (g in eval(args$gwaMethods)) {
        disp("  ",g," : ")
        eval(parse(text=paste(".dispListOpts(eval(args$",g,"Opts),lo=2)",
            sep="")))
    }
    
    disp("\nPRS candidate extraction options")
    disp("--------------------------------------------------------------------")
    disp("PRS extraction algorithms (prsMethods) : ",
        paste(eval(args$prsMethods),collapse=", "))
    disp("PRS algorithms options (prsMethods)    : ")
    for (p in eval(args$prsMethods)) {
        disp("  ",p," : ")
        eval(parse(text=paste(".dispListOpts(eval(args$",p,"Opts),lo=2)",
            sep="")))
    }
    
    disp("\nPRS construction and evaluation options")
    disp("--------------------------------------------------------------------")
    disp("Training set size (trainSize)              : ",
        paste0(100*eval(args$trainSize),"%"))
    disp("Number of PRS selection iterations (niter) : ",
        eval(args$niter))
    disp("PRS evaluation framework (evalWith)        : ",
        eval(args$evalWith))
    if (!is.null(eval(args$snpSelection)))
        disp("Provided candidate SNPs for PRS            : ",
            if (is.data.frame(eval(args$snpSelection))) nrow(eval(
                args$snpSelection)) else length(eval(
                args$snpSelection))," SNPs")
    
    disp("\nMiscellaneous options")
    disp("--------------------------------------------------------------------")
    disp("Local workspace directory (prsWorkspace) : ",
        ifelse(is.null(eval(args$prsWorkspace)),"auto",eval(args$prsWorkspace)))
    disp("Cleanup level (cleanup)                  : ",eval(args$cleanup))
    disp("Output type (output)                     : ",eval(args$output))
    disp("Run ID (runId)                           : ",
        ifelse(is.null(eval(args$runId)),"auto",eval(args$runId)))
    disp("Available cores fraction (rc)            : ",
        ifelse(is.null(eval(args$rc)),"single core",
            paste0(100*eval(args$rc),"%")))
}
