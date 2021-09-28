selectPrs <- function(metrics,snpSelection,gwe,method=c("maxima","elbow"),
    crit=c("prs_r2","prs_pvalue","prs_aic"),stat=c("mean","median","none"),
    r2type=c("adjusted","raw"),base=NULL) {
    method <- method[1]
    crit <- crit[1]
    stat <- stat[1]
    r2type <- r2type[1]
    .checkTextArgs("PRS selection method (method)",method,c("maxima","elbow"),
        multiarg=TRUE)
    .checkTextArgs("PRS selection criterion (crit)",crit,
        c("prs_r2","prs_pvalue","prs_aic"),multiarg=TRUE)
    .checkTextArgs("PRS selection statistic (stat)",stat,c("mean","median",
        "none"),multiarg=TRUE)
    .checkTextArgs("PRS R2 for selection (r2type)",r2type,c("adjusted","raw"),
        multiarg=TRUE)
    
    if (!requireNamespace("akmedoids"))
        stop("R package akmedoids is required!")
    
    # Is snpSelection sorted by frequency (required for elbow)
    if (is.unsorted(rev(snpSelection$freq)))
        snpSelection <- 
            snpSelection[order(snpSelection$freq,decreasing=TRUE),,drop=FALSE]
            
    # What type of metrics we have? Raw output from prsRegressionMetrics or
    # a summarized over iterations from prsSelection. If the former, mean/median
    # will not be available. If the latter, only mean/median available. It does
    # not matter, but there must be revert options.
    fromSel <- any(grepl("mean",colnames(metrics)))
    if (fromSel && stat=="none") # Silently fallback to mean
        stat <- "mean"
    else if (!fromSel && stat != "none") # Silently fallback to none
        stat <- "none"
    
    x <- metrics[,"n_snp"]
    switch(crit,
        prs_r2 = {
            switch(stat,
                mean = { sel <- "mean_prs_r2" },
                median = { sel <- "median_prs_r2" },
                none = { sel <- "prs_r2" }
            )
            y <- metrics[,sel]
        },
        prs_pvalue = {
            sel <- ifelse(fromSel,"reduced_p","prs_pvalue")
            y <- -log10(metrics[,sel])
        },
        prs_aic = {
            switch(stat,
                mean = { sel <- "mean_prs_aic" },
                median = { sel <- "median_prs_aic" },
                none = { sel <- "prs_aic" }
            )
            y <- -metrics[,sel]
        }
    )
    if (crit == "prs_r2" && r2type == "adjusted")
        y <- sqrt(y/log(x))
    
    if (method == "elbow") {
        # Now, find elbow...
        E <- elbow_point(x,y)
        # ...and get the SNPs as well as other metrics
        disp("The optimal number of markers for PRS based on ",crit," is ",
            round(E$x)," at ",crit,"=",round(E$y,5))
        # Make the selection
        #theSelection <- snpSelection[seq_len(round(E$x)),,drop=FALSE]
        snps <- rownames(snpSelection[seq_len(round(E$x)),,drop=FALSE])
        return(list(main=.constructOutput(snps,snpSelection,gwe),others=NULL))
    }
    else if (method == "maxima") {
        # This method will return a lot of points. They must be refined and
        # suggested as alternative results.
        im <- .localMaxima(y)
        
        # We essentially don't want the first (too many) and last (too few)
        # elements of im, unless of course there are no others
        if (length(im) > 2)
            im <- im[-c(1,length(im))]
        
        # If baseline R2 given and R2 is our metric, then use it to further
        # narrow down results
        if (!is.null(base) && crit == "prs_r2") {
            jm <- which(metrics[,sel] > mean(base))
            if (length(jm) > 0)
                im <- intersect(im,jm)
        }
        
        # Make the selection and construct list of SNP vectors
        fm <- metrics[im,c("n_snp","freq",sel),drop=FALSE]
        theBestIndex <- which(fm[,sel] == max(fm[,sel]))
        theBest <- im[theBestIndex]
        allSnps <- lapply(fm[,"freq"],function(f,D) {
            return(rownames(D[D[,"freq"] >= f,]))
        },snpSelection)
        
        disp("The optimal number of markers for PRS based on ",crit," is ",
            fm[theBestIndex,"n_snp"]," at ",crit,"=",
            round(fm[theBestIndex,sel],5))
        if (length(allSnps) > 1) {
            others <- setdiff(im,theBest)
            disp("Other possible options include:")
            for (o in others) {
                disp("  ",metrics[o,"n_snp"]," markers for PRS with ",crit,"=",
                    round(metrics[o,sel],5))
            }
        }
        
        main <- .constructOutput(allSnps[[theBestIndex]],snpSelection,gwe)
        others <- NULL
        if (length(allSnps) > 1)
            others <- lapply(allSnps[-theBestIndex],.constructOutput,
                snpSelection,gwe)
        
        return(list(main=main,others=others))
    }
}

prsRegressionMetrics <- function(snpSelection,gwe,response,covariates=NULL,
    pcs=FALSE,step=10,family=NULL,rc=NULL,...) {
    .canRunGwa(gwe)
    
    # Construct the glm data frame, check if response and covariate names are
    # valid
    p <- phenotypes(gwe)
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
        if (.hasPcaCovariates(gwe)) {
            pcov <- pcaCovariates(gwe)
            p <- cbind(p,pcov)
        }
        else
            warning("PC covariates requested in the model, but not calculated ",
                "PC covariates found! Ignoring...",immediate.=TRUE)
    }
    
    # We now need to calculate the PRS
    snps <- t(as(genotypes(gwe),"numeric"))
    
    # The idea is to bin (cumulatively in a forward/backward selection manner)
    # the snpSelection data frame and construct a PRS using the effects there
    indexList <- .makeStepList(step,nrow(snpSelection))
    
    disp("\nRegressing with GLM")
    disp("Trait(s)    : ",response)
    disp("Covariate(s): ",paste0(covariates,collapse=", "))
    disp("Regression  : ",family," (",fpredHelper,")")
    disp("Use PCs     : ",ifelse(pcs,"Yes","No"))
    
    # Full model   : phenoptype ~ covariates + PRS
    # Reduced model: phenoptype ~ covariates
    # Null model   : phenoptype ~ 1
    # We need to collect the following metrics:
    # - Number of SNPs used
    # - R^2 of the full model 
    # - R^2 of the reduced model
    # - R^2 difference (adjusted in PRSice) of full - reduced models
    # - p-value of the full model (F test of full against null)
    # - p-value of the reduced model (F test of reduced against null)
    # - p-value of the reduced against the full (F test)
    # - AIC of the full model
    # - AIC of the reduced model
    # - AIC of the difference (full - reduced), if < 0 then full better
    #https://stackoverflow.com/questions/17674148/one-p-value-for-glm-model
    
    # The reduced model
    disp("Creating the reduced model")
    dat <- p
    ii <- which(colnames(dat)==response)
    colnames(dat) <- make.names(colnames(dat))
    covs <- colnames(dat)[-ii]
    cres <- colnames(dat)[ii]
    if (length(covs) > 0)
        fr <- as.formula(paste(cres,paste0(covs,collapse="+"),sep="~"))
    else
        fr <- as.formula(paste(cres,1,sep="~"))
    disp("Reduced model formula is: ",level="full")
    disp(show(fr),level="full")
    redFit <- glm(fr,data=dat,family=family,...)
    redModel <- summary(redFit)
    
    # The null model
    disp("Creating the null model")
    fn <- as.formula(paste(cres,1,sep="~"))
    disp("Null model formula is: ",level="full")
    disp(show(fn),level="full")
    nullFit <- glm(fn,data=dat,family=family,...)
    nullModel <- summary(nullFit)
    
    # Stats for later use
    redR2 <- 1 - redModel$deviance/redModel$null.deviance
    redP <- anova(nullFit,redFit,test="F")[["Pr(>F)"]][2]
    
    disp("Testing full models")
    metricsList <- cmclapply(indexList,function(i,p,r,f,sdf,snps,rf,nf,...) {
        snpset <- sdf[i,,drop=FALSE]
        n <- rownames(snpset)
        freq <- snpset$freq[nrow(snpset)]
        
        disp("  testing PRS with SNPs from ",n[1]," to ",n[length(n)])
        thePrs <- .prs(snps[,n],sdf[n,"effect"])
        dat <- cbind(p,thePrs)
        colnames(dat)[ncol(dat)] <- "PRS"
        
        ii <- which(colnames(dat)==r)
        colnames(dat) <- make.names(colnames(dat))
        covs <- colnames(dat)[-ii]
        cres <- colnames(dat)[ii]
        fm <- as.formula(paste(cres,paste0(covs,collapse="+"),sep="~"))
        disp("Full model formula is: ",level="full")
        disp(show(fm),level="full")
        #fullFit <- glm(fm,data=dat,family=f,...)
        fullFit <- glm(fm,data=dat,family=f)
        fullModel <- summary(fullFit)
        
        # Number of SNPs
        nsnp <- length(n)
        
        # - R^2 of the full model
        fullR2 <- 1 - fullModel$deviance/fullModel$null.deviance
        
        # - p-value of the full model (F test of full against null)
        tmp <- anova(nf,fullFit,test="F")
        fullP <- tmp[["Pr(>F)"]][2]
        
        # - p-value of the reduced against the full (F test)
        tmp <- anova(rf,fullFit,test="F")
        diffP <- tmp[["Pr(>F)"]][2]
        
        # - AIC of the full model
        fullAIC <- fullModel$aic
        
        return(c(
            n_snp=nsnp,
            freq=freq,
            full_r2=fullR2,
            full_pvalue=fullP,
            prs_pvalue=diffP,
            full_aic=fullAIC
        ))
    },p,response,family,snpSelection,snps,redFit,nullFit,rc=rc)
    
    # Construct final metrics matrix
    metrics <- as.data.frame(do.call("rbind",metricsList))
    metrics$reduced_r2 <- rep(redR2,nrow(metrics))
    metrics$prs_r2 <- metrics$full_r2 - metrics$reduced_r2
    metrics$reduced_pvalue <- rep(redP,nrow(metrics))
    metrics$reduced_aic <- rep(redModel$aic,nrow(metrics))
    metrics$prs_aic <- metrics$full_aic - metrics$reduced_aic
    
    # data frame with metric
    return(metrics)
}

.constructOutput <- function(snps,snpSelection,gwe) {
    # Construct (the output
    obj <- gwe[snps,,drop=FALSE]
    df <- as.data.frame(gfeatures(obj))
    # Basic info
    out <- df[,c("chromosome","position","snp.name","allele.1","allele.2")]
    # Attach effects
    out$effect_weight <- snpSelection[snps,"effect"]
    out$OR <- exp(out$effect_weight)
    out$locus_name <- rep(NA,nrow(out))
    
    # Final alignment with the external API fetch outcomes from PGS catalog
    names(out)[c(3,4,5)] <- c("variant_id","risk_allele","reference_allele")
    out <- out[,c("chromosome","position","variant_id","risk_allele",
        "reference_allele","locus_name","effect_weight","OR")]
    gb <- genome(obj)
    if (is.null(gb))
        gb <- "nr"
    out$asm <- rep(gb,nrow(out))
    
    # Attach also the frequencies - revisit this later
    out$freq <- snpSelection[snps,"freq"]
    
    return(out)
}

#https://stackoverflow.com/questions/6836409/finding-local-maxima-and-minima
.localMaxima <- function(x) {
    y <- diff(c(-.Machine$integer.max, x)) > 0L
    y <- cumsum(rle(y)$lengths)
    y <- y[seq.int(1L, length(y), 2L)]
    if (x[[1]] == x[[2]])
        y <- y[-1]
    return(y)
}

.makeStepList <- function(step,N) {
    mo <- N%%step
    n <- (N-mo)/step
    indexList <- lapply(seq_len(n),function(i,s) {
        seq_len(i*s)
    },step)
    
    last <- NULL
    if (mo != 0)
        last <- seq_len(N)
    
    if (!is.null(last))
        return(c(indexList,list(last)))
    return(indexList)
}

#~ .evalGlmWorker <- function(dat,res,red,fam,...) {
#~     ii <- which(colnames(dat)==res)
#~     jj <- which(colnames(dat)==red)
#~     colnames(dat) <- make.names(colnames(dat))
#~     covs <- colnames(dat)[-ii]
#~     covsRed <- colnames(dat)[-c(ii,jj)]
#~     cres <- colnames(dat)[ii]
#~     ff <- as.formula(paste(cres,paste0(covs,collapse="+"),sep="~"))
#~     fr <- as.formula(paste(cres,paste0(covsRed,collapse="+"),sep="~"))
#~     fn <- as.formula(paste(cres,1,sep="~"))
#~     disp("Full model formula is: ",level="full")
#~     disp(show(ff),level="full")
#~     disp("Reduced model formula is: ",level="full")
#~     disp(show(fr),level="full")
#~     disp("Null model formula is: ",level="full")
#~     disp(show(fn),level="full")
#~     return(list(
#~         full=glm(ff,data=dat,family=fam,...),
#~         reduced=glm(fr,data=dat,family=fam,...),
#~         null=glm(fn,data=dat,family=fam,...)
#~     ))
#~ }
