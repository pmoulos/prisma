prsRegress <- function(snpSelection,gwe,response,covariates=NULL,pcs=FALSE,
    step=10,family=NULL,rc=NULL,...) {
    .canRunGwa(gwe)
    
    # Construct the glm data frame, check if response and covariate names are
    # valid
    p <- phenotypes(gwe)
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
    
    
    
    # Then regress with it
    
    disp("\nRegressing with GLM")
    disp("Trait(s)    : ",response)
    disp("Covariate(s): ",paste0(covariates,collapse=", "))
    disp("Regression  : ",family," (",fpredHelper,")")
    disp("Use PCs     : ",ifelse(pcs,"Yes","No"))
    
    metricsList <- cmclapply(indexList,function(i,p,r,f,sdf,snps,...) {
        snpset <- sdf[i,,drop=FALSE]
        if (!is.null(rownames(snpset)))
            n <- rownames(snpset)
        else
            n <- seq_along(nrow(snpset))
        
        disp("  testing PRS with SNPs from ",n[1]," to ",n[length(n)])
        thePrs <- .prs(snps[,n],sdf[n,"effect"])
        dat <- cbind(p,thePrs)
        colnames(dat)[ncol(dat)] <- "PRS"
        
        tmpFull <- .gwaGlmWorker("glm",dat,r,f)
        tmpNull <- .gwaGlmWorker("glm",p,r,f)
        
        r2Full <- 1 - tmpFull$deviance/tmpFull$null.deviance
        r2Null <- 1 - tmpNull$deviance/tmpNull$null.deviance
        aicFull <- tmpFull$aic
        aicNull <- tmpNull$aic
        
        # Calculate p-value
#https://stats.stackexchange.com/questions/129958/glm-in-r-which-pvalue-represents-the-goodness-of-fit-of-entire-model
        

    },p,response,family,sdf,snps,...,rc=rc)
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
    
    return(c(indexList,last))
}





