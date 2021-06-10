gwaGlm <- function(obj,response,covariates,pcs=TRUE,psig=0.05,rc=NULL,...) {
    .checkNumArgs("Association p-value",psig,"numeric",c(0,1),"both")
    # Construct the glm data frame, check if response and covariate names are
    # valid
    p <- phenotypes(obj)
    if (!(response %in% names(p)))
        stop("Response variable ",response," could not be found in the input ",
            "GWASExperiment phenotypes!")
    if (!all(covariates %in% names(p))) {
        nf <- covariates[!(covariates %in% names(p))]
        stop("Covariates ",paste0(nf,collapse=", ")," could not be found in ",
            "the input GWASExperiment phenotypes!")
    }
    
    # Reduce the phenotypes for testing and check GLM call integrity
    p <- p[,c(response,covariates)]
    glmArgs <- list(...)
    # We must ensure that if binomial requested, res is 0-1 and binary
    if ("family" %in% names(glmArgs) && glmArgs[["family"]] == "binomial")
        p[,response] <- .validateBinaryForBinomial(p[,response])
    
    if (pcs) { # Include robust PCs in the model
        pcov <- pcaCovariates(obj)
        p <- cbind(p,pcov)
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
    disp("Use PCs     : ",ifelse(pcs,"Yes","No"))
    assocList <- cmclapply(splits,function(n,p,r,...) {
        disp("  testing SNPs from ",n[1]," to ",n[length(n)])
        batch <- lapply(n,function(s,p,r,...) {
            disp("    testing ",s,level="full")
            dat <- cbind(p,t(as(snps[s,],"numeric")))
            tmp <- coefficients(.gwaGlmWorker(dat,r,...))
            rownames(tmp)[nrow(tmp)] <- s
            return(tmp[s,,drop=FALSE])
        },p,r,...)
        batch <- do.call("rbind",batch)
        rownames(batch) <- n
        return(batch)
    },p,response,...,rc=rc)
    
    assocResult <- do.call("rbind",assocList)
    disp("Found ",length(which(assocResult[,4]<psig))," associations ",
        "(uncorrected p-values).")
    
    return(assocResult)
}

# ... is other options passed to GLM
.gwaGlmWorker <- function(dat,res,...) {
    # We assume that dat contains only the deserved covariates
    ii <- which(colnames(dat)==res)
    # Before constructing formula, all the covariates must have proper names
    colnames(dat) <- make.names(colnames(dat))
    covs <- colnames(dat)[-ii]
    cres <- colnames(dat)[ii]
    f <- as.formula(paste(cres,paste0(covs,collapse="+"),sep="~"))
    return(summary(glm(f,data=dat,...)))
}

gwaBlup <- function(obj,response,covariates,pcs=TRUE,psig=0.05,rc=NULL) {
    # 1. Prepare the pheno/geno objects from obj for input to rrBLUP::GWAS
    # 2. Prepare the covariate names and include PCs (indepently of rrBLUP)
    # 3. Run and aggregate?
}
