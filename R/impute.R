.internalImputeWithSnpStats <- function(obj,mode=c("single","split"),rc=NULL) {
    mode <- mode[1]
    
    # Imputation per chromosome
    map <- gfeatures(obj)
    if ("chromosome" %in% names(map))
        parts <- split(obj,map$chromosome)
    else {
        splitFactor <- .splitFactorForParallel(nrow(obj),rc)
        parts <- split(obj,splitFactor)
    }
    
    disp("\nStarting imputation analysis in ",length(parts)," chunks")
    O <- lapply(names(parts),function(x,prt,rc) {
        disp("\n========== Imputing chromosome/part ",x)
        o <- prt[[x]]
        return(.internalImputeWithSnpStatsWorker(o,mode,rc))
        disp("========================================")
    },parts,rc)
    
    disp("\nImputation finished, re-merging the output")
    return(do.call("rbind",O))
}

# mode "single" for internal imputation on the same dataset, "split" for
# for splitting in train (samples with non-missing SNPs) and imputable (samples)
# with missing SNPs. If the number of samples is not big, "split" rather
# difficult to achieve
.internalImputeWithSnpStatsWorker <- function(obj,mode=c("single","split"),
    rc=NULL) {
    mode <- mode[1]
    x <- assay(obj,1)
    if (!any(is.na(x)))
        return(obj)
    
    # Find samples that do not have any missing SNPs (there should be some!)
    nai <- which(is.na(x),arr.ind=TRUE)
    smIndMiss <- unique(nai[,"col"])
    
    if (mode == "split" && length(smIndMiss) == ncol(x)) { # Then only scrime
        disp("No samples with non-missing values found... Using only scrime ",
            "for now, otherwise consider setting mode='single'.")
        y <- as(x,"numeric")
        y <- .internalImputeKnn(y)
        assay(obj,1) <- SnpMatrix(y)
        return(obj)
    }
    
    ############################################################################
    #! Because of a bug in snpStats leading to segmentation fault in imputation
    #! process, iterative imputation cannot take place... So we do it once and
    #! then impute the remaining with scrime...
    ## If found, split the dataset to derive imputation rules and interate
    #iter <- nMissCurr <- 0
    #nMissPrev <- 1
    #while (nrow(nai) > 0 || nMissCurr == nMissPrev) {
    #    iter <- iter + 1
    #    message("  imputation iteration ",iter)
    #    train <- x[,-smIndMiss]
    #    missSnps <- sort(unique(nai[,"row"]))
    #    #missSnps <- rownames(train)[sort(unique(nai[,"row"]))]
    #    noMissSnps <- setdiff(seq_len(nrow(train)),missSnps)
    #    #noMissSnps <- setdiff(rownames(train),missSnps)
    #    nMissPrev <- length(missSnps)
    #    
    #    missing <- train[missSnps,,drop=FALSE]
    #    present <- train[noMissSnps,,drop=FALSE]
    #    posMiss <- gfeatures(obj)[missSnps,"position"]
    #    posPres <- gfeatures(obj)[noMissSnps,"position"]
    #    
    #    # Define the rules
    #    rules <- snp.imputation(t(present),t(missing),posPres,posMiss)
    #    
    #    # Split the missing value matrix per sample to avoid for  
    #    # unnecessary replacements
    #    # CHECK: Potential parallelization
    #    h <- split(rownames(nai),nai[,"col"])
    #    for (i in smIndMiss) {
    #        simp <- impute.snps(rules,t(x[,i]),as.numeric=FALSE)
    #        snpmiss <- h[[as.character(i)]]
    #        x[snpmiss,i] <- simp[,snpmiss]
    #    }
    #    
    #    # Recheck missing values and continue
    #    nai <- which(is.na(x),arr.ind=TRUE)
    #    smIndMiss <- unique(nai[,"col"])
    #    nMissCurr <- length(unique(nai[,"row"]))
    #}
    ############################################################################
    
    # Define the rules
    disp("Creating imputation rules")
    if (mode == "single")
        rules <- .makeSingleImputationRules(x,gfeatures(obj),nai,smIndMiss)
    else if (mode == "split")
        rules <- .makeSplitImputationRules(x,gfeatures(obj),nai,smIndMiss)

    if (!any(can.impute(rules)))
        disp(length(rules)," imputation rules created but no SNP can be ",
            "imputed... Will skip and use kNN...")
    else {
        disp(length(rules)," imputation rules sucessfully created! With these ",
            length(which(can.impute(rules)))," genotypes can be imputed.")
        # Split the missing value matrix per sample to avoid unnecessary
        # replacements
        # CHECK: Potential parallelization
        disp("Imputing...")
        h <- split(rownames(nai),nai[,"col"])
        for (i in smIndMiss) {
            snpmiss <- h[[as.character(i)]]
            disp("  for sample ",colnames(x)[i]," ",length(snpmiss)," SNPs")
            disp("    ",paste(snpmiss,collapse="\n    "),level="full")
            simp <- impute.snps(rules,t(x[,i]),as.numeric=FALSE)
            x[snpmiss,i] <- simp[,snpmiss]
        }
    }

    # If NAs remain, scrime
    if (any(is.na(x))) {
        disp(length(which(is.na(x)))," missing values remaining... Will ",
            "use genotype kNN imputation.")
        y <- as(x,"numeric")
        y <- .internalImputeKnn(y)
        # We do not want to see the coercion message
        tmp <- capture.output({
            assay(obj,1) <- SnpMatrix(y)
        })
        return(obj)
    }
    
    assay(obj,1) <- x
    return(obj)
}

.makeSplitImputationRules <- function(x,feat,nai,smIndMiss) {
    train <- x[,-smIndMiss]
    missSnps <- sort(unique(nai[,"row"]))
    #missSnps <- rownames(train)[sort(unique(nai[,"row"]))]
    noMissSnps <- setdiff(seq_len(nrow(train)),missSnps)
    #noMissSnps <- setdiff(rownames(train),missSnps)
    nMissPrev <- length(missSnps)
    
    missing <- train[missSnps,,drop=FALSE]
    present <- train[noMissSnps,,drop=FALSE]
    posMiss <- feat[missSnps,"position"]
    posPres <- feat[noMissSnps,"position"]
    
    tmp <- capture.output({
        rules <- snp.imputation(t(present),t(missing),posPres,posMiss)
    })
    disp(paste(tmp,collapse="\n"))
    
    disp(length(smIndMiss)," samples have missing genotypes in ",
        length(missSnps)," SNPs in total.")
    disp(ncol(train)," samples with ",length(noMissSnps)," SNPs with complete ",
     "presence will be used to train the impute model.")
    
    return(rules)
}

.makeSingleImputationRules <- function(x,feat,nai,smIndMiss) {
    missSnps <- unique(nai[,"row"])
    posX <- feat[,"position"]
    
    tmp <- capture.output({
        rules <- snp.imputation(X=t(x),pos.X=posX,pos.Y=posX)
    })
    disp(paste(tmp,collapse="\n"))
    
    disp(length(smIndMiss)," samples have missing genotypes in ",
        length(missSnps)," SNPs in total.")
    
    return(rules)
}

.internalImputeKnn <- function(x) {
    if (!requireNamespace("scrime"))
        stop("R package scrime is required!")
    
    # Finalize genotypes
    x <- .toIntMat(x)
    #ximp <- knncatimputeLarge(x+1L,nn=5,
    #   verbose=prismaVerbosity() %in% c("normal","full"))
    # Still someting may go wrong, if yes, resort to mean
    ximp <- tryCatch({
        knncatimputeLarge(x+1L,nn=5,
            verbose=prismaVerbosity() %in% c("normal","full"))
    },error=function(e) {
        disp("Caught error during kNN imputation: ",e$message)
        disp("Imputing with simple average genotype...")
        ival <- round(mean(x,na.rm=TRUE))
        x[is.na(x)] <- ival
        return(x)
    })
    colnames(ximp) <- colnames(x)
    return(ximp)
}
