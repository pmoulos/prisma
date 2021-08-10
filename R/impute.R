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
    O <- lapply(names(parts),function(x,prt,m,rc) {
        disp("\n========== Imputing chromosome/part ",x)
        o <- prt[[x]]
        return(.internalImputeWithSnpStatsWorker(o,m,rc))
    },parts,mode,rc)
    
    disp("\nImputation finished, re-merging the output")
    
    # Until we find what is going on with metadata rbinding...
    m <- metadata(obj)
    O <- do.call("rbind",O)
    metadata(O) <- m
    
    return(O)
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
        y <- .internalImputeKnn(.integerifyImputed(y))
        assay(obj,1) <- SnpMatrix(y)
        return(obj)
    }
    
    iter <- nMissCurr <- 0
    nMissPrev <- 1
    while (nrow(nai) > 0 && nMissCurr != nMissPrev) {
        iter <- iter + 1
        nMissPrev <- length(sort(unique(nai[,"row"])))
        
        message("-----> Imputation iteration ",iter)
        
        disp("Creating imputation rules")
        if (mode == "single")
            rules <- .makeSingleImputationRules(x,gfeatures(obj),nai,smIndMiss)
        else if (mode == "split")
            rules <- .makeSplitImputationRules(x,gfeatures(obj),nai,smIndMiss)
        
        # Split the missing value matrix per sample to avoid for  
        # unnecessary replacements
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
        
        # Recheck missing values and continue
        nai <- which(is.na(x),arr.ind=TRUE)
        smIndMiss <- unique(nai[,"col"])
        nMissCurr <- length(unique(nai[,"row"]))
        
        tmp <- capture.output({
            x <- SnpMatrix(.integerifyImputed(x,na=FALSE))
        })
    }
    ############################################################################

    # If NAs remain, scrime
    if (any(is.na(x))) {
        disp(length(which(is.na(x)))," missing values remaining... Will ",
            "use genotype kNN imputation.")
        # Genotypes must be in range 1-3, as(x,"integer") collapses to vector
        #! Yikes, this converts a 1.8 imputed genotype to 1...
        #y <- as(x,"numeric") + 1 
        #!
        y <- .internalImputeKnn(.integerifyImputed(x))
        
        if (any(is.na(y)))
            disp(length(which(is.na(y)))," missing values remaining... This ",
                "should not happen, have you removed SNPs missing in all ",
                "samples? See filterGWAS()")
        # We do not want to see the coercion message
        tmp <- capture.output({
            #assay(obj,1) <- SnpMatrix(.integerifyImputed(x,na=FALSE))
            assay(obj,1) <- SnpMatrix(y)
        })
        return(obj)
    }
    
    assay(obj,1) <- x
    return(obj)
}

#~ .internalImputeWithSnpStatsWorker <- function(obj,mode=c("single","split"),
#~     rc=NULL) {
#~     mode <- mode[1]
#~     x <- assay(obj,1)
#~     if (!any(is.na(x)))
#~         return(obj)
    
#~     # Find samples that do not have any missing SNPs (there should be some!)
#~     nai <- which(is.na(x),arr.ind=TRUE)
#~     smIndMiss <- unique(nai[,"col"])
    
#~     if (mode == "split" && length(smIndMiss) == ncol(x)) { # Then only scrime
#~         disp("No samples with non-missing values found... Using only scrime ",
#~             "for now, otherwise consider setting mode='single'.")
#~         y <- as(x,"numeric")
#~         y <- .internalImputeKnn(y)
#~         assay(obj,1) <- SnpMatrix(y)
#~         return(obj)
#~     }
    
#~     ############################################################################
#~     #! Because of a bug in snpStats leading to segmentation fault in imputation
#~     #! process, iterative imputation cannot take place... So we do it once and
#~     #! then impute the remaining with scrime...
#~     ## If found, split the dataset to derive imputation rules and interate
#~     #iter <- nMissCurr <- 0
#~     #nMissPrev <- 1
#~     #while (nrow(nai) > 0 || nMissCurr == nMissPrev) {
#~     #    iter <- iter + 1
#~     #    message("  imputation iteration ",iter)
#~     #    train <- x[,-smIndMiss]
#~     #    missSnps <- sort(unique(nai[,"row"]))
#~     #    #missSnps <- rownames(train)[sort(unique(nai[,"row"]))]
#~     #    noMissSnps <- setdiff(seq_len(nrow(train)),missSnps)
#~     #    #noMissSnps <- setdiff(rownames(train),missSnps)
#~     #    nMissPrev <- length(missSnps)
#~     #    
#~     #    missing <- train[missSnps,,drop=FALSE]
#~     #    present <- train[noMissSnps,,drop=FALSE]
#~     #    posMiss <- gfeatures(obj)[missSnps,"position"]
#~     #    posPres <- gfeatures(obj)[noMissSnps,"position"]
#~     #    
#~     #    # Define the rules
#~     #    rules <- snp.imputation(t(present),t(missing),posPres,posMiss)
#~     #    
#~     #    # Split the missing value matrix per sample to avoid for  
#~     #    # unnecessary replacements
#~     #    # CHECK: Potential parallelization
#~     #    h <- split(rownames(nai),nai[,"col"])
#~     #    for (i in smIndMiss) {
#~     #        simp <- impute.snps(rules,t(x[,i]),as.numeric=FALSE)
#~     #        snpmiss <- h[[as.character(i)]]
#~     #        x[snpmiss,i] <- simp[,snpmiss]
#~     #    }
#~     #    
#~     #    # Recheck missing values and continue
#~     #    nai <- which(is.na(x),arr.ind=TRUE)
#~     #    smIndMiss <- unique(nai[,"col"])
#~     #    nMissCurr <- length(unique(nai[,"row"]))
#~     #}
#~     ############################################################################

#~     # Define the rules
#~     disp("Creating imputation rules")
#~     if (mode == "single")
#~         rules <- .makeSingleImputationRules(x,gfeatures(obj),nai,smIndMiss)
#~     else if (mode == "split")
#~         rules <- .makeSplitImputationRules(x,gfeatures(obj),nai,smIndMiss)

#~     if (!any(can.impute(rules)))
#~         disp(length(rules)," imputation rules created but no SNP can be ",
#~             "imputed... Will skip and use kNN...")
#~     else {
#~         disp(length(rules)," imputation rules sucessfully created! With these ",
#~             length(which(can.impute(rules)))," genotypes can be imputed.")
#~         # Split the missing value matrix per sample to avoid unnecessary
#~         # replacements
#~         # CHECK: Potential parallelization
#~         disp("Imputing...")
#~         h <- split(rownames(nai),nai[,"col"])
#~         for (i in smIndMiss) {
#~             snpmiss <- h[[as.character(i)]]
#~             disp("  for sample ",colnames(x)[i]," ",length(snpmiss)," SNPs")
#~             disp("    ",paste(snpmiss,collapse="\n    "),level="full")
#~             simp <- impute.snps(rules,t(x[,i]),as.numeric=FALSE)
#~             x[snpmiss,i] <- simp[,snpmiss]
#~         }
#~     }

#~     # If NAs remain, scrime
#~     if (any(is.na(x))) {
#~         disp(length(which(is.na(x)))," missing values remaining... Will ",
#~             "use genotype kNN imputation.")
#~         # Genotypes must be in range 1-3, as(x,"integer") collapses to vector
#~         #! Yikes, this converts a 1.8 imputed genotype to 1...
#~         #y <- as(x,"numeric") + 1 
#~         #!
#~         y <- round(as(x,"numeric")) + 1 
#~         mode(y) <- "integer"
#~         y <- .internalImputeKnn(y)
        
#~         if (any(is.na(y)))
#~             disp(length(which(is.na(y)))," missing values remaining... This ",
#~                 "should not happen, please file a bug.")
#~         # We do not want to see the coercion message
#~         tmp <- capture.output({
#~             assay(obj,1) <- SnpMatrix(y)
#~         })
#~         return(obj)
#~     }
    
#~     assay(obj,1) <- x
#~     return(obj)
#~ }

.makeSplitImputationRules <- function(x,feat,nai,smIndMiss) {
    train <- x[,-smIndMiss]
    missSnps <- sort(unique(nai[,"row"]))
    #missSnps <- rownames(train)[sort(unique(nai[,"row"]))]
    noMissSnps <- setdiff(seq_len(nrow(train)),missSnps)
    #noMissSnps <- setdiff(rownames(train),missSnps)
    
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
    
    ximp <- tryCatch({
        ..knncatimputeLarge(x,nn=5,
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

.integerifyImputed <- function(x,na=TRUE) {
    y <- round(as(x,"numeric")) + 1 
    mode(y) <- "integer"
    if (any(is.na(y)) && !na)
        y[is.na(y)] <- 0L
    return(y)
}

################ Borrowed from scrime to deal with monomorphisms ###############

..knncatimputeLarge <- function(data,mat.na=NULL,fac=NULL,fac.na=NULL,nn=3,
    distance=c("smc","cohen","snp1norm","pcc"),n.num=100,use.weights=TRUE,
    verbose=FALSE) {
    if (is.null(mat.na)) {
        rs <- rowSums(is.na(data))
        ids.na <- which(rs > 0)
        if (length(ids.na) == 0)
            stop("There are no missing values in data.")
        rn <- rownames(data)
        mat.na <- data[ids.na, , drop = FALSE]
        data <- data[-ids.na, , drop = FALSE]
        if (!is.null(fac)) {
            fac.na <- fac[ids.na]
            fac <- fac[-ids.na]
        }
    }
    else {
        rs <- rowSums(is.na(mat.na))
        if (any(rs == 0))
            stop("At least one of the rows of mat.na does not contain missing values.")
        rs <- rowSums(is.na(data))
        if (any(rs > 0))
            stop("At least one of the rows of data contains missing values.")
        ids.na <- NULL
    }
    if (nn < 1)
        stop("nn must be at least 1.")
    if (nn > nrow(data))
        stop("nn must be smaller than or equal to the number of rows of data.")
    if ((is.null(fac) & !is.null(fac.na)) | (!is.null(fac) &
        is.null(fac.na)))
        stop("Either both or none of fac and fac.na has to be specified.")
    n.cat <- ..checkX1X2(data, mat.na)
    #heck4Monomorphism(data)
    #check4Monomorphism(mat.na)
    if (is.null(fac)) {
        fac <- rep(1, nrow(data))
        fac.na <- rep(1, nrow(mat.na))
        verbose <- FALSE
    }
    if (length(fac) != nrow(data))
        stop("The length of fac must be equal to the number of columns of data.")
    if (length(fac.na) != nrow(mat.na))
        stop("The length of fac.na must be equal to the number of columns of mat.na.")
    vec.split <- sort(unique(fac.na))
    tmp.split <- unique(fac)
    if (any(!vec.split %in% tmp.split))
        stop("At least one of the values in fac.na is not in fac.")
    distance <- match.arg(distance)
    distance <- paste(distance, "2Mats", sep = "")
    for (i in 1:length(vec.split)) {
        if (verbose)
            cat("Now considering factor ", vec.split[i], ".",
                sep = "")
        tmp.mat <- data[fac == vec.split[i], , drop = FALSE]
        tmp.matna <- mat.na[fac.na == vec.split[i], , drop = FALSE]
        out <- ..replaceNAs(tmp.mat, tmp.matna, nn = nn, distance = distance,
            n.num = n.num, use.weights = use.weights, n.cat = n.cat)
        mat.na[fac.na == vec.split[i], ] <- out
        if (verbose)
            cat(" Done.\n", sep = "")
    }
    if (is.null(ids.na))
        return(mat.na)
    mat <- matrix(0, length(rs), ncol(mat.na))
    mat[ids.na, ] <- mat.na
    mat[-ids.na, ] <- data
    rownames(mat) <- rn
    return(mat)
}

..replaceNAs <- function(mat,mat.na,nn=3,distance="smc2Mats",n.num=100,
    use.weights=TRUE,n.cat=NULL) {
    n.row <- nrow(mat.na)
    FUN <- match.fun(distance)
    sets <- c(seq(1, n.row, n.num), n.row + 1)
    sets <- unique(sets)
    for (i in 1:(length(sets) - 1)) {
        consider <- sets[i]:(sets[i + 1] - 1)
        tmp <- mat.na[consider, , drop = FALSE]
        mat.dist <- FUN(mat, tmp, n.cat = n.cat)
        colS <- colSums(mat.dist < 10^-8)
        ids1 <- which(colS > 0)
        ids2 <- which(colS == 0)
        for (j in ids1) {
            tmp.ids <- which(mat.dist[, j] < 10^-8)[1]
            tmp[j, ] <- mat[tmp.ids, ]
        }
        for (j in ids2) {
            tmp.dist <- mat.dist[, j]
            names(tmp.dist) <- 1:length(tmp.dist)
            tmp.dist <- sort(tmp.dist)[1:nn]
            tmp.ids <- as.numeric(names(tmp.dist))
            tmp.mat <- mat[tmp.ids, , drop = FALSE]
            tmp.ids2 <- which(is.na(tmp[j, ]))
            for (k in tmp.ids2) tmp[j, k] <- if (use.weights)
                ..weightMode(tmp.dist^-1, tmp.mat[, k])
            else ..modeDist(tmp.mat[, k])
        }
        mat.na[consider, ] <- tmp
    }
    return(mat.na)
}

..checkX1X2 <- function (x1,x2,impute=TRUE) {
    namex2 <- ifelse(impute,"mat.na","newdata")
    if (!is.matrix(x1) || !is.matrix(x2))
        stop("Both data and ",namex2," must be matrix objects.",call.=FALSE)
    if (ncol(x1) != ncol(x2))
        stop("data and ",namex2," must have the same number of columns.",
            call. = FALSE)
    m1 <- max(x1,na.rm=TRUE)
    m2 <- max(x2,na.rm=TRUE)
    if (m1!=m2)
        stop("data and ",namex2," must consist of the same numbers of ", 
            "categories.",call.=FALSE)
    ..checkCatMat(x1,m1)
    ..checkCatMat(x2,m2,matname=namex2)
    return(m1)
}

..checkCatMat <- function (data,nCat,matname="data") {
    rangeData <- range(data,na.rm=TRUE)
    if (rangeData[1] != 1 || nCat != round(rangeData[2]))
        stop(matname," must consist of integers between 1 and ",round(nCat),".",
            call.=FALSE)
    if (any(!data[!is.na(data)] %in% seq_len(nCat)))
        stop(matname," must consist of integers between 1 and ",nCat, ".",
            call.=FALSE)
    if (any(!(seq_len(nCat)) %in% data))
        stop("Some of the values between 1 and ",nCat," are not in ",matname,
            ".", call.=FALSE)
}

..weightMode <- function(w,v,num=TRUE,samp=FALSE) {
    tab <- by(w, v, sum, na.rm = TRUE)
    if (!samp)
        out <- names(tab)[tab == max(tab)]
    else out <- if (length(tab) == 1)
        names(tab)
    else sample(names(tab), 1, prob = tab)
    if (length(out) > 1)
        out <- sample(out, 1)
    if (num)
        out <- as.numeric(out)
    return(out)
}

..modeDist <- function (x,num=TRUE) {
    tab <- table(x)
    out <- names(tab)[tab == max(tab)]
    if (length(out) > 1)
        out <- sample(out, 1)
    if (num)
        out <- as.numeric(out)
    return(out)
}
