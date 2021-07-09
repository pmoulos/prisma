.filterWithSnpStats <- function(obj,filters,imputeMissing=TRUE,rc=NULL) {
    disp("Performing basic filtering",level="normal")
    objBas <- .filterWithSnpStatsBasic(obj,filters)
    
    disp("Performing sample IBD filtering and classic PCA")
    objIbd <- .filterWithSnpStatsIbd(objBas,filters,rc)
    
    if (imputeMissing) {
        disp("\nImputing missing values with snpStats rules and scrime kNN")
        objIbd <- .internalImputeWithSnpStats(objIbd,rc)
    }
    
    disp("\nPerforming robust PCA filtering for automatic outlier detection")
    objPca <- .filterWithSnpStatsRobustPca(objIbd,filters)
    
    # IF any samples removed due to robust pca, recalculate PCs with LD
    if (ncol(objPca) < ncol(objIbd)) {
        disp("\nReperforming LD and PCA analysis after outlier sample removal")
        objPca <- .wrapPcaWithSnpStatsLd(objPca,filters,rc)
    }

    if (prismaVerbosity() %in% c("normal","full"))
        .filterReport(objPca)
    
    return(objPca)
}

.filterWithSnpStatsBasic <- function(obj,filters) {
    # Remember, SnpMatrix is inputed transposed! So:
    # SNP metrics at row.summmary instead of col.summary
    # Sample metrics are at col.summmary instead of row.summary
    # ATTENTION: The above does not work! We must internally transpose the
    # SnpMatrix again and output the correct indices.
    x <- t(assay(obj,1))
    snpSum <- col.summary(x)
    sampleSum <- row.summary(x)
    
    # SNP filters: call rate
    filteredSnpCall <- snpSum$Call.rate < filters$snpCallRate
    # SNP filters: call rate
    filteredMaf <- snpSum$MAF < filters$maf
    # SNP filters: Hardy-Weinberg p-value
    filteredHwe <- abs(snpSum$z.HWE) >= abs(qnorm(filters$hwe/2))
    
    # Sample filters: call rate
    filteredSampleCallRate <- sampleSum$Call.rate < filters$sampleCallRate
    # Sample filters: heterozygosity
    if (is.na(filters$heteroHard)) {
        if (filters$heteroStat == "mean") {
            loc <- mean(sampleSum$Heterozygosity,na.rm=TRUE)
            sca <- sd(sampleSum$Heterozygosity,na.rm=TRUE)
        }
        else if (filters$heteroStat == "median") {
            loc <- median(sampleSum$Heterozygosity,na.rm=TRUE)
            sca <- IQR(sampleSum$Heterozygosity,na.rm=TRUE)
        }
                
        filteredHetero <- 
            sampleSum$Heterozygosity < loc - filters$heteroFac*sca | 
                loc - sampleSum$Heterozygosity > filters$heteroFac*sca
    }
    else
        filteredHetero <- sampleSum$Heterozygosity < filters$heteroHard
    # Sample filters: inbreeding coefficient
    if (!is.na(filters$inbreed)) {
        hetF <- .calcInbreedFromSnpMatrix(x,snpSum,sampleSum)
        filteredInbreed <- hetF > filters$inbreed
    } 
    
    # Report filtered entities with metadata
    filteredSnpCallInd <- which(filteredSnpCall)
    filteredMafInd <- which(filteredMaf)
    filteredHweInd <- which(filteredHwe)
    filteredSampleCallRateInd <- which(filteredSampleCallRate)
    filteredHeteroInd <- which(filteredHetero)
    filteredInbreedInd <- which(filteredInbreed)
    
    fpar <- c("snpCallRate","maf","hwe","sampleCallRate","heteroHard","inbreed")
    filtDf <- data.frame(
        parameter=fpar,
        name=c("SNP call rate","Minor Allele Frequency",
            "Hardy-Weinberg equilibrium","Sample call rate",
            "Heterozygosity","Inbreed"),
        value=c(filters$snpCallRate,filters$maf,filters$hwe,
            filters$snpCallRate,filters$heteroHard,filters$inbreed),
        type=c(rep("SNP",3),rep("Sample",3)),
        filtered=c(length(filteredSnpCallInd),length(filteredMafInd),
            length(filteredHweInd),length(filteredSampleCallRateInd),
            length(filteredHeteroInd),length(filteredInbreedInd)),
        row.names=fpar
    )
    currFiltDf <- filterRecord(obj)
    if (nrow(currFiltDf) == 0)
        filterRecord(obj) <- filtDf
    else
        filterRecord(obj) <- rbind(currFiltDf,filtDf)
    
    # Finally
    snpInd <- Reduce("union",list(filteredSnpCallInd,filteredMafInd,
        filteredHweInd))
    sampleInd <- Reduce("union",list(filteredSampleCallRateInd,
        filteredHeteroInd,filteredInbreedInd))
    
    if (any(is.na(snpInd)))
        snpInd <- snpInd[-which(is.na(snpInd))]
    if (any(is.na(sampleInd)))
        sampleInd <- sampleInd[-which(is.na(sampleInd))]
    
    out <- obj
    if (length(snpInd) > 0)
        out <- out[-snpInd,,drop=FALSE]
    if (length(sampleInd) > 0)
        out <- out[,-sampleInd,drop=FALSE]
    
    return(out)
}

.filterWithSnpStatsIbd <- function(obj,filters,rc=NULL,.testing=FALSE) {
    m <- metadata(obj)
    
    # Read GDS
    disp("Reading GDS file ",m$gdsfile,"...")
    # New way 1st, fallback if fail
    gdsHandle <- tryCatch(snpgdsOpen(m$gdsfile,readonly=FALSE),
        error=function(e) { openfn.gds(m$gdsfile,readonly=FALSE) },
        finally="")
    
    gdsIds <- read.gdsn(index.gdsn(gdsHandle,"sample.id"))
    gdsIds <- sub("-1","",gdsIds)
    add.gdsn(gdsHandle,"sample.id",gdsIds,replace=TRUE)
    
    # LD pruning
    if (!is.na(filters$LD)) {
        disp("Performing LD pruning...\n")
        x <- t(assay(obj,1))
        snpsetIbd <- .ldPruningWithSnpRelate(gdsHandle,filters$LD,rownames(x),
            colnames(x),rc,.testing)
        disp("LD pruning finished! ",length(snpsetIbd)," SNPs will be used ",
            "for IBD analysis.")
    }
    else
        snpsetIbd <- rownames(obj)
    
    # IBD analysis
    ir <- NULL
    sr <- rownames(x)
    if (!is.na(filters$IBD)) {
        disp("Performing IBD calculation...")
        if (.testing)
            ibd <- snpgdsIBDMoM(gdsHandle,kinship=TRUE,sample.id=gdsIds,
                snp.id=snpsetIbd,num.thread=.coresFrac(rc),maf=0.1,
                autosome.only=FALSE)
        else
            ibd <- snpgdsIBDMoM(gdsHandle,kinship=TRUE,sample.id=gdsIds,
                snp.id=snpsetIbd,num.thread=.coresFrac(rc))
        disp("Performing IBD selection...")
        ibdCoeff <- snpgdsIBDSelection(ibd)
        
        # Are there related samples?
        ibdCoeff <- ibdCoeff[ibdCoeff$kinship>=filters$IBD,,drop=FALSE]
        # If yes, mark for removal
        if (nrow(ibdCoeff) > 0) {
            sr <- unique(c(ibdCoeff$ID1,ibdCoeff$ID2))
            ir <- match(sr,colnames(obj))
        }
        
        # Update filters
        filtDf <- data.frame(parameter="IBD",name="IBD",value=filters$IBD,
            type="Sample",filtered=length(ir),row.names="IBD")
        currFiltDf <- filterRecord(obj)
        if (nrow(currFiltDf) == 0)
            filterRecord(obj) <- filtDf
        else
            filterRecord(obj) <- rbind(currFiltDf,filtDf)
    }
           
    m <- metadata(obj)
    m$LDsnps <- snpsetIbd
    metadata(obj) <- m
    
    # Close GDS
    tryCatch(snpgdsClose(gdsHandle),error=function(e) {closefn.gds(gdsHandle)},
        finally="")
    
    if (!is.null(ir))
        return(obj[,-ir,drop=FALSE])
    else
        return(obj)
}

.filterWithSnpStatsRobustPca <- function(obj,filters) {
    if (!requireNamespace("rrcov"))
        stop("R package rrcov is required for robust PCA!")
    
    if (!filters$pcaOut || .isEmpty(filters$pcaRobust))
        return(obj)
    
    pco <- .robustPcaWithSnpStats(obj,method=filters$pcaRobust,npc=filters$nPC,
        transpose=TRUE)
    m <- metadata(obj)
    m$pcaRob <- pco
    metadata(obj) <- m
    jj <- which(!m$pcaRob$flag)
    
    # Update filters
    filtDf <- data.frame(parameter="pcaOut",name="Robust PCA",
        value=filters$pcaRobust,type="Sample",filtered=length(jj),
        row.names="pcaOut")
    currFiltDf <- filterRecord(obj)
    if (nrow(currFiltDf) == 0)
        filterRecord(obj) <- filtDf
    else
        filterRecord(obj) <- rbind(currFiltDf,filtDf)
    
    # Return possibly filtered object
    if (length(jj) > 0)
        return(obj[,-jj,drop=FALSE])
    else
        return(obj)
}

.ldPruningWithSnpRelate <- function(handle,ldCut,samples,snps,rc=NULL,
    .testing=FALSE) {
    mustClose <- FALSE
    if (!is(handle,"gds.class")) {
        mustClose <- TRUE
        handle <- snpgdsOpen(handle,readonly=FALSE)
    }
    
    if (.testing)
        snpSub <- snpgdsLDpruning(handle,ld.threshold=ldCut,maf=0.1,
            autosome.only=FALSE,sample.id=samples,snp.id=snps)
    else
        # Default is start.pos="random" but this requires a seed for reprod
        snpSub <- snpgdsLDpruning(handle,ld.threshold=ldCut,sample.id=samples,
            snp.id=snps,start.pos="first")
    #else {
    #    # Default is start.pos="random" but this requires a seed for reprod
    #    tmp <- capture.output({
    #        snpSub <- snpgdsLDpruning(handle,ld.threshold=filters$LD,
    #            sample.id=samples,snp.id=snps,start.pos="first")
    #    })
    #    disp(paste(tmp,collapse="\n"))
    #}
    
    if (mustClose)
        snpgdsClose(handle)
    
    return(unlist(snpSub,use.names=FALSE))
}

.wrapPcaWithSnpStatsLd <- function(obj,filters,rc=NULL,.testing=FALSE) {
    m <- metadata(obj)
    
    # Read GDS
    disp("Reading GDS file ",m$gdsfile,"...")
    gdsHandle <- snpgdsOpen(m$gdsfile,readonly=FALSE)
    gdsIds <- read.gdsn(index.gdsn(gdsHandle,"sample.id"))
    gdsIds <- sub("-1","",gdsIds)
    add.gdsn(gdsHandle,"sample.id",gdsIds,replace=TRUE)
    
    # LD pruning
    if (!is.na(filters$LD)) {
        disp("Performing LD pruning...\n")
        x <- t(assay(obj,1))
        snpsetIbd <- .ldPruningWithSnpRelate(gdsHandle,filters$LD,rownames(x),
            colnames(x),rc,.testing)
        disp("LD pruning finished! ",length(snpsetIbd)," SNPs will be used ",
            "for IBD analysis.")
    }
    else
        snpsetIbd <- rownames(obj)
           
    # Finally, perform a PCA with SNPRelate using IBD filtered samples
    disp("\nPerforming PCA...")
    if (!.isEmpty(filters$pcaRobust))
        pco <- .robustPcaWithSnpStats(obj,filters$pcaRobust,snps=snpsetIbd,
            npc=filters$nPC)
    else
        pco <- .pcaWithSnpRelate(gdsHandle,samples=colnames(obj),snps=snpsetIbd,
            npcs=filters$nPC,rc=rc,.testing=.testing)
    
    m <- metadata(obj)
    m$pcaCov <- pco
    metadata(obj) <- m
    
    # Close GDS
    snpgdsClose(gdsHandle)
    
    return(obj)
}

.robustPcaWithSnpStats <- function(obj,method=c("grid","hubert"),snps=NULL,
    npc=0,transpose=FALSE) {
    if (.isEmpty(snps)) {
        m <- metadata(obj)    
        if ("LDsnps" %in% names(m)) {
            snps <- m$LDsnps
            if (is.na(npc) || npc == 0)
                npc <- 32
        }
        else
            snps <- rownames(obj)
    }
    else # Sanity check for no crashing
        snps <- intersect(snps,rownames(obj))
    
    x <- assay(obj,1)
    y <- as(x,"numeric")
    if (any(is.na(y))) # Not imputed, quick kNN
        y <- .internalImputeKnn(y)
    y <- y[snps,,drop=FALSE]
    
    if (.isEmpty(npc))
        npc <- 0
    if (transpose)
        y <- t(y)
    if (method == "grid") {
        disp("  with grid search method")
        P <- PcaGrid(y,k=npc)
    }
    else if (method == "hubert") {
        disp("  with Hubert method")
        kmax <- ifelse(npc>10,npc,10)
        P <- PcaHubert(y,k=npc,kmax=kmax)
    }
    
    return(P)
}

.pcaWithSnpRelate <- function(handle,samples,snps,npcs=NULL,rc=NULL,
    .testing=FALSE) {
    if (.testing)
        pco <- snpgdsPCA(handle,maf=0.1,autosome.only=FALSE,sample.id=samples,
            snp.id=snps)
    else {
        #tmp <- capture.output({
        if (!.isEmpty(npcs))
            pco <- snpgdsPCA(handle,sample.id=samples,snp.id=snps,
                eigen.cnt=npcs,num.thread=.coresFrac(rc))
        else
            pco <- snpgdsPCA(handle,sample.id=samples,snp.id=snps,
                num.thread=.coresFrac(rc))
        #})
    }
    return(pco)
}
    
.internalImputeWithSnpStats <- function(obj,rc=NULL) {
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
        return(.internalImputeWithSnpStatsWorker(o,rc))
        disp("========================================")
    },parts,rc)
    
    disp("\nImputation finished, re-merging the output")
    return(do.call("rbind",O))
}

.internalImputeWithSnpStatsWorker <- function(obj,rc=NULL) {
    x <- assay(obj,1)
    if (!any(is.na(x)))
        return(obj)
    
    # Find samples that do not have any missing SNPs (there should be some!)
    nai <- which(is.na(x),arr.ind=TRUE)
    smIndMiss <- unique(nai[,"col"])
    
    if (length(smIndMiss) == ncol(x)) { # Then only scrime
        disp("No samples with non-missing values found... Using only scrime")
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
    
    train <- x[,-smIndMiss]
    missSnps <- sort(unique(nai[,"row"]))
    #missSnps <- rownames(train)[sort(unique(nai[,"row"]))]
    noMissSnps <- setdiff(seq_len(nrow(train)),missSnps)
    #noMissSnps <- setdiff(rownames(train),missSnps)
    nMissPrev <- length(missSnps)
    
    missing <- train[missSnps,,drop=FALSE]
    present <- train[noMissSnps,,drop=FALSE]
    posMiss <- gfeatures(obj)[missSnps,"position"]
    posPres <- gfeatures(obj)[noMissSnps,"position"]
    
    # Some verbosity
    disp(length(smIndMiss)," samples have missing genotypes in ",
        length(missSnps)," SNPs in total.")
    disp(ncol(train)," samples with ",length(noMissSnps)," SNPs with complete ",
     "presence will be used to train the impute model.")
    
    # Define the rules
    disp("Creating imputation rules")
    tmp <- capture.output({
        rules <- snp.imputation(t(present),t(missing),posPres,posMiss)
    })
    disp(paste(tmp,collapse="\n"))
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

.calcInbreedFromSnpMatrix <- function(x,snpSum=NULL,sampleSum=NULL) {
    # x is an original (i.e. not transposed SnpMatrix)
    if (is.null(snpSum))
        snpSum <- col.summary(x)
    if (is.null(sampleSum))
        sampleSum <- row.summary(x)
    
    calls <- !is.na(x)
    hetExp <- calls %*% (2*snpSum$MAF*(1-snpSum$MAF))
    hetObs <- sampleSum$Heterozygosity*ncol(x)*sampleSum$Call.rate
    
    return(1-(hetObs/hetExp))
}

.checkFilters <- function(f) {
    # Allowed and given values
    defaults <- getDefaults("filters")
    allowed <- names(defaults)
    given <- names(f)
    
    # Check if illegal filter names have been provided
    check <- given %in% allowed
    if (!any(check))
        stop("No valid filter name found!")
    if (!all(check)) {
        warning("The following filter names are invalid and will be ",
            "ignored:\n",paste(given[!check],collapse=", "))
        given <- given[check]
    }
    
    # heteroStat goes with heteroFac
    if (!is.null(f$heteroStat) && is.null(f$heteroFac)
        || is.null(f$heteroStat) && !is.null(f$heteroFac)) {
        warning("heteroStat and heteroFac must be provided together! Assuming ",
            "defaults for both...")
        f$heteroStat <- "median"
        f$heteroFac <- 3
    }
    # heteroHard and heteroStat, heteroFac are mutually exclusive
    if (!is.null(f$heteroHard) && !is.na(f$heteroHard) 
        && !is.null(f$heteroStat)) {
        warning("heteroHard and heteroStat filters are mutually exclusive! ",
            "Ignoring heteroHard...")
        f$heteroHard <- NULL
    }
    
    # Check remaining filter values
    if (!.isEmpty(f$snpCallRate))
        .checkNumArgs("snpCallRate filter",f$snpCallRate,"numeric",c(0,1),
            "both")
    if (!.isEmpty(f$sampleCallRate))
        .checkNumArgs("sampleCallRate filter",f$sampleCallRate,"numeric",c(0,1),
            "both")
    if (!.isEmpty(f$maf))
        .checkNumArgs("maf filter",f$maf,"numeric",c(0,1),"botheq")
    if (!.isEmpty(f$hwe))
        .checkNumArgs("hwe filter",f$hwe,"numeric",c(0,1),"both")
    if (!.isEmpty(f$LD))
        .checkNumArgs("LD filter",f$LD,"numeric",0,"gte")
    if (!.isEmpty(f$heteroFac))
        .checkNumArgs("heteroFac filter",f$heteroFac,"numeric",0,"gte")
    if (!.isEmpty(f$heteroHard) && !is.na(f$heteroHard))
        .checkNumArgs("heteroHard filter",f$heteroHard,"numeric",0,"gte")
    if (!.isEmpty(f$IBD))
        .checkNumArgs("IBD filter",f$IBD,"numeric",0,"gte")
    if (!.isEmpty(f$inbreed))
        .checkNumArgs("inbreed filter",f$inbreed,"numeric",0,"gte")
    if (!.isEmpty(f$heteroStat))
        .checkTextArgs("heteroStat filter",f$heteroStat,c("mean","median"),
            multiarg=FALSE)
    if (!.isEmpty(f$pcaRobust))
        .checkTextArgs("pcaRobust filter",f$pcaRobust,c("none","grid","hubert"),
            multiarg=FALSE)
    if (!.isEmpty(f$pcaOut) && !is.logical(f$pcaOut))
        stop("pcaOut filter should be TRUE or FALSE!")
    if (!.isEmpty(f$nPC))
        .checkNumArgs("nPC filter",f$nPC,"numeric",0,"gte")
    
    # Replace the defaults after value checking
    return(.setArg(defaults,f))
}

.filterReport <- function(obj) {
    filtDf <- filterRecord(obj)
    if (nrow(filtDf) > 0) {
        message("\n----- Filtering report: -----")
        message("Filtered SNPs:")
        if ("snpCallRate" %in% rownames(filtDf))
            message("  SNP call rate             : ",
                filtDf["snpCallRate","filtered"])
                
        if ("maf" %in% rownames(filtDf))
            message("  Minor Allele Frequency    : ",filtDf["maf","filtered"])
        
        if ("hwe" %in% rownames(filtDf))
            message("  Hardy-Weinberg equilibrium: ",filtDf["hwe","filtered"])
        
        message("\nFiltered samples:")
        if ("sampleCallRate" %in% rownames(filtDf))
            message("  Sample call rate          : ",
                filtDf["sampleCallRate","filtered"])
        
        if ("heteroHard" %in% rownames(filtDf))
            message("  Heterozygosity            : ",
                filtDf["heteroHard","filtered"])
        
        if ("inbreed" %in% rownames(filtDf))
            message("  Inbreed                   : ",
                filtDf["inbreed","filtered"])
        
        if ("IBD" %in% rownames(filtDf))
            message("  Identity By Descent (IBD) : ",filtDf["IBD","filtered"])
        
        if ("pcaOut" %in% rownames(filtDf))
            message("  Robust PCA                : ",
                filtDf["pcaOut","filtered"])
    }
}
