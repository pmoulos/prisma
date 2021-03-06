#~ # No other way, SNPRelate is a strange package...
#~ add.gdsn <- utils::getFromNamespace("add.gdsn","SNPRelate")
#~ closefn.gds <- utils::getFromNamespace("closefn.gds","SNPRelate")
#~ index.gdsn <- utils::getFromNamespace("index.gdsn","SNPRelate")
#~ openfn.gds <- utils::getFromNamespace("openfn.gds","SNPRelate")
#~ read.gdsn <- utils::getFromNamespace("read.gdsn","SNPRelate")

.filterWithSnpStats <- function(obj,filters,imputeMissing=TRUE,
    imputeMode=c("single","split"),withPlink=TRUE,rc=NULL) {
    imputeMode <- imputeMode[1]
    
    disp("\nPerforming basic filtering",level="normal")
    if (withPlink)
        objBas <- .filterWithPlinkBasic(obj,filters)
    else
        objBas <- .filterWithSnpStatsBasic(obj,filters)
        
    
    #if (!is.na(filters$LD) || !is.na(filters$IBD)) {
    if (!is.na(filters$IBD)) {
        disp("Performing sample IBD filtering and LD calculation")
        objIbd <- .filterWithSnpStatsIbd(objBas,filters,rc)
    }
    else
        objIbd <- objBas
    
    if (imputeMissing && any(is.na(genotypes(obj)))) {
        disp("\nImputing missing values with snpStats rules and scrime kNN")
        objIbd <- .internalImputeWithSnpStats(objIbd,imputeMode,rc=rc)
    }
    
    if (filters$pcaOut) {
        disp("\nPerforming robust PCA for automatic outlier detection")
        objPca <- .filterWithSnpStatsRobustPca(objIbd,filters)
    }
    else
        objPca <- objIbd
    
    # IF any samples removed due to robust pca, recalculate PCs with LD
    if (ncol(objPca) < ncol(objIbd)) {
        disp("\nReperforming LD and PCA analysis after outlier sample removal")
        objPca <- .wrapPcaWithSnpStatsLd(objPca,filters,rc)
    }

    if (prismaVerbosity() %in% c("normal","full"))
        .filterReport(objPca)
    
    return(objPca)
}

.filterWithPlinkBasic <- function(obj,filters) {
    # The process is to write filtered PLINK files and re-read them to
    # to temporary minimal GWASExperiment object just for the rownames and
    # colnames... Too costly to track all else in the initial object and 
    # reattach everything to a new one.
    
    # If no filters, nothing to do
    if (.emptyFilters(filters))
        return(obj)
    
    disp("  Exporting temporary PLINK files for the calculations")
    plinkBase <- tempfile()
    writePlink(obj,outBase=plinkBase)
    
    # Basic commands
    outBase <- tempfile()
    plink <- .getToolPath("plink")
    baseCommand <- paste(
       paste0(plink," \\"),
       paste0("  --bfile ",plinkBase," \\"),
       paste0("  --out ",outBase),
       sep="\n"
    )
    
    # Prepare for heterozygosity and inbreed filtering - latter not yet
    if (!.isEmpty(filters$inbreed)) {
        warning("Inbreed coefficient filter is not yet supported when ",
            "filtering with PLINK. Ignoring...",immediate.=TRUE)
        filters$inbreed <- NA
    }
    
    goodHet <- NULL
    if (!.isEmpty(filters$heteroHard) || !.isEmpty(filters$heteroStat)) {
        hetCommand <- paste0(baseCommand," --het")
        message("\nExecuting:\n",hetCommand)
        out <- tryCatch({
            system(hetCommand,ignore.stdout=TRUE,ignore.stderr=TRUE)
            if (prismaVerbosity() == "full") {
                logfile <- paste0(outBase,".log")
                log <- .formatSystemOutputForDisp(readLines(logfile))
                disp("\nPLINK output is:\n")
                disp(paste(log,collapse="\n"),"\n")
            }
            0L
        },error=function(e) {
            message("Caught error: ",e$message)
            return(1L)
        },finally="")
        
        if (out != 1L) {
            hetData <- read.table(paste0(outBase,".het"),header=TRUE,
                check.names=FALSE)
            rownames(hetData) <- hetData$FID
            heterozygosity <- 1 - hetData[,3]/hetData[,5]
            names(heterozygosity) <- rownames(hetData)
            
            # Filter based on what is requested
            if (.isEmpty(filters$heteroHard)) {
                if (filters$heteroStat == "mean") {
                    loc <- mean(heterozygosity,na.rm=TRUE)
                    sca <- sd(heterozygosity,na.rm=TRUE)
                }
                else if (filters$heteroStat == "median") {
                    loc <- median(heterozygosity,na.rm=TRUE)
                    sca <- IQR(heterozygosity,na.rm=TRUE)
                }
                        
                remainHetero <- 
                    (heterozygosity > loc - filters$heteroFac*sca) & 
                        (loc - heterozygosity < filters$heteroFac*sca)
            }
            else
                remainHetero <- heterozygosity < filters$heteroHard
            goodHet <- names(which(remainHetero))   
        }
        else {
            warning("There was an error in calculating heterozygosity! ",
                "Ignoring...",immediate.=TRUE)
            filters$heteroStat <- filters$heteroFac <- filters$heteroHard <- NA
        }
    }
    
    # Basic filtering
    command <- paste0(baseCommand," \\\n","  --make-bed")
    # SNP filters: call rate
    if (!.isEmpty(filters$snpCallRate)) {
        toadd <- paste0("  --geno ",
            format(1 - filters$snpCallRate,scientific=FALSE))
        command <- paste0(command," \\\n",toadd)
    }
    # SNP filters: MAF
    if (!.isEmpty(filters$maf)) {
        toadd <- paste0("  --maf ",format(filters$maf,scientific=FALSE))
        command <- paste0(command," \\\n",toadd)
    }
    # SNP filters: Hardy-Weinberg p-value
    if (!.isEmpty(filters$hwe)) {
        toadd <- paste0("  --hwe ",format(filters$hwe,scientific=FALSE))
        command <- paste0(command," \\\n",toadd)
    }
    # Sample filters: call rate
    if (!.isEmpty(filters$sampleCallRate)) {
        toadd <- paste0("  --mind ",
            format(1 - filters$sampleCallRate,scientific=FALSE))
        command <- paste0(command," \\\n",toadd)
    }
    
    # Now execute to get the outcome
    message("\nExecuting:\n",command)
    out <- tryCatch({
        system(command,ignore.stdout=TRUE,ignore.stderr=TRUE)
        if (prismaVerbosity() == "full") {
            logfile <- paste0(outBase,".log")
            log <- .formatSystemOutputForDisp(readLines(logfile))
            disp("\nPLINK output is:\n")
            disp(paste(log,collapse="\n"),"\n")
        }
        0L
    },error=function(e) {
        message("Caught error: ",e$message)
        return(1L)
    },finally="")
    
    if (out != 1L) {
        # Read-in the resulting object and get its row- and col-names
        tmpIn <- list(
            bed=paste0(outBase,".bed"),
            bim=paste0(outBase,".bim"),
            fam=paste0(outBase,".fam")
        )
        tmpObj <- importGWAS(tmpIn,writeGds=FALSE,gdsOverwrite=FALSE)
        snpsKeep <- rownames(tmpObj)
        samsKeep <- colnames(tmpObj)
        # Unify with hetero
        samsKeep <- union(samsKeep,goodHet)
    }
    else
        stop("Filtering with PLINK has failed! Please use another filtering ",
            "backend...")
    
    # If all ok
    disp(length(snpsKeep)," SNPs and ",length(samsKeep)," samples remaining")
    return(obj[snpsKeep,samsKeep])
}

.filterWithSnpStatsBasic <- function(obj,filters) {
    # Remember, SnpMatrix is inputed transposed! So:
    # SNP metrics at row.summmary instead of col.summary
    # Sample metrics are at col.summmary instead of row.summary
    # ATTENTION: The above does not work! We must internally transpose the
    # SnpMatrix again and output the correct indices.
    
    disp("  Calculating SNP summaries for filtering")
    x <- t(assay(obj,1))
    snpSum <- col.summary(x)
    sampleSum <- row.summary(x)
    
    # SNP filters init
    disp("  SNP filters")
    filteredSnpCall <- filteredMaf <- filteredHwe <- logical(nrow(obj))
    
    # SNP filters: call rate
    if (!.isEmpty(filters$snpCallRate)) {
        disp("    Call rate")
        filteredSnpCall <- snpSum$Call.rate < filters$snpCallRate
    }
    # SNP filters: MAF
    if (!.isEmpty(filters$maf)) {
        disp("    Minor Allele Frequency")
        filteredMaf <- snpSum$MAF < filters$maf
    }
    # SNP filters: Hardy-Weinberg p-value
    if (!.isEmpty(filters$hwe)) {
        disp("    Hardy-Weinberg equilibrium")
        filteredHwe <- abs(snpSum$z.HWE) >= abs(qnorm(filters$hwe/2))
    }
    
    # Sample filters init
    disp("  Sample filters")
    filteredSampleCallRate <- filteredHetero <- filteredInbreed <- 
        logical(ncol(obj))
    
    # Sample filters: call rate
    if (!.isEmpty(filters$sampleCallRate)) {
        disp("    Call rate")
        filteredSampleCallRate <- sampleSum$Call.rate < filters$sampleCallRate
    }
    # Sample filters: heterozygosity
    if (!.isEmpty(filters$heteroHard) || !.isEmpty(filters$heteroStat)) {
        disp("    Heterozygosity")
        if (.isEmpty(filters$heteroHard)) {
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
            #filteredHetero <- sampleSum$Heterozygosity < filters$heteroHard
            filteredHetero <- sampleSum$Heterozygosity > filters$heteroHard
    }
    
    # Sample filters: inbreeding coefficient
    if (!.isEmpty(filters$inbreed)) {
        disp("    Inbreeding coefficient")
        hetF <- .calcInbreedFromSnpMatrix(x,snpSum,sampleSum)
        filteredInbreed <- hetF > filters$inbreed
    }
    
    # Report filtered entities with metadata
    disp("  Summarizing SNP filter results")
    filteredSnpCallInd <- which(filteredSnpCall)
    filteredMafInd <- which(filteredMaf)
    filteredHweInd <- which(filteredHwe)
    disp("  Summarizing sample filter results")
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
            ibd <- tryCatch({
                snpgdsIBDMoM(gdsHandle,kinship=TRUE,sample.id=gdsIds,
                    snp.id=snpsetIbd,num.thread=.coresFrac(rc))
                },error=function(e) {
                disp("Caught error ",e$message)
                disp("Closing GDS file connection")
                snpgdsClose(gdsHandle)
            },interrupt=function(i) {
                disp("Caught keyboard interruption!")
                disp("Closing GDS file connection")
                snpgdsClose(gdsHandle)
            },finally="")
        
        disp("Performing IBD selection...")
        ibdCoeff <- snpgdsIBDSelection(ibd)
        
        # Are there related samples?
        ibdCoeff <- ibdCoeff[ibdCoeff$kinship>=filters$IBD,,drop=FALSE]
        # If yes, mark for removal
        if (nrow(ibdCoeff) > 0) {
            sr <- unique(c(ibdCoeff$ID1,ibdCoeff$ID2))
            z <- match(sr,colnames(obj))
            ir <- z[!is.na(z)]
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
    tryCatch(snpgdsClose(gdsHandle),error=function(e) {
        closefn.gds(gdsHandle)
    },finally="")
    
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
        snpSub <- tryCatch({
            snpgdsLDpruning(handle,ld.threshold=ldCut,sample.id=samples,
                snp.id=snps,start.pos="first")
            #tmp <- capture.output({
            #   snpSub <- snpgdsLDpruning(handle,ld.threshold=filters$LD,
            #   sample.id=samples,snp.id=snps,start.pos="first")
            #})
            #disp(paste(tmp,collapse="\n"))
        },error=function(e) {
            disp("Caught error ",e$message)
            disp("Closing GDS file connection")
            snpgdsClose(handle)
        },interrupt=function(i) {
            disp("Caught keyboard interruption!")
            disp("Closing GDS file connection")
            snpgdsClose(handle)
        },finally="")
    
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
            "for PCA.")
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
        pco <- tryCatch({
            #tmp <- capture.output({
            if (!.isEmpty(npcs))
                snpgdsPCA(handle,sample.id=samples,snp.id=snps,
                    eigen.cnt=npcs,num.thread=.coresFrac(rc))
            else
                pco <- snpgdsPCA(handle,sample.id=samples,snp.id=snps,
                    num.thread=.coresFrac(rc))
            #})
        },error=function(e) {
            disp("Caught error ",e$message)
            disp("Closing GDS file connection")
            snpgdsClose(handle)
        },interrupt=function(i) {
            disp("Caught keyboard interruption!")
            disp("Closing GDS file connection")
            snpgdsClose(handle)
        },finally="")
    }
    return(pco)
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
        f <- f[given[check]]
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

.emptyFilters <- function(f) {
    empty <- unlist(lapply(f,.isEmpty))
    return(all(empty))
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
        
        message("")
    }
}
