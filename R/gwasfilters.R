# input may be a prefix of bim, bed, fam or a list with each
importGWAS <- function(input,phenos=NULL,backend=c("snpStats","bigsnpr"),
    selection=NULL,genome=NA_character_,
    gdsfile=ifelse(backend=="snpStats",tempfile(),NA)) {
    backend <- backend[1]
    
    if (backend == "snpStats") {
        if (!requireNamespace("snpStats"))
            stop("Bioconductor package snpStats is required!")
        if (!requireNamespace("SNPRelate"))
            stop("Bioconductor package SNPRelate is required!")
        if (is.na(gdsfile) || !is.character(gdsfile)) {
            warning("A valid path to store GDS file for filtering is required!",
                " Assuming default (tempfile())...")
            gdsfile <- tempfile()
        }
    }
    if (backend == "bigsnpr" && !requireNamespace("snpStats"))
        stop("Bioconductor package bigsnpr is required!")
    
    selection <- .checkSelection(selection)
    
    if (is.character(input)) {
        if (!dir.exists(basename(input)))
            stop("The input directory containing PLINK files does not exist!")
        input <- list(bed=paste0(input,".bed"),bim=paste0(input,".bim"),
            bed=paste0(input,".fam"))
    }
    plExist <- vapply(input,file.exists,logical(1))
    if (!all(plExist))
        stop("Input PLINK file(s) ",paste(input[!plExist],collapse=", "),
            " do not exist!")
    
    if (backend == "snpStats") {
        message("Reading PLINK files with snpStats framework")
        snpObj <- read.plink(input$bed,input$bim,input$fam,
            select.subjects=selection$samples,select.snps=selection$snps)
        message("Reading PLINK files with SNPRelate framework and storing ",
            "output to ",gdsfile)
        snpgdsBED2GDS(input$bed,input$fam,input$bim,gdsfile,family=TRUE)
        # We also need to read a GDS file for LD and IBD filtering
        return(GWASExperiment(
            genotypes=t(snpObj$genotypes),
            features=snpObj$map,
            samples=snpObj$fam,
            phenotypes=phenos,
            metadata=list(
                genome=genome,
                backend=backend,
                filters=.initFilterInfo(),
                gdsfile=gdsfile
            )
        ))
    }
    else if (backend == "bigsnpr") {
        # TODO: Add bigsnpr options
        #snpObj <- ... snp_readSth
    }
        
}

filterGWAS <- function(obj,filters=getDefaults("filters"),verbose=TRUE,
    rc=NULL) {
    filters <- .checkFilters(filters)
    
    # Later input type may be more native, e.g. read from CSV files
    if (is(assay(obj,1),"SnpMatrix"))
        return(.filterWithSnpStats(obj,filters,verbose,rc))
    else if (is(assay(obj,1),"bigsnp"))
        return(.filterWithBigSnpr(x,filters,verbose))
}

.filterWithSnpStats <- function(obj,filters,verbose=TRUE,rc=NULL) {
    message("Performing basic filtering")
    objBas <- .filterWithSnpStatsBasic(obj,filters)
    
    message("\nPerforming sample IBD filtering and classic PCA")
    objIbd <- .filterWithSnpStatsIbd(objBas,filters,rc)
    
    message("\nPerforming robust PCA filtering for automatic outlier detection")
    objPca <- .filterWithSnpStatsRobustPca(objIbd,filters)
    
    if (verbose) {
        filtDf <- filters(objIbd)
        if (nrow(filtDf) > 0) {
            message("\nFiltered SNPs:")
            message("  SNP call rate             : ",
                filtDf["snpCallRate","filtered"])
            message("  Minor Allele Frequency    : ",filtDf["maf","filtered"])
            message("  Hardy-Weinberg equilibrium: ",filtDf["hwe","filtered"])
            
            message("\nFiltered samples:")
            message("  Sample call rate          : ",
                filtDf["sampleCallRate","filtered"])
            message("  Heterozygosity            : ",
                filtDf["heteroHard","filtered"])
            message("  Inbreed                   : ",
                filtDf["inbreed","filtered"])
            message("  Identity By Descent (IBD) : ",filtDf["IBD","filtered"])
            message("  Robust PCA                : ",
                filtDf["pcaOut","filtered"])
        }
    }
    
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
    currFiltDf <- filters(obj)
    if (nrow(currFiltDf) == 0)
        filters(obj) <- filtDf
    else
        filters(obj) <- rbind(currFiltDf,filtDf)
    
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
    message("  reading GDS file ",m$gdsfile,"...")
    gdsHandle <- openfn.gds(m$gdsfile,readonly=FALSE)
    gdsIds <- read.gdsn(index.gdsn(gdsHandle,"sample.id"))
    gdsIds <- sub("-1","",gdsIds)
    add.gdsn(gdsHandle,"sample.id",gdsIds,replace=TRUE)
    
    # LD pruning
    if (!is.na(filters$LD)) {
        message("  performing LD pruning...\n")
        x <- t(assay(obj,1))
        if (.testing)
            snpSub <- snpgdsLDpruning(gdsHandle,ld.threshold=filters$LD,maf=0.1,
                autosome.only=FALSE,sample.id=rownames(x),snp.id=colnames(x))
        else
            # Default is start.pos="random" but this requires a seed for reprod
            snpSub <- snpgdsLDpruning(gdsHandle,ld.threshold=filters$LD,
                sample.id=rownames(x),snp.id=colnames(x),start.pos="first")
                
        snpsetIbd <- unlist(snpSub,use.names=FALSE)
        message("  LD pruning finished! ",length(snpsetIbd)," SNPs will be ",
            "used for IBD analysis.")
    }
    else
        snpsetIbd <- colnames(x)
    
    # IBD analysis
    ir <- NULL
    sr <- rownames(x)
    if (!is.na(filters$IBD)) {
        message("  performing IBD calculation...")
        if (.testing)
            ibd <- snpgdsIBDMoM(gdsHandle,kinship=TRUE,sample.id=gdsIds,
                snp.id=snpsetIbd,num.thread=.coresFrac(rc),maf=0.1,
                autosome.only=FALSE)
        else
            ibd <- snpgdsIBDMoM(gdsHandle,kinship=TRUE,sample.id=gdsIds,
                snp.id=snpsetIbd,num.thread=.coresFrac(rc))
        message("  performing IBD selection...")
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
        currFiltDf <- filters(obj)
        if (nrow(currFiltDf) == 0)
            filters(obj) <- filtDf
        else
            filters(obj) <- rbind(currFiltDf,filtDf)
    }
           
    # Finally, perform a PCA with SNPRelate using IBD filtered samples
    message("  performing PCA...")
    if (.testing)
        pca <- snpgdsPCA(gdsHandle,maf=0.1,autosome.only=FALSE,sample.id=sr,
            snp.id=snpsetIbd)
    else
        pca <- snpgdsPCA(gdsHandle,sample.id=sr,snp.id=snpsetIbd,
            num.thread=.coresFrac(rc))
    
    m <- metadata(obj)
    m$pcaOut <- pca
    metadata(obj) <- m
    
    # Close GDS
    closefn.gds(gdsHandle)
    
    if (!is.null(ir))
        return(obj[,-ir,drop=FALSE])
    else
        return(obj)
    
    #pctab <- data.frame(sampleId=pca$sample.id,PC1=pca$eigenvect[,1],
    #   PC2=pca$eigenvect[,2],stringsAsFactors=FALSE)
}

.filterWithSnpStatsRobustPca <- function(obj,filters) {
    if (!requireNamespace("rrcov"))
        stop("R package rrcov is required for robust PCA!")
    
    if (!filters$pcaOut)
        return(obj)
    
    method <- filters$pcaRobust
    npc <- filters$nPC
    obj <- .robustPcaWithSnpStats(obj,method,npc)
    m <- metadata(obj)
    jj <- which(!m$pcaRob$flag)
    
    # Update filters
    filtDf <- data.frame(parameter="pcaOut",name="Robust PCA",
        value=filters$pcaRobust,type="Sample",filtered=length(jj),
        row.names="pcaOut")
    currFiltDf <- filters(obj)
    if (nrow(currFiltDf) == 0)
        filters(obj) <- filtDf
    else
        filters(obj) <- rbind(currFiltDf,filtDf)
    
    # Return possibly filtered object
    if (length(jj) > 0)
        return(obj[,-jj,drop=FALSE])
    else
        return(obj)
}

.robustPcaWithSnpStats <- function(obj,method=c("grid","hubert"),npc=0) {
    m <- metadata(obj)
    if ("pcaOut" %in% names(m)) {
        snps <- m$pcaOut$snp.id
        if (npc == 0)
            npc <- ncol(m$pcaOut$eigenvect)
    }
    else
        snps <- rownames(obj)
    
    x <- assay(obj,1)
    y <- as(x,"numeric")
    if (any(is.na(y))) # Not imputed, assume 0
        y[is.na(x)] <- 0
    y <- y[snps,,drop=FALSE]
    
    if (method == "grid") {
        message("  with grid search method")
        P <- PcaGrid(t(y),k=npc)
    }
    else if (method == "hubert") {
        message("  with Hubert method")
        kmax <- ifelse(npc>10,npc,10)
        P <- PcaGrid(t(y),k=npc,kmax=kmax)
    }
    
    m$pcaRob <- P
    metadata(obj) <- m
    
    return(obj)
}

.filterWithBigSnpr <- function(obj,filters,verbose=TRUE) {
    # Stub
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

.checkSelection <- function(s) {
    if (is.null(s))
        return(s)
    
    if (!is.list(s))
        stop("The selection argument must be a list!")
    if (!(all(names(s) %in% c("samples","snps"))))
        stop("The selection argument must be a named list with members ",
            "'samples' and 'snps'!")
            
    problemFound <- vapply(s,function(x) {
        if (!is.null(x)) {
            if (is.numeric(x) && all(x>0))
                return(FALSE)
            else
                return(TRUE)
        }
        else
            return(FALSE)
    },logical(1))
            
    if (any(problemFound))
        stop("All selection list members must be non-negative integers ",
            "denoting SNP matrix rows\nand/or columns! ",
            paste(names(s)[problemFound],collapse=", ")," are not...")
    
    return(s)
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

