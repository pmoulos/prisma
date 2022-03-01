# We must adjust the procedure for data derived from *exactly* the same arrays
# and manifest files etc.

mergeGWAS <- function(gwe1,gwe2,output=c("common","all"),mafResolve=NULL,
    gdsfile=NULL,writegds=FALSE) {
    # Input must be GWASExperiment
    if (!is(gwe1,"GWASExperiment") || !is(gwe2,"GWASExperiment"))
        stop("Both objects to be merged must be GWASExperiment objects!")
    
    # Check simple allele reverse resolution - reference based is safe when
    # data are derived e.g. from the *same* arrays
    if (!is.null(mafResolve) && !is.numeric(mafResolve))
        stop("mafResolve must be either NULL or numeric!")
    if (is.numeric(mafResolve) && !(mafResolve %in% c(1,2)))
        stop("When numeric, mafResolve must be either 1 or 2!")
    
    # Check if gdsfile location provided...
    if (writegds) {
        if (is.null(gdsfile)) {
            gdsfile <- tempfile()
            warning("It is advised to provide a location for the creation ",
                "of the new GDS file for better control...\nPlacing in ",
                "temporary ",gdsfile,immediate.=TRUE)
        }
        else if (!is.character(gdsfile))
            stop("gdsfile argument must be a valid path!")
    }
    
    output <- output[1]
    .checkTextArgs("Output SNPs (output)",output,c("common","all"),
        multiarg=FALSE)
    
    # Check if input objects have genomes
    g1 <- genome(gwe1)
    g2 <- genome(gwe2)
    if (is.na(g1) && is.na(g2)) {
        warning("A genome version is not provided and neither input ",
            "GWASExperiment object has an associated human genome ",
            "version!\n I will now try to guess...",immediate.=TRUE)
        g1 <-  guessHumanGenomeVersion(gwe1)
        g2 <-  guessHumanGenomeVersion(gwe2)
        if (!is.na(g1))
            genome(gwe1) <- g1
        if (!is.na(g2))
            genome(gwe2) <- g2
    }
    else if (!is.na(g1) && is.na(g2)) {
        warning("A genome version is not provided the second ",
            "GWASExperiment object does not have an associated human ",
            "genome version!\n I will now try to guess...",immediate.=TRUE)
        g2 <-  guessHumanGenomeVersion(gwe2)
        if (!is.na(g2))
            genome(gwe2) <- g2
    }
    else if (is.na(g1) && !is.na(g2)) {
        warning("A genome version is not provided the first ",
            "GWASExperiment object does not have an associated human ",
            "genome version!\n I will now try to guess...",immediate.=TRUE)
        g1 <-  guessHumanGenomeVersion(gwe1)
        if (!is.na(g1))
            genome(gwe1) <- g1
    }
    
    # Finally, must stop if not the same...
    if (g1 != g2)
        stop("When both input GWASExperiments have an associated ",
            "genome assembly, this must be the same!")
    
    # Check alleles on common SNPs
    disp("Now trying to merge features (markers, SNPs)")
    feat1 <- gfeatures(gwe1)
    feat2 <- gfeatures(gwe2)

    # Common markers to investigate alleles
    com <- intersect(rownames(feat1),rownames(feat2))
    disp("  Found ",length(com)," common SNPs")
    feat1Check <- feat1[com,]
    feat2Check <- feat2[com,]

    # Make investigation variables
    disp("  Checking common SNPs for reversed strands or alleles")
    invList <- .makeInvestigationVars(feat1Check,feat2Check)
    nm <- invList$nomatch
    investigate <- invList$investigate
    alleles <- invList$alleles

    # Can all be solved by simple strand flipping?
    disp("  Found ",length(nm)," SNPs with no mathcing alleles")
    disp("  Found ",length(investigate)," SNPs fixable by strand-flip")
    toResolve <- setdiff(seq_along(nm),investigate)
    
    # If not investigate further as there are major/minor allele discrepancies
    if (length(toResolve) > 0) { # Further investigation
        disp("  Found ",length(toResolve)," SNPs with potentially reversed ",
            "minor/major alleles")
        if (is.null(mafResolve)) {
            disp("    Trying to resolve online...")
            toReverse <- .resolveMinorAlleles(rownames(alleles)[toResolve],
                alleles)
        }
        else {
            disp("    Resolving using dataset ",mafResolve," as reference...")
            toReverse <- matrix(TRUE,length(rownames(alleles)[toResolve]),2)
            rownames(toReverse) <- rownames(alleles)[toResolve]
            toReverse[,mafResolve] <- FALSE
        }
        
        # Essentially the content of toReverse says:
        # If column 1 is TRUE, then reverse alleles in dataset 1 - not strand!
        # If column 2 is TRUE, then reverse alleles in dataset 2 - not strand!
        # We will also reverse genotypes later... This can be done with the same
        # matrix as it is named.

        # Begin allele/strand harmonization process...
        # First, we reverse the alleles that must be reversed in each dataset
        if (any(toReverse)) {
            disp("  ",length(which(apply(toReverse,1,any)))," SNPs can be ",
                "resolved, excluding ambiguous")
            a1ind <- rownames(toReverse)[toReverse[,1]]
            feat1Check <- .revMapAllele(feat1Check,a1ind)
            a2ind <- rownames(toReverse)[toReverse[,2]]
            feat2Check <- .revMapAllele(feat2Check,a2ind)
            # toResolve will be later used also for genotype reverse

            # Then proceed with 1st round of strand-flip. Which dataset should  
            # be flipped? The one with the most reversed alleles OR the one
            # which is not reference
            if (is.null(mafResolve)) {
                tab <- apply(toReverse,2,table)
                flipInd <- which(tab["TRUE",]==max(tab["TRUE",]))
            }
            else
                flipInd <- ifelse(mafResolve==1,2,1)
        }
        else # If cannot be decided like that, reverse the smaller
            flipInd <- ifelse(nrow(gwe2) > nrow(gwe1),1,2)
        disp("  Based on reverse allele majority or dataset SNP content, I ",
            "will strand-flip common SNPs from input dataset #",flipInd)
        
        disp("  Strand-flipping")
        if (flipInd == 1)
            feat1Check <- .revMapAllele(feat1Check,rownames(alleles),"strand")
        else
            feat2Check <- .revMapAllele(feat2Check,rownames(alleles),"strand")
        
        # After allele reversing and 1st round of strand-flipping, we must
        # probably re-investigate which alleles can be further matched by
        # reverse strand-fliping as it may be an effect of allele reversing
        # in 1st place. This should be done if not identical feat1Check,
        # feat2Check. In theory, most of them can now safely be flipped
        # according to flipInd and maybe a very few ambiguous are left.
        if (!identical(feat1Check,feat2Check)) {
            disp("  Fixing SNPs that were strand-flipped twice")
            invList2 <- .makeInvestigationVars(feat1Check,feat2Check,ambig=TRUE)
            nm2 <- invList2$nomatch
            investigate2 <- invList2$investigate
            alleles2 <- invList2$alleles
            
            # 2nd round
            if (flipInd == 1)
                feat1Check <- .revMapAllele(feat1Check,rownames(alleles2),
                    "strand")
            else
                feat2Check <- .revMapAllele(feat2Check,rownames(alleles2),
                    "strand")
        }
    }

    # If still not identical and alleles are different, we cannot do many 
    # things... These must be really few and are dropped
    disp("  Checking for non-resolvable alleles")
    bad <- which(feat1Check$allele.1 != feat2Check$allele.1 |
        feat1Check$allele.2 != feat2Check$allele.2)
    theBad <- NULL
    if (length(bad) > 0) {
        theBad <- rownames(feat1Check)[bad]
        disp("    Found ",length(bad),". These are:\n",paste(theBad,
            collapse=", "))
        disp("  They may be multiallelic SNPs mapped different on each set.")
        disp("  You may investigate further but I will have to drop them...")
        feat1Check <- feat1Check[-bad,]
        feat2Check <- feat2Check[-bad,]
    }

    # They should now be identical... If not, there may be slight differrences 
    # in positions... We have even seen different chromosomes! In this case we 
    # are using rsLocsFromEnsembl as rsnps does not support older assemblies
    if (!identical(feat1Check,feat2Check)) {
        z <- which(feat1Check$chromosome != feat2Check$chromosome
            | feat1Check$position != feat2Check$position)
        if (length(z) > 0) {
            disp("  Found ",length(z)," SNPs with non-conformable locations")
            disp("    Trying to resolve from dbSNP...")
            missAnn <- rownames(feat1Check)[z]
            # gv will be taken from GWASExperiment object
            missAnnInfo <- rsLocsFromEnsembl(missAnn,gv=genome(gwe1))
            # Have we retrieved NA positions? If yes they cannot be used
            missAnnInfo <- 
                missAnnInfo[!is.na(missAnnInfo$chrom_start),,drop=FALSE]
            # Rearrange
            rownames(missAnnInfo) <- missAnnInfo$refsnp_id
            missAnnInfo <- missAnnInfo[missAnn,,drop=FALSE]
            # It's possible that we now have NAs
            missAnnInfo <- 
                missAnnInfo[!is.na(missAnnInfo$chrom_start),,drop=FALSE]
            
            # Now harmonize
            isDF <- FALSE
            if (is(feat1Check,"DataFrame")) {
                isDF <- TRUE
                feat1Check <- as.data.frame(feat1Check)
                feat2Check <- as.data.frame(feat2Check)
            }
                
            feat1Check[rownames(missAnnInfo),"chromosome"] <- 
                missAnnInfo$chr_name
            feat2Check[rownames(missAnnInfo),"chromosome"] <- 
                missAnnInfo$chr_name
            feat1Check[rownames(missAnnInfo),"position"] <- 
                missAnnInfo$chrom_start
            feat2Check[rownames(missAnnInfo),"position"] <- 
                missAnnInfo$chrom_start
            
            if (isDF) {
                feat1Check <- DataFrame(feat1Check)
                feat2Check <- DataFrame(feat2Check)
            }
        }
    }

    # If still not identical, then purge SNPs that cannot be mapped in any way
    if (!identical(feat1Check,feat2Check)) {
        z <- which(feat1Check$chromosome != feat2Check$chromosome
                | feat1Check$position != feat2Check$position)
        if (length(z) > 0) { # Should be...
            disp("  Common SNP set still not identical! ",length(z)," SNPs ",
                "cannot be resolved by any means! Purging...")
            theBad <- c(theBad,feat1Check$snp.name[z])
            feat1Check <- feat1Check[-z,]
            feat2Check <- feat2Check[-z,]
        }
    }

    # If still not identical, then stop
    if (!identical(feat1Check,feat2Check))
        stop("Cannot merge data without a common reference! Sorry...")

    # Dataset-only SNPs... No really need to flip anything as we are creating a 
    # new dataset... Or flip strand for consistency? No because some SNPs may 
    # or may not be flipped across array designs
    if (output == "all") {
        disp("  Retrieving SNPs specific to each dataset")
        feat1Only <- feat1[setdiff(rownames(feat1),rownames(feat1Check)),]
        feat2Only <- feat2[setdiff(rownames(feat2),rownames(feat2Check)),]
    }

    # Now, we must merge and order the two feature sets
    disp("  Merging features")
    if (output == "all")
        featMerged <- rbind(feat1Only,feat1Check,feat2Only)
    else
        featMerged <- feat1Check
    featMerged <- featMerged[order(featMerged$chromosome,featMerged$position),]
    
    # Merge samples... should be easy... PLINK writing will show
    disp("Now trying to merge samples")
    sampMerged <- rbind(gsamples(gwe1),gsamples(gwe2))
    
    # Now, we must move on constructing the unified genotypes... toReverse 
    # matrix will be used for genotype reversing
    disp("Now trying to construct unified genotypes")
    geno1 <- genotypes(gwe1)
    geno2 <- genotypes(gwe2)
    # transpose for use snpStats::switch.alleles
    geno1Check <- geno1[rownames(feat1Check),]
    geno2Check <- geno2[rownames(feat2Check),]
    # Switch and retranspoe
    if (any(toReverse)) {
        disp("  Reversing alleles that were resolved during SNP matching")
        geno1Check <- t(switch.alleles(t(geno1Check),intersect(a1ind,
            rownames(geno1Check))))
        geno2Check <- t(switch.alleles(t(geno2Check),intersect(a2ind,
            rownames(geno2Check))))
    }
    # Merge the commons
    disp("  Merging genotypes")
    genoCommon <- cbind(geno1Check,geno2Check)
    if (output == "all") {
        # Now the rest from each dataset - only genotypes + NAs for the 
        # rest samples
        geno1Only <- geno1[rownames(feat1Only),]
        geno2Only <- geno2[rownames(feat2Only),]
        geno1Add <- SnpMatrix(matrix(0,nrow(geno1Only),ncol(geno2Only)))
        rownames(geno1Add) <- rownames(geno1Only)
        colnames(geno1Add) <- colnames(geno2Only)
        geno2Add <- SnpMatrix(matrix(0,nrow(geno2Only),ncol(geno1Only)))
        rownames(geno2Add) <- rownames(geno2Only)
        colnames(geno2Add) <- colnames(geno1Only)
        geno1Both <- cbind(geno1Only,geno1Add)
        geno2Both <- cbind(geno2Add,geno2Only)
        # Merge!
        genoMerged <- rbind(geno1Both,genoCommon,geno2Both)
    }
    else
        genoMerged <- genoCommon
        
    # Align...
    if (!identical(rownames(featMerged),rownames(genoMerged)))
        genoMerged <- genoMerged[rownames(featMerged),]
    if (!identical(colnames(featMerged),rownames(sampMerged)))
        genoMerged <- genoMerged[,rownames(sampMerged)]
    
    # And now finally return a new GWASExperiment
    disp("Done!")
    gwem <- GWASExperiment(
        genotypes=genoMerged,
        features=featMerged,
        samples=sampMerged,
        metadata=list(
            genome="hg19", # Should be taken from gwe1,2
            backend="snpStats",
            filters=.initFilterInfo(),
            gdsfile=gdsfile,
            alleleOrder="plink"
        )
    )
    
    # It's possible somehow that the bad SNPs have escaped...
    if (length(theBad) > 0) {
        #print(theBad)
        nexpr <- paste0("^(",paste(theBad,collapse="|"),")$")
        wtf <- grep(nexpr,rownames(gwem),perl=TRUE)
        if (length(wtf) > 0)
            gwem <- gwem[-wtf,]
    }
    
    if (writegds) {
        disp("Writing also GDS file...")
        GWASExperiment2GDS(gwem)
        disp("Done!")
    }
    disp(" ")
    
    return(gwem)
}

.revMapAllele <- function(x,n,w=c("allele","strand")) {
    w <- w[1]
    
    inDF <- ifelse(is(x,"DataFrame"),TRUE,FALSE)
    if (inDF)
        x <- as.data.frame(x)
    
    if (w == "allele") {
        a1 <- x[n,"allele.1"]
        a2 <- x[n,"allele.2"]
        x[n,"allele.1"] <- a2
        x[n,"allele.2"] <- a1
    }
    else if (w == "strand") {
        x[n,"allele.1"] <- unname(.revBase(x[n,"allele.1"]))
        x[n,"allele.2"] <- unname(.revBase(x[n,"allele.2"]))
    }
    
    if (inDF)
        return(DataFrame(x))
    else
        return(x)
}

# Then how we decide on flip of the rest cases? Based on the flipping of the
# common cases... In this case, we mostly flip the NAFLD array

# Minor allele flip resolution function
# 1. Query SNPs to resolve with biomaRt - genome version not important as we
#    are not interested on locations, only alleles.
# 2. If minor allele is available from Ensembl, we make the flip decision based
#    on it. If not (e.g. outdated/merged SNPs, we query dbSNP with rsnps 
#    package) and fill the Ensembl result.
# 3. Some SNPs may have been removed from dbSNP. Nothing can be done. Flip for
#    these will be based on the majority of the rest.
# 4. Some found SNPs may not have information on MAF and minor alleles. We
#    requery dbSNP with rsnps and again we are trying to fill missing info.
#    Again, some will not be resolved and will be flipped according to majority.
.resolveMinorAlleles <- function(snps,alleles) {
    if (!requireNamespace("rsnps"))
        stop("R package rsnps is required!")
    
    # We do not care about the genome as we are not looking for locations, so
    # we leave it at hg38 (default)
    disp("      using Ensembl...")
    snpInfo <- rsLocsFromEnsembl(snps)
    rownames(snpInfo) <- snpInfo$refsnp_id
    
    toReverse <- matrix(FALSE,length(snps),2)
    rownames(toReverse) <- snps
    
    # Some SNPs may be outdated - try rsnps package
    if (nrow(snpInfo) < length(snps)) {
        disp("      using dbSNP...")
        miss <- setdiff(snps,rownames(snpInfo))
        missHits <- tryCatch({
            suppressWarnings(rsnps::ncbi_snp_query(miss))
        },error=function(e) {
            disp("Caught error: ",e$message)
            disp("Probably SNP does not currentl exist in dbSNP! Will drop...")
            return(data.frame())
        },finally="")
        if (nrow(missHits) > 0) { # Construct a df to add to snpInfo
            toAdd <- data.frame(
                chr_name=missHits$chromosome,
                chrom_start=missHits$bp,
                refsnp_id=missHits$query,
                allele=gsub(",","/",missHits$alleles),
                minor_allele=.minorFromRsnpsHits(missHits,"allele"),
                minor_allele_freq=.minorFromRsnpsHits(missHits,"freq"),
                chrom_strand=1,
                row.names=missHits$query
            )
            snpInfo <- rbind(snpInfo,toAdd)
        }
    }
    
    # Check if we are ok, otherwise we need to exclude SNPs that cannot be
    # matched anywhere
    if (nrow(snpInfo) < length(snps)) {
        stillMiss <- setdiff(snps,rownames(snpInfo))
        jj <- which(snps %in% stillMiss)
        snps <- snps[-jj]
    }
    
    # Reorder
    snpInfo <- snpInfo[snps,,drop=FALSE]
    
    # Now look, if we have still MAFs missing so as to query dbSNP again
    noMafInd <- which(snpInfo$minor_allele == "")
    if (length(noMafInd) > 0) {
        miss <- rownames(snpInfo)[noMafInd]
        missHits <- suppressWarnings(rsnps::ncbi_snp_query(miss))
        if (nrow(missHits) > 0) {
            toReplace <- data.frame(
                chr_name=missHits$chromosome,
                chrom_start=missHits$bp,
                refsnp_id=missHits$query,
                allele=gsub(",","/",missHits$alleles),
                minor_allele=.minorFromRsnpsHits(missHits,"allele"),
                minor_allele_freq=.minorFromRsnpsHits(missHits,"freq"),
                chrom_strand=1,
                row.names=missHits$query
            )
            snpInfo[rownames(toReplace),] <- toReplace
        }
    }
    
    # We 've done all that we could... Time to try and resolve
    tmp <- lapply(snps,function(x) {
        # Get the recorded minor alleles from two datasets
        dma <- alleles[x,c(1,3)]
        # Get the "official" minor alleles from dbSNP
        oma <- snpInfo[x,"minor_allele"]
        
        # The other one than the one below, should be allele-flipped - therefore
        # negation
        return(which(!(dma %in% c(oma,.revBase(oma)))))
    })
    names(tmp) <- snps
    
    # Some may remain unresolvable... They will be flipped per the majority
    unres <- which(lengths(tmp)>1)
    
    # Fill resolution table
    tmp <- unlist(tmp[lengths(tmp)==1])
    toReverse[names(tmp)[tmp==1],1] <- TRUE
    toReverse[names(tmp)[tmp==2],2] <- TRUE
    
    # Decide on unresolvable flips
    if (length(unres) > 0) {
        tab <- apply(toReverse,2,table)
        toReverse[names(unres),which(tab["TRUE",]==max(tab["TRUE",]))] <- TRUE
    }
    
    return(toReverse)
}

.minorFromRsnpsHits <- function(x,out=c("allele","freq")) {
    out <- out[1]
    if (out == "allele") {
        o <- x$minor
        # May be NAs so logical indexing fails
        ii <- which(x$maf > 0.5)
        if (length(ii) > 0)
            o[ii] <- x$ref_seq[ii]
    }
    else if (out == "freq") {
        o <- x$maf
        ii <- which(o > 0.5)
        if (length(ii) > 0)
            o[ii] <- 1 - o[ii]
    }
    return(o)
}

.revBase <- function(x) {
    return(c(
        C="G",
        G="C",
        A="T",
        `T`="A"
    )[x])
}

.makeInvestigationVars <- function(x1,x2,ambig=FALSE) {
    # Let's find where alleles are not the same
    nm <- which((x1$allele.1 != x2$allele.1) 
        | (x1$allele.2 != x2$allele.2))

    # Investigate - init a matrix for later use
    inv <- data.frame(
        gn.a1=x1$allele.1[nm],
        gn.a2=x1$allele.2[nm],
        gt.a1=x2$allele.1[nm],
        gt.a2=x2$allele.2[nm]
    )
    rownames(inv) <- rownames(x1)[nm]

    # Collapse inv
    colgts <- paste(inv[,1],inv[,2],inv[,3],inv[,4],sep="")
    names(colgts) <- rownames(inv)

    # Also initialize a logical matrix with all non-same SNPs which should be
    # strand flipped or genotype reversed.


    # These will have to be strand-flipped only, which will be flipped will
    # depend the dataset with less flips of the process below. The other one
    # will be fliped.
    investigate <- list()
    # If A->C become T->G OK
    investigate$actg <- which(colgts=="ACTG")
    # If A->G become T->C OK
    investigate$agtc <- which(colgts=="AGTC")
    # If C->A become G->T OK
    investigate$cagt <- which(colgts=="CAGT")
    # If C->T become G->A OK
    investigate$ctga <- which(colgts=="CTGA")
    # If G->A become C->T OK
    investigate$gact <- which(colgts=="GACT")
    # If G->T become C->A OK
    investigate$gtca <- which(colgts=="GTCA")
    # If T->C become A->G OK
    investigate$tcag <- which(colgts=="TCAG")
    # If T->G become A->C OK
    investigate$tgac <- which(colgts=="TGAC")
    
    if (ambig) {
        # If A->T become T->A OK? Ambiguous
        investigate$atta <- which(colgts=="ATTA")
        # If C->G become G->C OK? Ambiguous
        investigate$cggc <- which(colgts=="CGGC")
        # If G->C become C->G OK? Ambiguous
        investigate$gccg <- which(colgts=="GCCG")
        # If T->A become A->T OK? Ambiguous
        investigate$taat <- which(colgts=="TAAT")
    }
    
    # Merge them
    investigate <- Reduce("union",investigate)
    
    return(list(nomatch=nm,investigate=investigate,alleles=inv))
}
