# I think that for large queries we must download and cache the raw tab
# delimited file... maybe not... takes some but gwasrapidd works!
prismaLookup <- function(prismaOut) {
    # Required packages
    if (!requireNamespace("gwasrapidd"))
        stop("R package gwasrapidd is required!")
    
    # Determine if we have a prisma output complete object or a list of 
    # candidate PRSs
    isPrismaOut <- is.list(prismaOut) && 
        all(c("params","results") %in% names(prismaOut))
    if (!isPrismaOut) # Check the list of PRSs as it may be manual
        prismaOut <- .checkPrismaSelectionsOutput(prismaOut)
    
    if (isPrismaOut) {
        lookupOut <- vector("list",length(prismaOut$results))
        names(lookupOut) <- names(prismaOut$results)
        for (m in names(prismaOut$results)) {
            disp("-----> Querying for ",m)
            lookupOut[[m]] <- 
                .prismaLookupWorker(prismaOut$results[[m]]$candidates)
        }
        return(lookupOut)
    }
    else
        return(list(.prismaLookupWorker(prismaOut)))
}

.prismaLookupWorker <- function(selections) {
    # Check the format of the input list
    selections <- .checkPrismaSelectionsOutput(selections)
    
    # Empty trait data frame
    ..emptyLookupDf <- function() {
        nil <- data.frame(variant_id="v",risk_allele="N",association_id="1",
            gene_name="N",efo_id="E",trait="T",freq="1")
        return(nil[-1,])
    }
    
    # Find the largest and get SNP names
    size <- as.numeric(names(selections))
    forLook <- rownames(selections[[which(size == max(size))]])
    
    # Make the API call - turn off warnings when SNP not found
    disp("Querying GWAS Catalog with the largest PRS")
    A <- gwasrapidd::get_associations(variant_id=forLook,warnings=FALSE)
    
    # We need to create a data frame which summarizes the findings... Then we
    # attach them to each smaller size PRS candidate by intersecting
    if (nrow(A@associations) > 0) {
        # Need to collapse gene names per association id
        ge <- as.data.frame(A@genes)
        ge <- split(ge,factor(ge$association_id,
            levels=unique(ge$association_id)))
        #ge <- unlist(lapply(ge,function(x) {
        #    paste0(x$gene_name,collapse=", ")
        #}))
        ge <- lapply(ge,function(x) {
            if (all(is.na(x$gene_name)))
                return(NA)
            else if (any(is.na(x$gene_name)))
                return(paste0(x$gene_name[!is.na(x$gene_name)],collapse=", "))
            else
                return(paste0(x$gene_name,collapse=", "))
        })
        
        # For each association, we need to run get_traits...
        tmpT <- tmpE <- vector("list",length(A@associations$association_id))
        names(tmpT) <- names(tmpE) <- A@associations$association_id
        for (a in A@associations$association_id) {
            disp("  Querying association ",a)
            tmp <- gwasrapidd::get_traits(association_id=a)
            tmpE[[a]] <- paste0(tmp@traits$efo_id,collapse=", ")
            tmpT[[a]] <- paste0(tmp@traits$trait,collapse=", ")
        }
        
        # The found hits
        #hits <- data.frame(
        #    association_id=A@associations$association_id,
        #    variant_id=A@risk_alleles$variant_id,
        #    risk_allele=A@risk_alleles$risk_allele,
        #    gene_name=ge,
        #    efo_id=unlist(tmpE,use.names=FALSE),
        #    trait=unlist(tmpT,use.names=FALSE)
        #)
        
        tmpdf <- data.frame(
            variant_id=A@risk_alleles$variant_id,
            risk_allele=A@risk_alleles$risk_allele,
            association_id=A@risk_alleles$association_id
        )
        tmps <- split(tmpdf,factor(tmpdf$association_id,
            levels=unique(tmpdf$association_id)))
        tmps <- lapply(names(tmps),function(n,D) {
            if (!is.null(ge[[n]]))
                D[[n]]$gene_name <- rep(ge[[n]],nrow(D[[n]]))
            else
                D[[n]]$gene_name <- rep(NA,nrow(D[[n]]))
            
            if (!is.null(tmpE[[n]]))
                D[[n]]$efo_id <- rep(tmpE[[n]],nrow(D[[n]]))
            else
                D[[n]]$efo_id <- rep(NA,nrow(D[[n]]))
            
            if (!is.null(tmpE[[n]]))
                D[[n]]$trait <- rep(tmpT[[n]],nrow(D[[n]]))
            else
                D[[n]]$trait <- rep(NA,nrow(D[[n]]))
            
            return(D[[n]])
        },tmps)
        
        hits <- do.call("rbind",tmps)
        
        # Now we have to split the hits according to each candidate PRS
        partHits <- lapply(selections,function(x) {
            m <- which(hits$variant_id %in% rownames(x))
            if (length(m) == 0)
                return(..emptyLookupDf())
            else {
                y <- hits[m,,drop=FALSE]
                z <- x[y$variant_id,"freq"]
                y$freq <- z[!is.na(z)]
                return(y)
            }
        })
        
        return(partHits)
    }
    else {
        o <- list(..emptyLookupDf())
        names(o) <- names(selections)
        return(o)
    }
}

getGWASVariants <- function(efoId=NULL,efoTrait=NULL,removeUnknownRisk=TRUE,
    retries=5) {
    # Package checking, meaningless to continue
    if (!requireNamespace("gwasrapidd"))
        stop("R package gwasrapidd is required!")
    if (!requireNamespace("jsonlite"))
        stop("R package jsonlite is required!")
    
    # Input checking and sanitization
    if ((missing(efoId) || is.null(efoId)) 
        && (missing(efoTrait) || is.null(efoTrait)))
        stop("At least one of GWAS EFO ID (efoId) or Trait (efoTrait) must be ",
            "provided!")
    
    # Multiple EFO ids may have been provided... Check them
    if (!is.null(efoId))
        .checkEFOFormat(efoId)
    # Check also if traits are characters... No other checks are the terms are
    # specific and case sensitive
    if (!is.null(efoTrait) && !is.character(efoTrait))
        stop("efoTrait must be a character vector!")
    
    # Calls
    disp("Scheduling ",length(efoId) + length(efoTrait)," GWAS API calls...")
    
    # EFO id calls
    if (!is.null(efoId))
        resultId <- .iterGWASVariantsAPICall(input=efoId,type="id",
            retries=retries)
    else
        resultId <- .emptyVariantsDf()
    
    # EFO trait calls
    if (!is.null(efoTrait))
        resultTr <- .iterGWASVariantsAPICall(input=efoTrait,type="trait",
            retries=retries)
    else
        resultTr <- .emptyVariantsDf()
    
    # Join the results
    results <- rbind(do.call("rbind",resultId),do.call("rbind",resultTr))
    
    # Remove duplicates if any
    results <- results[!duplicated(results),,drop=FALSE]
    rownames(results) <- results$variant_id
    
    # Remove SNPs where association with risk allele is not defined (if any)
    if (removeUnknownRisk)
        results <- results[!is.na(results$risk_allele),,drop=FALSE]
        
    return(results)
}

.iterGWASVariantsAPICall <- function(input,type=c("id","trait"),retries=5) {
    type <- type[1]
    N <- length(input)
    result <- vector("list",N)
    names(result) <- input
    
    suff <- ifelse(type=="id"," EFO ID(s)","EFO trait(s)")
    disp("  Requesting ",N," GWAS API calls with ",suff,"...")
    
    complete <- FALSE
    times <- 1
    while(!complete && times<=retries) {
        failed <- logical(N)
        names(failed) <- input
        calls <- 0
        
        for (id in input) {
            calls <- calls + 1
            disp("    Call ",calls,": ",id)
            tmp <- .GWASVariantsWorker(input=id,type=type)
            if (tmp$success) # Both API calls successful
                result[[id]] <- tmp$data
            else # Which id failed?
                failed[id] <- TRUE
        }
        
        complete <- !any(failed)
        if (!complete) {
            input <- names(failed)[failed]
            disp("  Failed API calls for ",suff,":")
            disp("    ",paste(input,collapse=", "))
            disp("  Will retry ",retries-times," times...")
        }
        times <- times + 1
    }
    
    return(result)
}

.GWASVariantsWorker <- function(input=NULL,type=c("id","trait")) {
    input <- input[1] # Protection from accidental pass more ids/traits
    type <- type[1]
    
    # Code templates
    if (type == "id") {
        exv <- paste0('get_variants(efo_id="',input,'",warning=FALSE)')
        exa <- paste0('get_associations(efo_id="',input,'",warning=FALSE)')
    }
    else if (type == "trait") {
        exv <- paste0('get_variants(efo_trait="',input,'",warning=FALSE)')
        exa <- paste0('get_associations(efo_trait="',input,'",warning=FALSE)')
    }
    
    # Get the variants
    efoVariants <- tryCatch({
        eval(parse(text=exv))
    },error=function(e) {
        disp("Possible connection failure! Marking...")
        disp("Caught error: ",e$message)
        return(toJSON(list(type=type,input=input),auto_unbox=TRUE))
    })
    
    # Firstly, we need the get_variants call... If failed, we don't proceed to
    # get_associations and we mark the whole call as failed
    if (!is(efoVariants,"json")) {
        if (nrow(efoVariants@variants) == 0) # Nothing found - OK
            return(list(success=TRUE,data=.emptyVariantsDf()))
        else # Call the get_assocations API
            efoAssoc <- tryCatch({
                eval(parse(text=exa))
            },error=function(e) {
                disp("Possible connection failure! Marking...")
                disp("Caught error: ",e$message)
                return(toJSON(list(type=type,input=input),auto_unbox=TRUE))
            },finally="")
        
        # Now, if get_associations call failed, again we mark the trial failed
        # otherwise construct the final output
        if (!is(efoAssoc,"json")) {
            V <- data.frame(
                chromosome=efoVariants@variants$chromosome_name,
                position=efoVariants@variants$chromosome_position,
                variant_id=efoVariants@variants$variant_id
            )
            rownames(V) <- V$variant_id
            V <- V[order(V$chromosome,V$position),]
            A <- data.frame(
                variant_id=efoAssoc@risk_alleles$variant_id,
                risk_allele=efoAssoc@risk_alleles$risk_allele,
                risk_frequency=efoAssoc@risk_alleles$risk_frequency
            )
            # Remove duplicated rows (if any) caused by the requested data
            A <- A[!duplicated(A),,drop=FALSE]
            rownames(A) <- A$variant_id
            A <- A[rownames(V),,drop=FALSE]
            
            # Healthy output
            return(list(
                success=TRUE,
                data=cbind(V,A[,c("risk_allele","risk_frequency")]))
            )
        }
        else # Problem in get_associations
            return(list(success=FALSE,data=fromJSON(efoAssoc)))
    }
    else # Problem in get_variants
        return(list(success=FALSE,data=fromJSON(efoVariants)))
}

# validateLoc: in case of duplicates remaining after basic position/id/risk
# allele deduplication, when TRUE, performs location validation, as in some
# cases, duplicate rs have different locations (especially when the genome
# assemble is not reported (nr) and potentially indicating very old genome
# assemblies (<hg18). The process relies on the package biomaRt and uses it to
# query the duplicate rs ids and cross-validate the locations retrieved from
# Ensembl variation databases. SNPs with locations that cannot be validated are
# dropped and deduplicated. As the process involves calling Ensembl APIs, it may
# be slow. Therefore, if the intended use of the retrieved PRS is to plugin the
# SNPs to other GWAS data, validateLoc can be FALSE to speed up the process and
# reduce failures, as the unique SNP ids are used with local genotypes. In other
# cases where strict validation is required (e.g. when ORs must be used or
# reported, then it should be TRUE.
getPGSScores <- function(pgsId=NULL,efoId=NULL,pubmedId=NULL,base=NULL,
    retries=5,validateLoc=FALSE) {
    # Package checking, meaningless to continue
    if (!requireNamespace("quincunx"))
        stop("R package quincunx is required!")
    if (!requireNamespace("jsonlite"))
        stop("R package jsonlite is required!")
    
    # Input checking and sanitization
    if ((missing(pgsId) || is.null(pgsId)) 
        && (missing(efoId) || is.null(efoId))
        && (missing(pubmedId) || is.null(pubmedId)))
        stop("At least one of PGS Id (pgsId) or GWAS EFO ID (efoId) or PubMed ",
            "ID (pubmedId) must be provided!")
    
    # Multiple PGS ids may have been provided... Check them
    if (!is.null(pgsId))
        .checkPGSFormat(pgsId)
    # Same with EFO
    if (!is.null(efoId))
        .checkEFOFormat(efoId)
    # Same with PMID
    if (!is.null(pubmedId))
        .checkPMIDFormat(pubmedId)
        
    # Is base a valid local directory?
    if (!is.null(base) && is.character(base)) {
        if (!dir.exists(base)) {
            warning("The provided local path does not exist! Will use the ",
                "remote path provided by PGS Catalog...",immediate.=TRUE)
            base <- NULL
        }
    }
    
    # Calls
    disp("Scheduling ",length(pgsId) + length(efoId) + length(pubmedId),
        " PGS API calls...")
    
    # PGS id calls
    if (!is.null(pgsId))
        resultPgs <- .iterPGSScoreAPICall(input=pgsId,type="pgs",base=base,
            retries=retries)
    else
        resultPgs <- list(.emptyVariantsDf("pgs"))
    
    # EFO id calls
    if (!is.null(efoId))
        resultEfo <- .iterPGSScoreAPICall(input=efoId,type="efo",base=base,
            retries=retries)
    else
        resultEfo <- list(.emptyVariantsDf("pgs"))
    
    # PMID id calls
    if (!is.null(pubmedId))
        resultPmid <- .iterPGSScoreAPICall(input=pubmedId,type="pmid",base=base,
            retries=retries)
    else
        resultPmid <- list(.emptyVariantsDf("pgs"))
    
    # Join the results
    results <- rbind(do.call("rbind",resultPgs),do.call("rbind",resultEfo),
        do.call("rbind",resultPmid))
    
    # At this point, we must perform a more educated duplicate removal (if any)
    rownames(results) <- NULL
    
    disp("Dealing with duplicates (if any)")
    
    # 0. First remove complete duplicates (chromosome, position, rs if exists,
    # risk_allele). We do not put OR/effect weight, it's from different
    # studies!
    if (any(grepl("^rs",results$variant_id)))
        criteria <- c("chromosome","position","variant_id","risk_allele")
    else
        criteria <- c("chromosome","position","risk_allele")
    if (any(duplicated(results[,criteria]))) {
        disp("  removing exact duplicates (location, risk allele, id ",
            "if available)")
        results <- results[!duplicated(results[,criteria]),,drop=FALSE]
    }
    
    # 1. If rs ids are available duplicates remain in rs id, ensure the proper 
    # SNP locations as there are cases with duplicate rs after joining the 
    # scores and wrong positions. We assume that only if duplicate rs are found,
    # one is mispositioned so we check only these... If other mispositioned 
    # exist, no other choice really. But if all scores are accompanied by rs ids
    # and locations, then no action is taken and we do not know if it's hg38...
    # All these if validateLoc as we have another external API call.
    if (validateLoc && any(grepl("^rs",results$variant_id))) {
        if (any(duplicated(results$variant_id))) {
            results <- results[order(results$variant_id),]
            dup <- which(duplicated(results$variant_id))
            # Not correct as it may catch more rows if an elements is 
            # tuple-icated, but should survive location check
            dup <- unique(unlist(lapply(dup,function(i) { return((i-1):i) })))
            # A function to check only the dupes and return the wrong index so
            # that we can use the initial dupes to remove offending
            #disp("Duplicates found in results. Beginning cleaning process...")
            disp("  validating SNP locations in remaining duplicates")
            rem <- .checkSnpLocs(results[dup,])
            if (length(rem) > 0)
                # This must be removed
                results <- results[-dup[rem],,drop=FALSE]
        }
        
        # 2. If duplicates remain after position check, we check duplicates 
        # for same risk allele. If the same, collapse by averaging OR/weight
        # TODO if needed...
    }
    
    # 3. If rs ids are not available, then only duplicated positions can be 
    # removed, irrespectively of allele order?
    if (!any(grepl("^rs",results$variant_id))) {
        if (any(duplicated(results[,c("chromosome","position",
            "risk_allele")]))) {
            results <- results[-which(duplicated(results[,c("chromosome",   
                "position","risk_allele")])),,drop=FALSE]
        }
    }
    
    #rownames(results) <- results$variant_id
    return(results)
}

.iterPGSScoreAPICall <- function(input,type=c("pgs","efo","pmid"),base=NULL,
    retries=5) {
    type <- type[1]
    N <- length(input)
    result <- vector("list",N)
    names(result) <- input
    
    suff <- ifelse(type=="pgs","PGS ID(s)",ifelse(type=="efo","EFO ID(s)",
        "PubMed ID(s)"))
    disp("  Requesting ",N," PGS Catalog API calls with ",suff,"...")
    
    complete <- FALSE
    times <- 1
    while(!complete && times<=retries) {
        failed <- logical(N)
        names(failed) <- input
        calls <- 0
        
        for (id in input) {
            calls <- calls + 1
            disp("    Call ",calls,": ",id)
            tmp <- .PGSScoreWorker(input=id,type=type,base=base)
            if (tmp$success) # Both API calls successful
                result[[id]] <- tmp$data
            else # Which id failed?
                failed[id] <- TRUE
        }
        
        complete <- !any(failed)
        if (!complete) {
            input <- names(failed)[failed]
            disp("  Failed API calls for ",suff,":")
            disp("    ",paste(input,collapse=", "))
            disp("  Will retry ",retries-times," times...")
        }
        times <- times + 1
    }
    
    return(result)
}

# base is either NULL for default scoring file URL retrieved from quincunx API
# call, or an absolute local path which replaces the
# "http://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores" part of the API retrieved
# path, but expects the same EBI FTP folder structure, that is, if the local
# path is called "/data/resources/PGS", the file structures inside the local
# path should be:
#
# LOCAL_PATH/ 
#   PGSXXXXXX/
#     ScoringFiles/PGSXXXXXX.txt.gz
#     Metadata/PGSXXXXXX_metadata_(extensions)
#   PGSYYYYYY/
#     ScoringFiles/PGSYYYYYY.txt.gz
#     Metadata/PGSYYYYYY_metadata_(extensions)
#
# This can be achieved e.g. by
#
# mkdir -p LOCAL_PATH && cd LOCAL_PATH
# wget -m ftp://anonymous:@ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/
#
# Then base <-  LOCAL_PATH/ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/
# If the lonk path is not desired:
#
# cd LOCAL_PATH/ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/
# mv * ../../../../../../
# cd LOCAL_PATH
# rm -r ./ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/
.PGSScoreWorker <- function(input=NULL,type=c("pgs","efo","pmid"),base=NULL) {
    input <- input[1] # Protection from accidental pass more ids/traits
    type <- type[1]
    
    # Code templates
    arg <- ifelse(type=="pgs","pgs_id",ifelse(type=="efo","efo_id","pubmed_id"))
    exs <- paste0('get_scores(',arg,'="',input,'",warning=FALSE)')
    
    # Get the PGS information
    scoreStru <- tryCatch({
        eval(parse(text=exs))
    },error=function(e) {
        disp("Possible connection failure! Marking...")
        disp("Caught error: ",e$message)
        return(toJSON(list(type=type,input=input),auto_unbox=TRUE))
    })
    
    # Firstly, we need the get_scores call... If failed, we don't proceed to
    # actually fetch the scores from FTP and we mark the whole call as failed
    if (!is(scoreStru,"json")) {
        if (nrow(scoreStru@scores) == 0 || (nrow(scoreStru@scores) == 1 
            && is.na(scoreStru@scores$pgs_id[1]))) # Nothing found - OK
            return(list(success=TRUE,data=.emptyVariantsDf("pgs")))
        else { # Call local score retrieval function
            # Drop non-supported assemblies, e.g. hg18, NCBI36
            sfs <- scoreStru@scores$scoring_file
            pid <- scoreStru@scores$pgs_id
            asm <- scoreStru@scores$assembly
            bad <- .checkSuppAsms(tolower(asm))
            if (!is.null(bad)) {
                warning("Unsupported human genome assemblies found: ",
                    paste0(asm[bad],collapse=", "),".\nThe respective PGS ",
                    "scores (",paste0(pid[bad],collapse=", "),") will be ",
                    "dropped ",immediate.=TRUE)
                sfs <- sfs[-bad]
                pid <- pid[-bad]
                asm <- asm[-bad]
                names(sfs) <- names(asm) <- pid
            }
            
            # Correct sfs for local_path base
            if (!is.null(base)) { # Validation in previous wrapper function
                sfs <- unlist(lapply(sfs,function(x,b) {
                    s <- .splitPath(x)
                    return(file.path(b,s[3],s[2],s[1]))
                },base))
                names(sfs) <- names(asm) <- pid
            }
            
            theScores <- tryCatch({
                disp("Retrieving and enriching scores file(s) for ",input)
                tmpdf <- lapply(sfs,function(x,b) {
                    #Sys.sleep(5)
                    m <- ifelse(is.null(b),"remote","local")
                    disp("  ",m," file ",basename(x))
                    return(.retrieveScoreFile(x))
                },base)
                names(tmpdf) <- pid
                tmpdf
            },error=function(e) {
                disp("Possible connection failure! Marking...")
                disp("Caught error: ",e$message)
                return(toJSON(list(type=type,input=input),auto_unbox=TRUE))
            },finally="")
        }
        
        # Now, if enrichScoreFile call failed, again we mark the trial failed
        # otherwise construct the final output
        if (!is(theScores,"json")) {
            allScores <- lapply(names(theScores),function(n,A,S) {
                disp("  processing ",n)
                gb <- .assemblyToGv(tolower(A[n]))
                pgsScores <- enrichScoreFile(S[[n]],gb,clean=TRUE)
                names(mcols(pgsScores))[names(mcols(pgsScores))=="rsID"] <- 
                    "variant_id"
                names(mcols(pgsScores))[names(mcols(
                    pgsScores))=="effect_allele"] <- "risk_allele"
                    
                # Convert to data frame for compatibility with rest methods and
                # statistical modeling
                pgsScores <- as.data.frame(pgsScores)
                pgsScores <- pgsScores[,names(pgsScores)!="strand"]
                names(pgsScores)[c(1,2)] <- c("chromosome","position")
                pgsScores$asm <- gb
                return(pgsScores)
            },asm,theScores)
            
            pgsScores <- do.call("rbind",allScores)
            
            # Healthy output
            return(list(
                success=TRUE,
                data=pgsScores
            ))
        }
        else # Problem in get_associations
            return(list(success=FALSE,data=fromJSON(theScores)))
    }
    else # Problem in get_scores
        return(list(success=FALSE,data=fromJSON(scoreStru)))
}

.checkSnpLocs <- function(df) {
    if (!requireNamespace("biomaRt"))
        stop("Bioconductor package biomaRt is required!")
    
    # SNP attributes and filters for biomaRt
    filters <- c("snp_filter","variation_source")
    attrs <- c("chr_name","chrom_start","refsnp_id")
    
    # Init offending indices
    offend <- numeric(0)
    
    # Check assemblies
    asms <- unique(df$asm)
    for (a in asms) {
        disp("    assembly: ",a)
        
        subdf <- df[df$asm==a,,drop=FALSE]
        values <- list(snp_filter=unique(subdf$variant_id),
            variation_source="dbSNP")
            
        if (a %in% c("hg19","hg38")) {
            host <- ifelse(a=="hg38","www.ensembl.org","grch37.ensembl.org")
            #mart <- useMart(biomart="ENSEMBL_MART_SNP",host=host,
            #    dataset="hsapiens_snp")
            mart <- useEnsembl(biomart="snps",host=host,dataset="hsapiens_snp")
            snpInfo <- getBM(attributes=attrs,filters=filters,values=values,
                mart=mart)
            tmp <- which(!(subdf$position %in% snpInfo$chrom_start))
            # All wrong?
            if (length(tmp) != nrow(subdf))
                offend <- c(offend,tmp)
        }
        else { # We must check both assemblies...
            #mart1 <- useMart(biomart="ENSEMBL_MART_SNP",host="www.ensembl.org",
            #    dataset="hsapiens_snp")
            mart1 <- useEnsembl(biomart="snps",host="www.ensembl.org",
                dataset="hsapiens_snp")
            snpInfo1 <- getBM(attributes=attrs,filters=filters,values=values,
                mart=mart1)
            tmp1 <- which(!(subdf$position %in% snpInfo1$chrom_start))
            # If got, then there must be some matches
            if (length(tmp1) != nrow(subdf))
                offend <- c(offend,tmp1)
            
            #mart2 <- useMart(biomart="ENSEMBL_MART_SNP",
            #    host="grch37.ensembl.org",dataset="hsapiens_snp")
            mart2 <- useEnsembl(biomart="snps",host="grch37.ensembl.org",
                dataset="hsapiens_snp")
            snpInfo2 <- getBM(attributes=attrs,filters=filters,values=values,
                mart=mart2)
            tmp2 <- which(!(subdf$position %in% snpInfo2$chrom_start))
            # If got, then there must be some matches
            if (length(tmp2) != nrow(subdf))
                offend <- c(offend,tmp2)
        }
    }
    
    return(offend)

    #filters=c("snp_filter","variation_source")
    #values=list(snp_filter="rs6025",variation_source="dbSNP")
}

# Will return only good indexes
.checkSuppAsms <- function(x) {
    ch <- grepl("37",x) | grepl("19",x) | grepl("38",x) | grepl("nr",x)
    if (!all(ch))
        return(which(!ch))
    return(NULL)
}

.assemblyToGv <- function(x=NULL) {
    if (is.null(x))
        return("nr")
    # We expect human genomes
    if (grepl("37",x) || grepl("19",x))
        return("hg19")
    else if (grepl("38",x))
        return("hg37")
    else
        return("nr")
}

.retrieveScoreFile <- function(sf) {
    if (.isValidUrl(sf)) {
        dest <- tempfile()
        download.file(url=sf,destfile=dest,method="libcurl")
    }
    else { # Could be local file
        if (is.character(sf) && file.exists(sf))
            dest <- sf
    }
    return(read.delim(dest,comment.char="#"))
}

# Consider of add hoardr for file caching
# OR = exp(effect_weight)
# effect_weight (beta) = log(OR)
enrichScoreFile <- function(sf,gb=c("hg19","hg38","nr"),clean=FALSE) {
    gb <- gb[1]
    
    # Template fill all the possible columns
    if (is.null(sf$rsID) && is.null(sf$chr_name))
        stop("Both rsID and chromosomal location are missing! Are you sure ",
            "the provided input represents a valid PGS score?")
    if (is.null(sf$effect_allele))
        stop("Score effect allele is missing! Are you sure ",
            "the provided input represents a valid PGS score file?")
    
    if (is.null(sf$rsID)) { # Try and retrieve rs IDs
        disp("    trying to assign dbSNP (rs) ids",level="full")
        rownames(sf) <- paste0(sf$chr_name,"_",sf$chr_position,"_",
            sf$reference_allele,"_",sf$effect_allele)
        sf$id <- rownames(sf) # Assign names anyway
        sf <- .score2GPos(sf)
        if (gb != "nr") { # Else nothing can be done, we keep the defaults
            sf <- .tryAssignSNPrs(sf,from=gb,to="hg38")
            # Do we clean the dataset? Only mapped up-to-date locations
            # and SNP ids?
            if (clean)
                sf <- sf[grep("^rs",sf$id,perl=TRUE)]
        }
        else
            disp("    unable to assign dbSNP (rs) ids: unspecified genome ", 
                "assembly",level="full")
        sf$rsID <- sf$id
    }
    else
        sf$id <- sf$rsID
    
    # If this is valid, then sf is not a GPos yet and has only SNP ids
    if (!is(sf,"GPos") && is.null(sf$chr_name)) {
        # Files are a mess, there may be total duplicate lines
        disp("    removing duplicates and trying to assign genomic position",
            level="full")
        sf <- sf[!duplicated(sf),,drop=FALSE]
        rownames(sf) <- sf$id <- sf$rsID
        # We assign hg38 locations when possible
        sf <- .tryAssignSNPLocs(sf)
        # Do we want a clean dataset? Only up-to-date SNPs?
        if (clean)
            sf <- sf[!(sf$chr_name=="Z" & sf$chr_position==1L),,drop=FALSE]
    }
    
    # Convert to GPos if not yet
    if (!is(sf,"GPos"))
        sf <- .score2GPos(sf)
    
    # Check the rest variables required for a complete object
    if (is.null(sf$locus_name))
        sf <- .tryAssignLocusName(sf,gb)
        #sf$locus_name <- rep(NA,length(sf))
    
    if (is.null(sf$reference_allele)) {
        disp("    trying to infer reference allele",level="full")
        sf <- .tryInferRefAllele(sf,gb)
    }

    if (is.null(sf$effect_weight) && !is.null(sf$OR))
        sf$effect_weight <- log(sf$OR)
    
    if (!is.null(sf$effect_weight) && is.null(sf$OR))
        sf$OR <- exp(sf$effect_weight)
        
    if (is.null(sf$effect_weight) && is.null(sf$OR)) {
        sf$OR <- rep(0,nrow(sf))
        sf$effect_weight <- rep(1,nrow(sf))
    }
    
    # Harmonize the column names of the GPos object and remove the "id" helper
    # column
    mcols(sf) <- mcols(sf)[,c("rsID","effect_allele","reference_allele",
        "locus_name","effect_weight","OR")]
    
    return(unname(sf))
    #return(as(sf,"DataFrame"))
    
    # We also need some - PLINK to GPos thing...
}

# TxDb very limited... sitadela!!! WIP
.tryAssignLocusName <- function(gp,gv=c("hg19","hg38"),sitDb=NULL) {
    gv <- gv[1]
    
    gp$locus_name <- rep(NA,length(gp))
    if (!(gv %in% c("hg19","hg38"))) # Nothing can be done
        return(gp)
    
    if (!requireNamespace("sitadela")) { # Sitadela is required
        warning("Bioconductor package sitadela is required! Returning NAs...")
        return(gp)
    }
    
    # Let's use RefSeq annotation by default
    db <- ifelse(is.null(sitDb),sitadela::getDbPath(),sitDb)
    ann <- loadAnnotation(gv,"refseq",type="gene",version="auto",db=db)
    # We need the same seqlevels styles - for some reason seqlevelsStyles fail
    gptmp <- gp
    gptmp <- tryCatch({
        seqlevelsStyle(gptmp) <- seqlevelsStyle(ann)
    },error=function(e) {
        if (!grepl("chr",as.character(seqnames(gptmp)[1]))) {
            seqlevels(gptmp) <- paste0("chr",seqlevels(gptmp))
            seqnames(gptmp) <- paste0("chr",seqnames(gptmp))
        }
        return(gptmp)
    },finally="")
    
    # Find overlaps
    o <- findOverlaps(gptmp,ann)
    gp$locus_name[queryHits(o)] <- ann$gene_name[subjectHits(o)]
    
    return(gp)
}

.tryInferRefAllele <- function(gp,gv=c("hg19","hg38")) {
    gv <- gv[1]
    
    gp$reference_allele <- rep("N",length(gp))
    if (!(gv %in% c("hg19","hg38"))) # Nothing can be done
        return(gp)
    
    if (!requireNamespace("BSgenome"))
        stop("Bioconductor package BSgenome is required!")
    if (!requireNamespace("SNPlocs.Hsapiens.dbSNP151.GRCh38"))
        stop("Bioconductor package SNPlocs.Hsapiens.dbSNP151.GRCh38 is ",
            "required!")
    bsp <- paste0("BSgenome.Hsapiens.UCSC.",gv)
    if (!requireNamespace(bsp))
        stop("Bioconductor package ",bsp," is required!")
    
    genome <- BSgenome::getBSgenome(gv)
    
    res <- snpsById(
        SNPlocs.Hsapiens.dbSNP151.GRCh38::SNPlocs.Hsapiens.dbSNP151.GRCh38,
        gp$id,ifnotfound="drop")
    seqlevelsStyle(res) <- "UCSC"
    inf <- tryCatch(inferRefAndAltAlleles(res,genome),error=function(e) {
        disp("Error caught during allele inference: ",e$message)
        disp("Will return default (N) as reference allele")
        return(NA)
    },finally="")
    
    if (is.na(inf))
        return(gp)
    
    gp[res$RefSNP_id]$reference_allele <- inf$ref_allele
    
    return(gp)
}

.tryAssignSNPLocs <- function(df) {
    # Accept list of rs and BSgenome
    # Use BSgenome::snpsById
    if (!requireNamespace("BSgenome"))
        stop("Bioconductor package BSgenome is required!")
    if (!requireNamespace("SNPlocs.Hsapiens.dbSNP151.GRCh38"))
        stop("Bioconductor package SNPlocs.Hsapiens.dbSNP151.GRCh38 is ",
            "required!")
    
    res <- snpsById(
        SNPlocs.Hsapiens.dbSNP151.GRCh38::SNPlocs.Hsapiens.dbSNP151.GRCh38,
        df$id,ifnotfound="drop")
    
    # Initialize
    df$chr_name <- rep("Z",nrow(df))
    df$chr_position <- rep(1,nrow(df))
    # Assign found
    df[res$RefSNP_id,"chr_name"] <- as.character(seqnames(res))
    df[res$RefSNP_id,"chr_position"] <- start(res)
    
    return(df)
}

# test: gp[13] old id = "rs11583200", newid = "rs976463154"
# test: is(gp,"GPos")
.tryAssignSNPrs <- function(gp,from,to) {
    # We have locations, we don't have rs ids
    if (!requireNamespace("BSgenome"))
        stop("Bioconductor package BSgenome is required!")
    if (!requireNamespace("SNPlocs.Hsapiens.dbSNP151.GRCh38"))
        stop("Bioconductor package SNPlocs.Hsapiens.dbSNP151.GRCh38 is ",
            "required!")
    
    # Convert gp to hg38 so as to match SNP locations and retrieve SNP ids
    # (if required)
    # Input objects come named
    if (from == "hg19" && to == "hg38")
        gpl <- .liftOverSNPs(gp,from,to)
    else
        gpl <- gp
    
    # Now query dbSNP base on genomic positions
    res <- snpsByOverlaps(
        SNPlocs.Hsapiens.dbSNP151.GRCh38::SNPlocs.Hsapiens.dbSNP151.GRCh38,gpl)
    
    # Assign what is found...
    io <- paste0(seqnames(gpl),"_",start(gpl))
    names(io) <- names(gpl)
    im <- paste0(seqnames(res),"_",start(res))
    mIndex <- match(im,io)
    gp[names(io)[mIndex]]$id <- res$RefSNP_id
    
    return(gp)
}

# test: identical(out,gp) if from=to
# test: !any(identical(starts)) if from!=to
# liftOver not working anymore without list?
.liftOverSNPs <- function(gp,from,to) {
    if (!requireNamespace("rtracklayer"))
        stop("Bioconductor package rtracklayer is required!")

    if (from == to) # Nothing to do
        return(gp)
    
    if (from == "hg38" && to == "hg19")
        chainFile <- system.file(package="liftOver","extdata",
            "hg38ToHg19.over.chain")
    else if (from == "hg19" && to == "hg38") {
        #chainFile <- system.file(package="prisma","extdata",
        #   "hg19ToHg38.over.chain")
        #chgz <- "C:/software/prisma/inst/extdata/hg19ToHg38.over.chain.gz"
        #chgz <- "/media/sevenofnine/raid/software/prisma/inst/extdata/hg19ToHg38.over.chain.gz"
        chgz <- system.file(package="prisma","extdata",
           "hg19ToHg38.over.chain.gz")
        chainFile <- file.path(tempdir(),"hg19ToHg38.over.chain")
        if (!file.exists(chainFile))
            R.utils::gunzip(chgz,destname=chainFile,remove=FALSE)
    }
    ch <- import.chain(chainFile)
    
    orsls <- seqlevelsStyle(gp)
    seqlevelsStyle(gp) <- "UCSC"
    
    sp <- split(gp,seqnames(gp))
    tmp <- lapply(names(sp),function(n,S) {
        disp("  lifting ",n)
        x <- S[[n]]
        lo <- liftOver(x,ch)
        le <- lengths(lo)
        # If more than one mappings (rare) keep 1st
        mult <- which(le > 1)
        if (length(mult) > 0) {
            for (m in mult)
                lo[[m]] <- lo[[m]][1]
        }
        lout <- unlist(lo)
        # If mapping on different chromosome, drop position
        off <- which(seqnames(lout) != n)
        if (length(off) > 0) {
            seqnames(lout)[off] <- n
            start(lout)[off] <- 0
            end(lout)[off] <- 0
            dropsf <- seqlevels(lout)[seqlevels(lout) != n]
            lout <- dropSeqlevels(lout,dropsf)
        }
        return(as.data.frame(lout))
    },sp)
    lifted <- GRanges(do.call("rbind",tmp))
    
    seqlevelsStyle(lifted) <- orsls[1]
    nams <- names(lifted) # Names are not kept in GPos conversion
    out <- GPos(lifted,stitch=FALSE)
    mcols(out) <- mcols(lifted)
    names(out) <- nams
    
    return(out)
}

# test: is(out,GPos)
# test: !is.null(mcols(out))
# test: length(out) == nrow(df)
.score2GPos <- function(df,from=c("pgs","gwas")) {
    from <- from[1]
    # Assuming chr_name, chr_position and position exist
    if (from == "pgs") {
        names(df)[names(df)=="chr_name"] <- "chromosome"
        names(df)[names(df)=="chr_position"] <- "start"
    }
    else if (from == "gwas")
        names(df)[names(df)=="position"] <- "start"
    df$end <- df$start
    interm <- GRanges(df) # Automatically gather mcols
    out <- GPos(interm,stitch=FALSE)
    mcols(out) <- mcols(interm)
    if (!is.null(rownames(df)))
        names(out) <- rownames(df)
    return(out)
}

.emptyVariantsDf <- function(from=c("gwas","pgs")) {
    from <- from[1]
    if (from=="gwas")
        nil <- data.frame(chromosome="1",position=1,variant_id="v",
            risk_allele="N",risk_frequency=9)
    else if (from=="pgs")
        nil <- data.frame(chromosome="1",position=1,variant_id="v",
            risk_allele="N",reference_allele="N",effect_weight=0,
            weight_type="N",variant_description="N",id="N",OR=0)
    return(nil[-1,,drop=FALSE])
}

rsLocsFromEnsembl <- function(rs,gv=c("hg38","hg19"),canonical=TRUE) {
    gv <- gv[1]
    if (!requireNamespace("biomaRt"))
        stop("Bioconductor package biomaRt is required!")
    
    # SNP attributes and filters for biomaRt
    filters <- c("snp_filter","variation_source")
    attrs <- c("chr_name","chrom_start","refsnp_id","allele","minor_allele",
        "minor_allele_freq","chrom_strand")
    values <- list(snp_filter=unique(rs),variation_source="dbSNP")
    
    host <- ifelse(gv=="hg38","www.ensembl.org","grch37.ensembl.org")
    #mart <- useMart(biomart="ENSEMBL_MART_SNP",host=host,dataset="hsapiens_snp")
    mart <- useEnsembl(biomart="snps",host=host,dataset="hsapiens_snp")
    snpInfo <- getBM(attributes=attrs,filters=filters,values=values,mart=mart)
    
    # Keep only canonical chromosomes?
    if (canonical) {
        chrsExp <- paste0("^(",paste0(c(as.character(seq_len(22)),"X","Y"),
            collapse="|"),")$")
        return(snpInfo[grep(chrsExp,snpInfo$chr_name),,drop=FALSE])
    }
    return(snpInfo)
}
