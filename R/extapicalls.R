getGWASVariants <- function(efoId=NULL,efoTrait=NULL,removeUnknownRisk=TRUE,
    retries=5) {
    # Package checking, meaningless to continue
    if (!require(gwasrapidd))
        stop("R package gwasrapidd is required!")
    if (!require(jsonlite))
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

getPGSVariants <- function(pgsId=NULL,efoId=NULL,pubmedId=NULL,retries=5) {
    # Package checking, meaningless to continue
    if (!require(quincunx))
        stop("R package quincunx is required!")
    if (!require(jsonlite))
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
        .checkEFOFormat(pgsId)
    # Same with PMID
    if (!is.null(pubmedId))
        .checkPMIDFormat(pgsId)
    
    # Calls
    disp("Scheduling ",length(pgsId) + length(efoId) + length(pubmedId),
        " PGS API calls...")
    
    # PGS id calls
    if (!is.null(pgsId))
        resultPgs <- .iterPGSScoreAPICall(input=pgsId,type="pgs",
            retries=retries)
    else
        resultPgs <- .emptyVariantsDf("pgs")
    
    # EFO id calls
    if (!is.null(efoId))
        resultEfo <- .iterPGSScoreAPICall(input=efoId,type="efo",
            retries=retries)
    else
        resultEfo <- .emptyVariantsDf("pgs")
    
    # PMID id calls
    if (!is.null(pubmedId))
        resultPmid <- .iterPGSScoreAPICall(input=pubmedId,type="pmid",
            retries=retries)
    else
        resultPmid <- .emptyVariantsDf("pgs")
    
    # Join the results
    results <- rbind(do.call("rbind",resultPgs),do.call("rbind",resultEfo),
        do.call("rbind",resultPmid))
    
    # Remove duplicates if any - duplication is relevant... We provide some
    # potential combination methods of collapsing ORs and weights if the risk
    # allele is the same, otherwise they are not identical
    results <- results[!duplicated(results),,drop=FALSE]
    rownames(results) <- results$variant_id
        
    return(results)
}

.iterPGSScoreAPICall <- function(input,type=c("pgs","efo","pmid"),retries=5) {
    type <- type[1]
    N <- length(input)
    result <- vector("list",N)
    names(result) <- input
    
    suff <- ifelse(type=="pgs"," PGS ID(s)",ifelse(type=="efo"," EFO ID(s)",
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
            tmp <- .PGSScoreWorker(input=id,type=type)
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

.PGSScoreWorker <- function(input=NULL,type=c("pgs","efo","pmid")) {
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
        if (nrow(scoreStru@scores) == 0) # Nothing found - OK
            return(list(success=TRUE,data=.emptyVariantsDf("pgs")))
        else # Call local score retrieval function
            theScores <- tryCatch({
                disp("Retrieving and enriching scores file for ",input)
                .retrieveScoreFile(scoreStru@scores$scoring_file)
            },error=function(e) {
                disp("Possible connection failure! Marking...")
                disp("Caught error: ",e$message)
                return(toJSON(list(type=type,input=input),auto_unbox=TRUE))
            },finally="")
        
        # Now, if enrichScoreFile call failed, again we mark the trial failed
        # otherwise construct the final output
        if (!is(theScores,"json")) {
            gb <- .assemblyToGv(scoreStru@scores$assembly)
            pgsScores <- enrichScoreFile(theScores,gb,clean=TRUE)
            names(mcols(pgsScores))[names(mcols(pgsScores))=="rsID"] <- 
                "variant_id"
            names(mcols(pgsScores))[names(mcols(pgsScores))=="effect_allele"] <-
                "risk_allele"
            
            # Finally, convert to data frame for compatibility with rest methods
            # and statistical modeling
            pgsScores <- as.data.frame(pgsScores)
            pgsScores <- pgsScores[,names(pgsScores)!="strand"]
            names(pgsScores)[c(1,2)] <- c("chromosome","position")
            
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
        download.file(url=sf,destfile=dest)
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
            basename(scoreFile)," is a valid PGS score file?")
    if (is.null(sf$effect_allele))
        stop("Score effect allele is missing! Are you sure ",
            basename(scoreFile)," is a valid PGS score file?")
    
    if (is.null(sf$rsID)) { # Try and retrieve rs IDs
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
    }
    else
        sf$id <- sf$rsID
    
    # If this is valid, then sf is not a GPos yet and has only SNP ids
    if (!is(sf,"GPos") && is.null(sf$chr_name)) {
        # Files are a mess, there may be total duplicate lines
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
        # TODO: Finish .tryAssignLocusName
        # sf <- .tryAssignLocusName(sf,gb)
        sf$locus_name <- rep(NA,length(sf))
    
    if (is.null(sf$reference_allele))
        sf <- .tryInferRefAllele(sf,gb)

    if (is.null(sf$effect_weight) && !is.null(sf$OR))
        sf$effect_weight <- log(sf$OR)
    
    if (!is.null(sf$effect_weight) && is.null(sf$OR))
        sf$OR <- exp(sf$effect_weight)
        
    if (is.null(sf$effect_weight) && is.null(sf$OR)) {
        sf$OR <- rep(0,nrow(sf))
        sf$effect_weight <- rep(1,nrow(sf))
    }
    
    return(sf)
    #return(as(sf,"DataFrame"))
    
    # We also need some - PLINK to GPos thing...
}

# TxDb very limited... sitadela!!! WIP
.tryAssignLocusName <- function(gp,gv=c("hg19","hg38")) {
    gv <- gv[1]
    
    gp$locus_name <- rep(NA,length(gp))
    if (!(gv %in% c("hg19","hg38"))) # Nothing can be done
        return(gp)
    
    if (!requireNamespace("sitadela"))
        stop("Bioconductor package sitadela is required!")
    
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
    
    res <- snpsById(SNPlocs.Hsapiens.dbSNP151.GRCh38,gp$id,ifnotfound="drop")
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
    
    res <- snpsById(SNPlocs.Hsapiens.dbSNP151.GRCh38,df$id,ifnotfound="drop")
    
    # Initialize
    df$chr_name <- rep("Z",nrow(tmp))
    df$chr_position <- rep(1,nrow(tmp))
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
    res <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP151.GRCh38,gpl)
    
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
        chgz <- "/media/raid/software/prisma/inst/extdata/hg19ToHg38.over.chain.gz"
        chainFile <- file.path(tempdir(),"hg19ToHg38.over.chain")
        if (!file.exists(chainFile))
            R.utils::gunzip(chgz,destname=chainFile,remove=FALSE)
    }
    ch <- import.chain(chainFile)
    
    orsls <- seqlevelsStyle(gp)
    seqlevelsStyle(gp) <- "UCSC"
    lifted <- liftOver(gp,ch)
    lifted <- unlist(lifted)
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

#~ .prepareScores <- function(scores,gb=c("hg19","hg38")) {
#~     # Accept DataFrame object
#~     # Depending on its contents
#~     # .retrieveSNPLocation
#~     # .retrieveSNPrs
#~     # Add some metacolumn to indicate inverse allele (some thinking)
#~     # Return a GPos object with as much metadata as possible from GWAS or PGS
    
#~     if (!is(scores,"DataFrame"))
#~         stop("Input must be a DataFrame object!")
    
#~     # Check if we need to retrieve locations
#~     if (scores$tmp$chr_name)
#~ }
