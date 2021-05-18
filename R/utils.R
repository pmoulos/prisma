getDefaults <- function(what) {
    allowed <- c("filters")
    
    if (!(what %in% allowed))
        stop("what must be one of ",paste(allowed,collapse=", "))
    
    switch(what,
        filters = {
            return(list(
                snpCallRate=0.98,
                sampleCallRate=0.95,
                maf=0.05,
                hwe=1e-6,
                heteroStat="median",
                heteroFac=3,
                heteroHard=NA,
                pcaOut=TRUE,
                pcaRobust="grid",
                nPC=0,
                LD=0.2,
                IBD=0.1,
                inbreed=0.1
            ))
        }
    )
}

cmclapply <- function(...,rc) {
    if (suppressWarnings(!requireNamespace("parallel")) 
        || .Platform$OS.type!="unix")
        ncores <- 1
    else
        ncores <- .coresFrac(rc)
    if (ncores > 1)
        return(mclapply(...,mc.cores=ncores,mc.set.seed=FALSE))
    else
        return(lapply(...))
}

.coresFrac <- function(rc=NULL) {
    if (missing(rc) || is.null(rc))
        return(1)
    else {
        n <- parallel::detectCores()
        if (n > 1)
            return(ceiling(rc*ncores))
        else 
            return(1)
    }
}

.checkEFOFormat <- function(id) {
    nInput <- length(id)
    notPrefixCompliant <- notUnderscoreCompliant <- notSuffixCompliant <- NULL
    
    # Check prefix of EFO id
    checkPrefix <- grepl("^EFO",id,perl=TRUE)
    if (!all(checkPrefix)) {
        notPrefixCompliant <- paste(id[!checkPrefix],collapse=", ")
        id <- id[checkPrefix]
    }
    
    if (length(id) == 0)
        stop("All the provided EFO ids do not seem to be EFO ids: ",
            notPrefixCompliant)
    
    # Check if underscore exists
    test <- strsplit(id,"_")
    checkUnderscore <- lengths(test) == rep(2,length(id))
    if (!all(checkUnderscore)) {
        notUnderscoreCompliant <- paste(id[!checkUnderscore],collapse=", ")
        id <- id[checkUnderscore]
    }
    
    if (length(id) == 0)
        stop("All the provided EFO ids do not seem to be EFO ids: ",
            paste0(c(notPrefixCompliant,", ",notUnderscoreCompliant)))
    
    # After underscore test, accessing the second EFO id part should not fail
    checkSuffix <- unlist(lapply(test,function(x) {
        !is.na(suppressWarnings(as.numeric(x[2])))
    }))
    if (!all(checkSuffix)) {
        notSuffixCompliant <- paste(id[!checkSuffix],collapse=", ")
        id <- id[checkSuffix]
    }
    
    if (length(id) == 0)
        stop("All the provided EFO ids do not seem to be EFO ids: ",
            paste0(c(notPrefixCompliant,", ",notUnderscoreCompliant,", ",
            notSuffixCompliant)))
    
    # Warning message if problems found and ids remaining
    if (length(id) != nInput)
        warning("The following provided EFO ids do not seem to be EFO ids: ",
            paste0(c(notPrefixCompliant,", ",notUnderscoreCompliant,", ",
            notSuffixCompliant)),immediate.=TRUE)
    
    return(id)
}

.checkPGSFormat <- function(id) {
    nInput <- length(id)
    notPrefixCompliant <- notSuffixCompliant <- NULL
    
    # Check prefix of EFO id
    checkPrefix <- grepl("^PGS",id,perl=TRUE)
    if (!all(checkPrefix)) {
        notPrefixCompliant <- paste(id[!checkPrefix],collapse=", ")
        id <- id[checkPrefix]
    }
    
    if (length(id) == 0)
        stop("All the provided PGS ids do not seem to be PGS ids: ",
            notPrefixCompliant)
    
    # Second PGS id part should be integer of length 6
    test <- strsplit(id,"S")
    checkSuffix <- unlist(lapply(test,function(x) {
        !is.na(suppressWarnings(as.integer(x[2]))) && nchar(x[2])==6
    }))
    if (!all(checkSuffix))
        notSuffixCompliant <- paste(id[!checkSuffix],collapse=", ")
        id <- id[checkSuffix]
    
    if (length(id) == 0)
        stop("All the provided PGS ids do not seem to be PGS ids: ",
            paste0(c(notPrefixCompliant,", ",notSuffixCompliant)))
    
    # Warning message if problems found and ids remaining
    if (length(id) != nInput) {
        warning("The following provided PGS ids do not seem to be PGS ids: ",
            paste0(c(notPrefixCompliant,", ",notSuffixCompliant)),
            immediate.=TRUE)
    }
    
    return(id)
}

.checkPMIDFormat <- function(id) {
    nInput <- length(id)
    notCompliant <- NULL
     
    # Check PMID is integer of length up to 8
    check <- !(is.na(suppressWarnings(as.integer(id))) | nchar(id) > 8)
    if (!all(check)) {
        notCompliant <- paste(id[!check],collapse=", ")
        id <- id[check]
    }
    
    if (length(id) == 0)
        stop("All the provided PubMed ids do not seem to be PubMed ids: ",
            notCompliant)
    
    # Warning message if problems found and ids remaining
    if (length(id) != nInput)
        warning("The following provided PubMed ids do not seem to be PubMed ",
            "ids: ",notCompliant,immediate.=TRUE)
    
    return(id)
}

.emptyVariantsDf <- function() {
    nil <- data.frame(chromosome="1",position=1,variant_id="v",risk_allele="N",
        risk_frequency=9)
    return(nil[-1,,drop=FALSE])
}

.isValidUrl <- function(url) {
    regx <- paste0("((([A-Za-z]{3,9}:(?:\\/\\/)?)(?:[-;:&=\\+\\$,\\w]+@)",
        "?[A-Za-z0-9.-]+|(?:www.|[-;:&=\\+\\$,\\w]+@)[A-Za-z0-9.-]+)",
        "((?:\\/[\\+~%\\/.\\w\\-_]*)?\\??(?:[-\\+=&;%@.\\w_]*)#?(?:[\\w]*))?)")
    return(is.character(url) && grepl(regx,url,perl=TRUE))
}

.checkNumArgs <- function(argName,argValue,argType,argBounds,direction) {
    # First generic check so not to continue if fail
    if (!is(argValue,argType))
        stop("\"",argName,"\" parameter must be a(n) ",argType," value!")
    
    # Then, proceed with a lookup table to avoid repetition (suggested by
    # Marcel Ramos during package review)
    lookup <- list(
        both=list(
            fail=function(x) x<argBounds[1] || x>argBounds[2],
            cls=class,
            msg=function(x) paste0("larger than or equal to ",
                argBounds[1]," and smaller than or equal to ",
                argBounds[2])
        ),
        botheq=list(
            fail=function(x) x<=argBounds[1] || x>=argBounds[2],
            cls=class,
            msg=function(x) paste0("larger than ",argBounds[1],
                " and smaller than ",argBounds[2])
        ),
        gt=list(
            fail=function(x) x<=argBounds[1], 
            cls=class,
            msg=function(x) paste0("greater than ",argBounds[1])
        ),
        lt=list(
            fail=function(x) x>=argBounds[1], 
            cls=class,
            msg=function(x) paste0("lower than ",argBounds[1])
        ),
        gte=list(
            fail=function(x) x<argBounds[1], 
            cls=class,
            msg=function(x) paste0("greater than or equal to ",
                argBounds[1])
        ),
        lte=list(
            fail=function(x) x>argBounds[1], 
            cls=class,
            msg=function(x) paste0("lower than or equal to ",
                argBounds[1])
        )
    )
    
    check <- lapply(lookup[[direction]],function(f) f(argValue))
    if (check$fail || check$cls != argType)
        stop("\"",argName,"\""," parameter must be a(n) ",argType,
            " value ",check$msg,"!")
}

.checkTextArgs <- function(argName,argValue,argList,multiarg=FALSE) {
    if (multiarg) {
        argValue <- tolower(argValue)
        if (!all(argValue %in% argList))
            stop("\"",argName,"\""," parameter must be one or more of ",
                paste(paste("\"",argList,sep=""),collapse="\", "),"\"!")
    }
    else {
        argSave <- argValue[1]
        argValue <- tolower(argValue[1])    
        if (!(argValue %in% argList))
            stop("\"",argName,"\""," parameter must be one of ",
                paste(paste("\"",argList,sep=""),collapse="\", "),"\"!")
    }
}

.getArg <- function(argList,argName) {
    return(argList[argName])
}

.setArg <- function(argList,argName,argValue=NULL) {
    if (is.list(argName))
        argList[names(argName)] <- argName
    else if (is.character(argName)) {
        tmp <- vector("list",length(argName))
        names(tmp) <- argName
        i <- 0
        for (n in argName) {
            i <- i + 1
            tmp[[n]] <- argValue[i]
        }
        argList[argName] <- tmp
    }
    return(argList)
}

.gimmeTestFiles <- function() {
    fam <- system.file("extdata/sample.fam",package="snpStats")
    bim <- system.file("extdata/sample.bim",package="snpStats")
    bed <- system.file("extdata/sample.bed",package="snpStats")
    return(list(fam=fam,bim=bim,bed=bed))
}

.isEmpty <- function(x) {
    return(is.null(x) || is.na(x) || x == "" || length(x) == 0)
}


#~ getAPIBase <- function() {
#~     base <- getOption("rpgscat_base")
#~     if (!is.null(base) && .isValidUrl(base))
#~         return(base)
#~     else
#~         return(.defaultUrlBase())
#~ }

#~ setAPIBase <- function(url=NULL) {
#~     if (is.null(url))
#~         url <- .defaultUrlBase()
#~     else {
#~         if (!.isValidUrl(url)) {
#~             warning("The provided URL string does not seem to be a valid ",
#~                 "URL! Assuming default PGS Catalog base...",immediate.=TRUE)
#~             url <- .defaultUrlBase()
#~         }
#~     }
#~     options("rpgscat_base"=url)
#~ }

#~ .defaultUrlBase <- function() {
#~     return("https://www.pgscatalog.org/rest")
#~ }
