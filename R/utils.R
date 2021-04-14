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
