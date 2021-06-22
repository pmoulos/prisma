partitionGWAS <- function(obj,by,n=1,frac=0.5,ngrp=min(5,ncol(obj)),
    replace=FALSE,out=c("index","train","ttboth"),rc=NULL) {
    if (!is(obj,"GWASExperiment"))
        stop("obj must be a GWASExperiment object!")
    
    p <- phenotypes(obj)
    if (is.null(p) || nrow(p) == 0)
        stop("obj does not have associated phenotypes!")
    
    if (missing(by)) {
        warning("by argument is missing! Will use the first variable ",
            "from object's phenotypes...",immediate.=TRUE)
        by <- names(p)[1]
    }
    if (!(by %in% names(p)))
        stop("Unkonwn phenotype ",by,"! The provided phenotypes are: ",
            paste(names(p),collapse=", "))
    
    out <- out[1]
    .checkTextArgs("out",out,c("index","train","ttboth"),multiarg=FALSE)
    .checkNumArgs("n",as.integer(n),"integer",1L,"gte")
    .checkNumArgs("frac",frac,"numeric",c(0,1),"both")
    .checkNumArgs("ngrp",as.integer(ngrp),"integer",1L,"gte")
    if (!is.logical(replace))
        stop("The replace argument must be TRUE or FALSE")
    
    disp("Partitioning dataset into ",n," partitions by ",by)
    y <- p[,by]
    split <- createSplit(y,n=n,frac=frac,ngrp=ngrp,replace=replace,rc=rc)
    if (n == 1) {
        if (out == "train")
            return(obj[,split[[1]],drop=FALSE])
        else if (out == "ttboth")
            return(list(
                train=obj[,split[[1]],drop=FALSE],
                test=obj[,-split[[1]],drop=FALSE]
            ))
        else
            return(split)
    }
    else
        return(split)
}

createSplit <- function (y,n=1,frac=0.5,ngrp=min(5,length(y)),replace=FALSE,
    out=c("index","binary"),rc=NULL) {
    out <- out[1]
    .checkTextArgs("out",out,c("index","binary"),multiarg=FALSE)
    .checkNumArgs("n",as.integer(n),"integer",1L,"gte")
    .checkNumArgs("frac",frac,"numeric",c(0,1),"both")
    .checkNumArgs("ngrp",as.integer(ngrp),"integer",1L,"gte")
    if (!is.logical(replace))
        stop("The replace argument must be TRUE or FALSE")
    
    if (length(y) < 2) 
        stop("The class to be split must have at least 2 data points!")
    if (ngrp < 2) 
        ngrp <- 2
    if (is.numeric(y) || is.integer(y)) {
        y <- cut(y,unique(quantile(y,probs=seq(0,1,length=ngrp))), 
            include.lowest=TRUE)
    }
    else {
        xtab <- table(y)
        if (any(xtab == 0)) {
            warning(paste("Some classes have no records (",
                paste(names(xtab)[xtab==0],sep="",collapse= ", "),") and ",
                "these will be ignored"),immediate.=TRUE)
            y <- factor(as.character(y))
        }
        if (any(xtab == 1)) {
            warning(paste("Some classes have a single record (", 
                paste(names(xtab)[xtab==1],sep="",collapse=", "), ") and ",
                "these will be selected for the sample"),immediate.=TRUE)
        }
    }
    
    N <- seq_len(n)
    splits <- cmclapply(N,function(i,x,p,r) {
        tmp <- data.frame(x=x,index=seq_along(x))
        if (nrow(tmp) == 1)
            return(tmp$index)
        else
            return(sort(sample(tmp$index,size=ceiling(nrow(tmp)*p),replace=r)))
    },y,frac,replace,rc=rc)
    
    if (out=="binary") {
        splits <- cmclapply(splits,function(x,n) {
            z <- integer(n)
            z[x] <- 1
            return(z)
        },length(y))
    }
    
    names(splits) <- as.character(N)
    
    return(splits)
}

createSample <- function() {
    
}

getDefaults <- function(what) {
    allowed <- c("filters","gwargs")
    
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
                nPC=NA,
                LD=0.2,
                IBD=0.1,
                inbreed=0.1
            ))
        },
        gwargs = {
            return(list(
                glm=list(
                    family="gaussian"
                ),
                rrblup=list(
                    pcblup="auto",
                    npcs=NULL
                ),
                statgen=list(
                    # Nothing for now
                ),
                snptest=list(
                    test="frequentist",
                    model="additive"
                )
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

prismaVerbosity <- function(level=NULL) {
    if (missing(level) || is.null(level)) {
        v <- getOption("prisma_verbosity")
        if (!is.null(v))
            return(v)
        else
            level  <- "normal"
    }
    .checkTextArgs("verbosity level",level,c("silent","normal","full"),
        multiarg=FALSE)
    options("prisma_verbosity"=level)
}

disp <- function(...,level=c("normal","full")) {
    level <- level[1]
    v <- prismaVerbosity()
    if (v=="full" && level=="normal" || v=="full" && level=="full"
        || v=="normal" && level=="normal")
        message(...)
}

.splitFactorForParallel <- function(n,rc) {
    .checkNumArgs("The number of elements to split",n,"numeric",0,"gt")
    nc <- .coresFrac(rc)
    mo <- n%%nc
    if (mo == 0) {
        fl <- n/nc
        return(factor(rep(seq_len(nc),each=fl)))
    }
    else {
        fl <- (n-mo)/nc
        return(factor(c(rep(seq_len(nc),each=fl),rep(nc,mo))))
    }
}

.coresFrac <- function(rc=NULL) {
    if (missing(rc) || is.null(rc))
        return(1)
    else {
        n <- parallel::detectCores()
        if (n > 1)
            return(ceiling(rc*n))
        else 
            return(1)
    }
}

.toIntMat <- function(x) {
    x <- round(x,digits=1)
    for (i in seq_len(ncol(x)))
        x[,i] <- as.integer(x[,i])
    return(x)
}

.validateBinaryForBinomial <- function(x) {
    if (length(unique(x)) > 2)
        stop("The response variable cannot have more than two values when ",
            "GLM family is 'binomial'!")
    if (!is.factor(x)) {
        if (!all(unique(x) %in% c(0,1))) {
            warning("When GLM family is 'binomial', the response variable ",
                "must be either 0, 1 or a 2-level factor! Converting to ",
                "factor...")
            x <- as.factor(x)
        }
    }
    return(x)
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
    if (argType == "numeric" && check$cls == "integer")
        argType <- "integer"
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

.downloadSnptest <- function(ver=c("v2.5.6","v2.5.4","v2.5.2")) {
    # This downloads the latest version by default... Might not work for older
    # Linux operating systems
    ver <- ver[1]
    base <- "http://www.well.ox.ac.uk/~gav/resources/"
    if (ver == "v2.5.6")
        ext = "snptest_v2.5.6_CentOS_Linux7.8-x86_64_dynamic.tgz"
    else if (ver == "v2.5.4")
        ext <- "snptest_v2.5.4-beta3_linux_x86_64_static.tgz"
    else if (ver == "v2.5.2")
        ext <- "snptest_v2.5.2_linux_x86_64_static.tgz"
    src <- paste0(base,ext)
    des <- file.path(system.file(package="prisma"),"tools",ext)
    download.file(src,des)
    untar(des,exdir=file.path(system.file(package="prisma"),"tools"))
    unlink(des)
}

.randomString <- function(n=1) {
  a <- do.call(paste0,replicate(5,sample(LETTERS,n,replace=TRUE),FALSE))
  return(paste0(a,sprintf("%04d",sample(9999,n,replace=TRUE)),sample(LETTERS,n,
    replace=TRUE)))
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
