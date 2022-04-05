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

cmclapply <- function(...,rc,setseed=FALSE,preschedule=TRUE) {
    if (suppressWarnings(!requireNamespace("parallel")) 
        || .Platform$OS.type!="unix")
        ncores <- 1
    else
        ncores <- .coresFrac(rc)
    if (ncores > 1)
        return(mclapply(...,mc.cores=ncores,mc.set.seed=setseed,
            mc.preschedule=preschedule))
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
    .checkTextArgs("verbosity level",level,c("silent","normal","full","debug"),
        multiarg=FALSE)
    options("prisma_verbosity"=level)
}

disp <- function(...,level=c("normal","full","debug")) {
    level <- level[1]
    v <- prismaVerbosity()
    if ((v=="full" && level=="normal") || (v=="full" && level=="full")
        || (v=="normal" && level=="normal") || (v=="debug" && level=="normal")
        || (v=="debug" && level=="full") || (v=="debug" && level=="debug"))
        message(...)
}

.checkAndCorrectFactorsFormat <- function(x) {
    # First two are FID+IID - PRSice does not care
    cls <- unlist(lapply(x,class))[3:ncol(x)]
    if (!any(cls=="character")) # All OK
        return(x)
    
    crc <- which(cls=="character") + 2 # FID + IID removed!
    for (j in crc) {
        tmp <- x[,j]
        tmp <- as.factor(tmp)
        levels(tmp) <- seq_along(levels(tmp))
        x[,j] <- tmp
    }
    message("    variable(s) ",paste(names(x)[crc],sep=", ")," converted ",
        "to factor")
    
    return(x)
}

.splitFactorForParallel <- function(n,rc) {
    .checkNumArgs("The number of elements to split",n,"numeric",0,"gt")
    if (is.null(rc))
        return(factor(rep(1,n)))
    if (rc>0 && rc<=1) # Cores fraction,else chunks
        nc <- .coresFrac(rc)
    else
        nc <- rc
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
    if (!is.character(argValue))
        stop(argValue," must be a character scalar or vector!")
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

.validateWorkspacePath <- function(path,tool="tool") {
    # Various path options
    if (!is.null(path)) {
        if (!is.character(path)) {
            warning("The desired ",tool," workspace must be a valid string! ",
                "Setting to NULL for autocreation",immediate.=TRUE)
            path <- NULL
        }
        else {
            if (!dir.exists(path)) {
                #if (!interactive()) # Just create...
                    dir.create(path,recursive=TRUE)
                #else {
                #    # No disp here
                #    message("The desired ",tool," workspace does not exist!")
                #    confirm <- menu(c("Yes","No"),
                #        title=sprintf(paste0(tool," workspace is going to be ",
                #            "created at %s. Do you agree?"),path))
                #    if (confirm == 1)
                #        dir.create(path,recursive=TRUE)
                #    else
                #        stop(paste0(tool," cannot be executed without a valid ",
                #            "workspace.\nPlease agree or use the 'path' ",
                #            "argument in the respective function."))
                #}
            }
        }
    }
    else {
        path <- file.path(tempdir(),paste0(tool,"_",.randomString()))
        dir.create(path,recursive=TRUE,showWarnings=FALSE)
    }
    
    # Finally
    disp(tool," operations will be executed at ",path)
    return(path)
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

.randomString <- function(n=1,s=5) {
  a <- do.call(paste0,replicate(s,sample(LETTERS,n,replace=TRUE),FALSE))
  return(paste0(a,sprintf("%04d",sample(9999,n,replace=TRUE)),sample(LETTERS,n,
    replace=TRUE)))
}

.capFirst <- function(x) {
  z <- strsplit(x," ")[[1]]
  return(paste(toupper(substring(z,1,1)),substring(z,2),sep="",collapse=" "))
}

.splitPath <- function(x) {
    if (dirname(x) == x) 
        return(x) 
    else 
        return(c(basename(x),.splitPath(dirname(x))))
}

.whatIsMyName <- function(x) {
    return(deparse(substitute(x)))
}

.formatSystemOutputForDisp <- function(x) {
    clean1 <- gsub("^\\s*\\[\\d+\\]\\s+","",x,perl=TRUE)
    clean2 <- gsub('\\"','',clean1)
    return(paste(clean2,collapse="\n"))
}

.elap2human <- function(start.time) {
    start.time <- as.POSIXct(start.time)
    dt <- difftime(Sys.time(),start.time,units="secs")
    ndt <- as.numeric(dt)
    if (ndt<60)
        format(.POSIXct(dt,tz="GMT"),"%S seconds")
    else if (ndt>=60 && ndt<3600)
        format(.POSIXct(dt,tz="GMT"),"%M minutes %S seconds")
    else if (ndt>=3600 && ndt<86400)
        format(.POSIXct(dt,tz="GMT"),"%H hours %M minutes %S seconds")
    else if (ndt>=86400)
        format(.POSIXct(dt,tz="GMT"),"%d days %H hours %M minutes %S seconds")
}

.dispListOpts <- function(x,lo=1) {
    nm <- nchar(names(x))
    ns <- max(nm) + 1 - nm
    names(ns) <- names(x)
    ss <- lapply(ns,function(z) 
        paste(paste(rep(" ",z),collapse=""),": ",sep=""))
    
    if (lo > 1)
        off <- paste(rep("  ",lo),sep="")
    else
        off <- "  "
        
    for (n in names(x)) {
        y <- x[[n]]
        if (is.logical(y))
            disp(off,paste(n,ifelse(y,"Yes","No"),sep=ss[[n]]))
        else if (.isEmpty(y))
            disp(off,paste(n,"-",sep=ss[[n]]))
        else
            disp(off,paste(n,y,sep=ss[[n]]))
    }
}

.symbolBar <- function(sym,len) {
    return(paste(rep(sym,len),collapse=""))
}

.makeReportMessages <- function(lang) {
    switch(lang,
        en = {
            messages <- list(
                stat=list(
                    glm="GLM",
                    snptest="SNPTEST",
                    rrblup="rrBLUP",
                    statgen="statgenGWAS",
                    plink="PLINK"
                ),
                prs=list(
                    lassosum="lassosum",
                    prsice="PRSice-2"
                ),
                meta=list(
                    intersection="intersection of individual results",
                    union="union of individual results",
                    fisher="Fisher's method",
                    minp="minimum p-value across results",
                    maxp="maximum p-value across results",
                    pandora="PANDORA weighted p-value across results",
                    simes="Simes correction and combination method",
                    whitlock=paste("Whitlock's Z-transformation method",
                        "(Bioconductor package survcomp)"),
                    harmonic="Wilson's Harmonic mean of p-values",
                    none=paste("no meta-analysis, reported p-values from the",
                        "first supplied statistical algorithm")
                ),
                adjust=list(
                    holm="Holm FWER",
                    hochberg="Hochberg DFR",
                    hommel="Hommel FWER",
                    bonferroni="Bonferroni FWER",
                    bh="Benjamini-Hochberg FDR",
                    by="Benjamini-Yekutiely FDR",
                    fdr="Benjamini-Hochberg FDR",
                    none="no multiple test correction",
                    qvalue="Storey-Tibshirani FDR"
                ),
                plots=list(
                    mds="multidimensional scaling"
                ),
                export=list(
                    annotation="Annotation",
                    p_value="p-value"
                ),
                explain=list(
                    mds=paste(
                "Multidimensional Scaling (MDS) plots constitute a means",
                "of visualizing the level of similarity of individual cases",
                "of a dataset. It is similar to Principal Component Analysis",
                "(PCA), but instead of using the covariance matrix to find",
                "similarities between cases, MDS uses absolute distance",
                "metrics such as the classical Euclidean distance. Because",
                "of the relative linear relations between sequencing samples,",
                "it provides a more realistic clustering of samples. MDS",
                "serves quality control and it can be interpreted as follows:",
                "when the distance between samples of the same biological",
                "condition in the MDS space is small, this is an indication",
                "of high correlation and reproducibility between them. When",
                "this distance is larger or heterogeneous (e.g. the 3rd",
                "sample of a triplicate set is further from the other 2),",
                "this constitutes an indication of low correlation and",
                "reproducibility between samples. It can help exclude poor",
                "samples from further analysis.",collapse=" "
                    )
                ),
                references=list(
                    qc=paste("Moulos, P., Hatzis, P. (2015). Systematic",
                        "integration of RNA-Seq statistical algorithms for",
                        "accurate detection of differential gene expression",
                        "patterns. Nucleic Acids Research 43(4), e25."),
                    prs=list(
                        edaseq=paste("Risso, D., Schwartz, K., Sherlock, G.,",
                            "and Dudoit, S. (2011). GC-content normalization",
                            "for RNA-Seq data. BMC Bioinformatics 12, 480.")
                    ),
                    stat=list(
                        deseq=paste("Anders, S., and Huber, W. (2010).",
                            "Differential expression analysis for sequence",
                            "count data. Genome Biol 11, R106.")
                    ),
                    meta=list(
                        fisher=paste("Fisher, R.A. (1932). Statistical",
                            "Methods for Research Workers (Edinburgh, Oliver",
                            "and Boyd)."),
                        fperm=paste("Fisher, R.A. (1932). Statistical",
                            "Methods for Research Workers (Edinburgh, Oliver",
                            "and Boyd)."),
                        whitlock=c(
                            paste("Whitlock, M.C. (2005). Combining",
                                "probability from independent tests:",
                                "the weighted Z-method is superior to Fisher's",
                                "approach. J Evol Biol 18, 1368-1373."),
                            paste("Schroder, M.S., Culhane, A.C., Quackenbush,",
                                "J., and Haibe-Kains, B. (2011). survcomp:",
                                "an R/Bioconductor package for performance",
                                "assessment and comparison of survival",
                                "models. Bioinformatics 27, 3206-3208.")
                        ),
                        weight=paste("Genovese, C.R., Roeder, K., Wasserman,",
                            "L. (2006). False discovery control with p-value",
                            "weighting. Biometrika 93 (3): 509-524."),
                        simes=paste("Simes, R. J. (1986). An improved",
                            "Bonferroni procedure for multiple tests of",
                            "significance. Biometrika 73 (3): 751-754."),
                        harmonic=paste("Wilson, D. J. (2019). The harmonic",
                            "mean p-value for combining dependent tests.",
                            "PNAS 116 (4): 1995-1200."),
                        none=NULL
                    ),
                    multiple=list(
                        BH=paste("Benjamini, Y., and Hochberg, Y. (1995). ",
                            "Controlling the False Discovery Rate: A Practical",
                            "and Powerful Approach to Multiple Testing.",
                            "Journal of the Royal Statistical Society Series",
                            "B (Methodological) 57, 289-300."),
                        fdr=paste("Benjamini, Y., and Hochberg, Y. (1995). ",
                            "Controlling the False Discovery Rate: A Practical",
                            "and Powerful Approach to Multiple Testing.",
                            "Journal of the Royal Statistical Society Series",
                            "B (Methodological) 57, 289-300."),
                        BY=paste("Benjamini, Y., and Yekutieli, D. (2001). The",
                            "control of the false discovery rate in multiple",
                            "testing under dependency. Annals of Statistics",
                            "26, 1165-1188."),
                        bonferroni=paste("Shaffer, J.P. (1995). Multiple",
                            "hypothesis testing. Annual Review of",
                            "Psychology 46, 561-576."),
                        holm=paste("Holm, S. (1979). A simple sequentially",
                            "rejective multiple test procedure. Scandinavian",
                            "Journal of Statistics 6, 65-70."),
                        hommel=paste("Hommel, G. (1988). A stagewise rejective",
                            "multiple test procedure based on a modified",
                            "Bonferroni test. Biometrika 75, 383-386."),
                        hochberg=paste("Hochberg, Y. (1988). A sharper",
                            "Bonferroni procedure for multiple tests of",
                            "significance. Biometrika 75, 800-803."),
                        qvalue=paste("Storey, J.D., and Tibshirani, R. (2003).",
                            "Statistical significance for genomewide studies.",
                            "Proc Natl Acad Sci U S A 100, 9440-9445.")
                    ),
                    figure=list(
                        mds=paste("Planet, E., Attolini, C.S., Reina, O.,",
                            "Flores, O., and Rossell, D. (2012). htSeqTools:",
                            "high-throughput sequencing quality control,",
                            "processing and visualization in R. Bioinformatics",
                            "28, 589-590.")
                    )
                )
            )
        }
    )
    return(messages)
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

makeSimData <- function(nsnp=10000,nsam=100,nphe=3,csnp=100,
    psnp=c(0.01,0.05,0.1,0.15,0.20,0.25,0.3,0.4)) {
    # Checks
    if (!requireNamespace("PhenotypeSimulator"))
        stop("R package PhenotypeSimulator is required!")
    
    # Pseudo-FAM table
    proband <- seq(1,nsam,by=3)
    father <- mother <- affected <- sex <- rep(NA,nsam)
    father[proband] <- paste0("sample_",proband+1)
    mother[proband] <- paste0("sample_",proband+2)
    affected[proband] <- 2
    sex[proband] <- sample(c(1,2),length(proband),replace=TRUE)
    sex[proband+1] <- 1
    sex[proband+2] <- 2
    sex <- sex[seq_len(nsam)]
    pedigree=rep(seq_len(length(proband)),each=3)
    pedigree <- pedigree[seq_len(nsam)]
    fam <- data.frame(
        pedigree=pedigree,
        member=paste0("sample_",seq_len(nsam)),
        father=father,
        mother=mother,
        sex=sex,
        affected=affected,
        row.names=paste0("sample_",seq_len(nsam))
    )
    
    # Pseudo-MAP table
    nuc <- c("A","C","G","T")
    allele_a <- sample(nuc,nsnp,replace=TRUE)
    allele_b <- sample(nuc,nsnp,replace=TRUE)
    while(any(allele_a==allele_b)) {
        s <- which(allele_a==allele_b)
        replace=ifelse(length(s)<4,TRUE,FALSE)
        allele_a[s] <- sample(nuc,length(s),replace=TRUE)
    }
    map <- data.frame(
        chromosome=sample(seq_len(22),nsnp,replace=TRUE),
        snp.name=paste0("snp_",seq_len(nsnp)),
        position=sample(1000:1000000,nsnp),
        allele.1=allele_a,
        allele.2=allele_b,
        row.names=paste0("snp_",seq_len(nsnp))
    )

    # Simulated phenotypes and genotypes
    simObj <- .runSimulation(N=nsam,P=nphe,genVar=0.4,h2s=0.6,h2bg=0.4,rho=0.25,
        phi=0.25,delta=0.5,tNrSNP=nsnp,cNrSNP=csnp,SNPfrequencies=psnp)

    simpheno <- as.data.frame(simObj$phenoComponentsFinal$Y)
    rownames(simpheno) <- rownames(fam)
    geno <- simObj$rawComponents$genotypes$genotypes
    rownames(geno) <- rownames(fam)
    colnames(geno) <- rownames(map)
    
    return(list(
        snp=SnpMatrix(t(geno+1)),
        sample=fam,
        feature=map,
        pheno=simpheno
    ))
}

################################################################################

# Borrowed from PhenotypeSimulator
.runSimulation <- function(N,P,genVar=NULL,h2s=NULL,theta=0.8,h2bg=NULL,
    eta=0.8,noiseVar=NULL,rho=NULL,delta=NULL,gamma=0.8,phi=NULL,alpha=0.8,
    tNrSNP=5000,cNrSNP=20,SNPfrequencies=c(0.1,0.2,0.4),genotypefile=NULL,
    format='delim',genoFilePrefix=NULL,genoFileSuffix=NULL,genoDelimiter=",",
    skipFields=NULL,header=FALSE,probabilities=FALSE,chr=NULL,
    NrSNPsOnChromosome=NULL,NrChrCausal=NULL,kinshipfile=NULL,
    kinshipHeader=FALSE,kinshipDelimiter=",",standardise=TRUE,
    distBetaGenetic="norm",mBetaGenetic=0,sdBetaGenetic=1,
    pTraitsAffectedGenetics=1,pIndependentGenetic=0.4,
    pTraitIndependentGenetic=0.2,keepSameIndependentSNPs=FALSE,
    NrFixedEffects=1,NrConfounders=10,distConfounders="norm",mConfounders=0,
    sdConfounders=1,catConfounders=NULL,probConfounders=NULL,
    distBetaConfounders="norm",mBetaConfounders=0,sdBetaConfounders=1,
    pTraitsAffectedConfounders=1,pIndependentConfounders=0.4,
    pTraitIndependentConfounders=0.2,keepSameIndependentConfounders=FALSE,
    pcorr=0.8,corrmatfile=NULL,meanNoiseBg=0,sdNoiseBg=1,nonlinear=NULL,
    logbase=10,expbase=NULL,power=NULL,customTransform=NULL,transformNeg="abs",
    proportionNonlinear=0,sampleID="ID_",phenoID="Trait_",snpID="SNP_",
    seed=219453) {
    set.seed(seed)
    verbose <- FALSE
    model <- setModel(genVar=genVar,h2s=h2s,h2bg=h2bg,theta=theta,eta=eta,
        noiseVar=noiseVar,rho=rho,delta=delta,gamma=gamma,phi=phi,alpha=alpha,
        pcorr=pcorr,pIndependentConfounders=pIndependentConfounders,
        pTraitIndependentConfounders=pTraitIndependentConfounders,
        pIndependentGenetic=pIndependentGenetic,
        pTraitIndependentGenetic=pTraitIndependentGenetic,
        proportionNonlinear=proportionNonlinear,cNrSNP=cNrSNP,
        NrConfounders=NrConfounders,verbose=verbose)
    id_snps <- NULL
    id_samples <- NULL
    id_phenos <- paste(phenoID, 1:P, sep="")

    ### create simulated phenotypes
    # 1. Simulate genetic terms
    if (grepl('Fixed', model$modelGenetic)) {
        if (is.null(genoFilePrefix) && is.null(genotypefile)) {
            if (!grepl('Bg', model$modelGenetic) && tNrSNP != cNrSNP) {
                warning(paste0("The genetic model does not contain ",
                    "infinitesimal genetic effects but the total number of ",
                    "SNPs to simulate (tNrSNP: ",tNrSNP,") is larger than the ", 
                    "number of genetic variant effects SNPs (cNrSNP: ", cNrSNP,
                    "). If genotypes are not needed, consider setting tNrSNPs=",
                    "cNrSNPs to speed up computation"))
            }
            genotypes <- simulateGenotypes(N=N,NrSNP=tNrSNP, 
                frequencies=SNPfrequencies,sampleID=sampleID,snpID=snpID, 
                verbose=verbose)
            id_samples <- genotypes$id_samples
            id_snps <- genotypes$id_snps
        } 
        else if (! is.null(genotypefile)) {
            genotypes <- readStandardGenotypes(N=N,filename=genotypefile, 
                format=format,verbose=verbose,sampleID=sampleID,snpID=snpID, 
                delimiter=genoDelimiter)
            id_samples <- genotypes$id_samples
            id_snps <- genotypes$id_snps
        } 
        else
            genotypes <- NULL

        causalSNPs <- getCausalSNPs(N=N,NrCausalSNPs=cNrSNP,chr=chr, 
            NrChrCausal=NrChrCausal,NrSNPsOnChromosome=NrSNPsOnChromosome,
            genotypes=genotypes$genotypes,genoFilePrefix=genoFilePrefix, 
            genoFileSuffix=genoFileSuffix,format=format,
            probabilities=probabilities,skipFields=skipFields,header=header,
            delimiter=genoDelimiter,sampleID=sampleID,verbose=verbose)
        if (is.null(id_snps))
            id_snps <- colnames(causalSNPs)
        if (is.null(id_samples))
            id_samples <- rownames(causalSNPs)
            
        genFixed <- geneticFixedEffects(X_causal=causalSNPs,N=N,P=P, 
            pTraitsAffected=pTraitsAffectedGenetics,
            pIndependentGenetic=model$pIndependentGenetic, 
            pTraitIndependentGenetic=model$pTraitIndependentGenetic,
            keepSameIndependent=keepSameIndependentSNPs,
            distBeta=distBetaGenetic,mBeta=mBetaGenetic,sdBeta=sdBetaGenetic,
            id_phenos=id_phenos,id_samples=id_samples,phenoID=phenoID,
            verbose=verbose)

        var_genFixed_shared <- model$theta * model$h2s * model$genVar
        var_genFixed_independent <- (1 - model$theta) * model$h2s * model$genVar
        genFixed_shared_rescaled <- rescaleVariance(genFixed$shared, 
            var_genFixed_shared)
        genFixed_independent_rescaled <- rescaleVariance(genFixed$independent,
            var_genFixed_independent)
        Y_genFixed <- .addNonNulls(list(genFixed_shared_rescaled$component, 
            genFixed_independent_rescaled$component))
    } else {
        genFixed <- NULL
        genFixed_shared_rescaled <- NULL
        genFixed_independent_rescaled <- NULL
        genotypes <- NULL
        var_genFixed_shared <- 0
        var_genFixed_independent <- 0
        Y_genFixed <- NULL
        cNrSNP <- 0
    }
    if (grepl('Bg', model$modelGenetic)) {
        if (is.null(kinshipfile)) {
            if (is.null(genotypes) && !is.null(genotypefile)){
                genotypes <- readStandardGenotypes(genotypefile,format=format,
                    verbose=verbose,sampleID=sampleID,snpID=snpID, 
                    delimiter=genoDelimiter)
            }
            if (is.null(genotypes)) {
                genotypes <- simulateGenotypes(N=N,NrSNP=tNrSNP, 
                    frequencies=SNPfrequencies,sampleID=sampleID,
                    verbose=verbose)
            }
            id_samples <- genotypes$id_samples
            id_snps <- genotypes$id_snps
            
            kinship <- getKinship(X=genotypes$genotypes,N=N,
                standardise=standardise,id_samples=id_samples,sampleID=sampleID,
                verbose=verbose)
        } else {
            kinship <- getKinship(kinshipfile=kinshipfile,sep=kinshipDelimiter,
                N=N,header=kinshipHeader,verbose=verbose,sampleID=sampleID,
                id_samples=id_samples)
            if(!is.null(genotypes)) {
                if(colnames(kinship) != genotypes$id_samples) {
                    stop(paste("Sample names in kinship file do not match", 
                        "sample names of genotypes",sep=""))
                }
            } 
            else
                id_samples <- colnames(kinship)
        }

        if (eta == 1) {
            genBgShared <- TRUE
            genBgIndependent <- FALSE
        }
        else if (eta == 0) {
            genBgShared <- FALSE
            genBgIndependent <- TRUE
        } else {
            genBgShared <- TRUE
            genBgIndependent <- TRUE
        }

        genBg <- geneticBgEffects(N=N,P=P,kinship=kinship,shared=genBgShared, 
            independent=genBgIndependent,id_phenos=id_phenos)
        eval_kinship <- genBg$eval_kinship
        evec_kinship <- genBg$evec_kinship
        
        var_genBg_shared <- model$eta * model$h2bg * model$genVar
        var_genBg_independent <- (1 - model$eta) * model$h2bg * model$genVar

        genBg_shared_rescaled <- rescaleVariance(genBg$shared,var_genBg_shared)
        genBg_independent_rescaled <- rescaleVariance(genBg$independent,
            var_genBg_independent)

        Y_genBg <- .addNonNulls(list(genBg_shared_rescaled$component, 
            genBg_independent_rescaled$component))

        cov_genBg_shared <- genBg$cov_shared
        cov_genBg_shared_rescaled <- cov_genBg_shared * 
            genBg_shared_rescaled$scale_factor^2

        cov_genBg_independent <- genBg$cov_independent
        cov_genBg_independent_rescaled <- cov_genBg_independent * 
            genBg_independent_rescaled$scale_factor^2

        cov_genBg <- cov_genBg_shared_rescaled + cov_genBg_independent_rescaled
    } else {
        genBg <- NULL
        genBg_shared_rescaled <- NULL
        genBg_independent_rescaled <- NULL
        kinship <- NULL
        var_genBg_shared <- 0
        var_genBg_independent <- 0
        Y_genBg <- NULL
        cov_genBg <- NULL
        cov_genBg_shared <- NULL
        cov_genBg_independent <- NULL
        eval_kinship <- NULL
        evec_kinship <- NULL
    }
    # 1. Simulate noise terms
    if (grepl('Correlated', model$modelNoise))  {
        corr_mat <- NULL
        if (!is.null(corrmatfile)){
            corr_mat <- as.matrix(data.table::fread(corrmatfile,sep=",",
                data.table=FALSE,stringsAsFactors=FALSE))
        }
        correlatedBg <- correlatedBgEffects(N=N,P=P,corr_mat=corr_mat, 
            pcorr=model$pcorr,id_phenos=id_phenos,id_samples=id_samples,
            sampleID=sampleID,phenoID=phenoID)
        var_noiseCorrelated <- model$rho *  model$noiseVar
        correlatedBg_rescaled <- rescaleVariance(correlatedBg$correlatedBg, 
            var_noiseCorrelated)

        Y_correlatedBg <- correlatedBg_rescaled$component

        cov_correlatedBg <- correlatedBg$cov_correlated * 
            correlatedBg_rescaled$scale_factor^2
    } else {
        correlatedBg <- NULL
        var_noiseCorrelated <- 0
        Y_correlatedBg <- NULL 
        cov_correlatedBg <- NULL
    }
    if (grepl('Bg', model$modelNoise)) {
        if (alpha == 1) {
            noiseBgShared <- TRUE
            noiseBgIndependent <- FALSE
        }
        else if (alpha == 0) {
            noiseBgShared <- FALSE
            noiseBgIndependent <- TRUE
        } else {
            noiseBgShared <- TRUE
            noiseBgIndependent <- TRUE
        }
        noiseBg <- noiseBgEffects(N=N,P=P,mean=meanNoiseBg,sd=sdNoiseBg,
            shared=noiseBgShared,independent=noiseBgIndependent,
            id_phenos=id_phenos, id_samples=id_samples,sampleID=sampleID,
            phenoID=phenoID)

        var_noiseBg_shared <- model$alpha * model$phi * model$noiseVar
        var_noiseBg_independent <- (1 - model$alpha) * model$phi *model$noiseVar

        noiseBg_shared_rescaled <- rescaleVariance(noiseBg$shared, 
            var_noiseBg_shared)
        noiseBg_independent_rescaled <- rescaleVariance(noiseBg$independent, 
            var_noiseBg_independent)

        Y_noiseBg <- .addNonNulls(list(noiseBg_shared_rescaled$component, 
            noiseBg_independent_rescaled$component))

        cov_noiseBg_shared <- noiseBg$cov_shared
        cov_noiseBg_shared_rescaled <- cov_noiseBg_shared * 
            noiseBg_shared_rescaled$scale_factor^2

        cov_noiseBg_independent <- noiseBg$cov_independent
        cov_noiseBg_independent_rescaled <- cov_noiseBg_independent * 
            noiseBg_independent_rescaled$scale_factor^2

        cov_noiseBg <- cov_noiseBg_shared_rescaled + 
            cov_noiseBg_independent_rescaled

    } 
    else {
        noiseBg <- NULL
        noiseBg_shared_rescaled <- NULL
        noiseBg_independent_rescaled <- NULL
        var_noiseBg_shared <- 0
        var_noiseBg_independent <- 0
        Y_noiseBg <- NULL
        cov_noiseBg <- NULL
        cov_noiseBg_shared <- NULL
        cov_noiseBg_independent <- NULL
    }
    if (grepl('Fixed', model$modelNoise)) {
        noiseFixed <- noiseFixedEffects(P=P,N=N,NrFixedEffects=NrFixedEffects,
            NrConfounders=NrConfounders,
            pTraitsAffected=pTraitsAffectedConfounders,
            pIndependentConfounders=model$pIndependentConfounders,
            pTraitIndependentConfounders=model$pTraitIndependentConfounders,
            keepSameIndependent=keepSameIndependentConfounders,
            distConfounders=distConfounders,mConfounders=mConfounders, 
            sdConfounders=sdConfounders,catConfounders=catConfounders, 
            probConfounders=probConfounders,distBeta=distBetaConfounders, 
            mBeta=mBetaConfounders,sdBeta=sdBetaConfounders,id_phenos=id_phenos,
            id_samples=id_samples,sampleID=sampleID,phenoID=phenoID,
            verbose=verbose)

        var_noiseFixed_shared <- model$gamma * model$delta * model$noiseVar
        var_noiseFixed_independent <- (1 - model$gamma) * model$delta * 
            model$noiseVar

        noiseFixed_shared_rescaled <- rescaleVariance(noiseFixed$shared, 
            var_noiseFixed_shared)
        noiseFixed_independent_rescaled <- rescaleVariance(
            noiseFixed$independent,var_noiseFixed_independent)

        Y_noiseFixed <- .addNonNulls(list(noiseFixed_shared_rescaled$component, 
            noiseFixed_independent_rescaled$component))
    } else {
        noiseFixed <- NULL
        noiseFixed_shared_rescaled <- NULL
        noiseFixed_independent_rescaled <- NULL
        var_noiseFixed_shared <- 0
        var_noiseFixed_independent <- 0
        Y_noiseFixed <- NULL
    }

    # 3. Construct final simulated phenotype 
    components <- list(Y_genFixed,Y_genBg,Y_noiseFixed,Y_noiseBg,Y_correlatedBg)
    Y <- .addNonNulls(components)

    # 4. Transformation
    if (!is.null(nonlinear)) {
        Y_transformed <- transformNonlinear(Y,method=nonlinear, 
            alpha=proportionNonlinear,logbase=logbase,expbase=expbase,
            power=power,f=customTransform,transformNeg=transformNeg)
        Y_transformed <- scale(Y_transformed)
    } else
        Y_transformed <- NULL

    Y <- scale(Y)

    varComponents <- data.frame(genVar=model$genVar,h2s=model$h2s, 
        h2bg=model$h2bg,proportionNonlinear=model$proportionNonlinear,
        var_genFixed_shared=var_genFixed_shared, 
        var_genFixed_independent=var_genFixed_independent, 
        var_genBg_shared=var_genBg_shared, 
        var_genBg_independent=var_genBg_independent, 
        noiseVar=model$noiseVar,var_noiseFixed_shared=var_noiseFixed_shared, 
        var_noiseFixed_independent=var_noiseFixed_independent, 
        var_noiseBg_shared=var_noiseBg_shared, 
        var_noiseBg_independent=var_noiseBg_independent, 
        var_noiseCorrelated=var_noiseCorrelated)
    phenoComponentsFinal <- list(Y=Y,Y_genFixed=Y_genFixed,Y_genBg=Y_genBg, 
        Y_noiseFixed=Y_noiseFixed,Y_noiseBg=Y_noiseBg, 
        Y_correlatedBg=Y_correlatedBg,Y_transformed=Y_transformed,
        cov_genBg=cov_genBg,cov_noiseBg=cov_noiseBg,
        cov_correlatedBg=cov_correlatedBg)
    phenoComponentsIntermediate <- list(
        Y_genFixed_shared=genFixed_shared_rescaled$component,
        Y_genFixed_independent=genFixed_independent_rescaled$component,
        Y_noiseFixed_shared=noiseFixed_shared_rescaled$component,
        Y_noiseFixed_independent=noiseFixed_independent_rescaled$component,
        Y_genBg_shared=genBg_shared_rescaled$component,
        Y_genBg_independent=genBg_independent_rescaled$component,
        Y_noiseBg_shared=noiseBg_shared_rescaled$component,
        Y_noiseBg_independent=noiseBg_independent_rescaled$component,
        cov_genBg_shared=cov_genBg_shared, 
        cov_genBg_independent=cov_genBg_independent, 
        cov_noiseBg_shared=cov_noiseBg_shared,
        cov_noiseBg_independent=cov_noiseBg_independent, 
        genFixed=genFixed,
        noiseFixed=noiseFixed)

    setup <- list(P=P,N=N,NrCausalSNPs=cNrSNP,modelGenetic=model$modelGenetic, 
        modelNoise=model$modelNoise,id_samples=id_samples,id_phenos=id_phenos,
        id_snps=id_snps)
    rawComponents <- list(kinship=kinship,genotypes=genotypes,
        eval_kinship=eval_kinship,evec_kinship=evec_kinship,causal=causalSNPs)

    return(list(varComponents=varComponents, 
        phenoComponentsFinal=phenoComponentsFinal, 
        phenoComponentsIntermediate=phenoComponentsIntermediate, 
        setup=setup,rawComponents=rawComponents))
}

.addNonNulls <- function(compList) {
    if (!is.list(compList))
        stop("addNonNulls expects input of type list")
    nonNulls <- compList[!sapply(compList, is.null)]
    if (length(nonNulls) == 0) 
        return(NULL)
    if (length(unique(sapply(nonNulls, ncol))) != 1)
        stop("Column dimensions of list elements are different")
    if (length(unique(sapply(nonNulls, nrow))) != 1)
        stop("Row dimensions of list elements are different")
    Reduce("+",nonNulls)
}
