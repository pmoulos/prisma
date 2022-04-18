calcPcaCovar <- function(obj,ld=0.2,method=c("auto","snprel","grid","hubert"),
    npc=0,pval=0.05,rc=NULL) {
    if (!is(obj,"GWASExperiment"))
        stop("obj must be an object of class GWASExperiment!")
    
    method <- tolower(method[1])
    .checkTextArgs("PCA method",method,c("auto","snprel","grid","hubert"),
        multiarg=FALSE)
    .checkNumArgs("Number of PCs npc",npc,"numeric",0,"gte")
    if (!.isEmpty(ld))
        .checkNumArgs("Linkage disequilibrium",ld,"numeric",0,"gt")
    else
        ld <- NA
    
    disp("Performing PCA to select covariates to include in associations")
    switch(method,{
        auto = 
            # Check if pcaCov exists in obj
            # 1. Redo LD pruning with SNPRelate
            # 2. Robust PCA with the LD-pruned SNPs
            # We can use the .wrapPcaWithSnpStatsLd combined function
            pseudof <- .guessPcaParamsFromObject(obj,ld)
            obj <- .wrapPcaWithSnpStatsLd(obj,pseudof,rc)
        },
        snprel = {
            obj <- .wrapPcaWithSnpStatsLd(obj,list(LD=ld),rc)
        },
        grid = {
            obj <- .wrapPcaWithSnpStatsLd(obj,list(LD=ld,pcaRobust="grid"),rc)
        },
        hubert = {
            obj <- .wrapPcaWithSnpStatsLd(obj,list(LD=ld,pcaRobust="hubert"),rc)
        }
    )
    
    disp("Selecting number of PCs to include")
    m <- metadata(obj)
    pco <- m$pcaCov
    if (npc == 0) {
        disp("  with Tracy-Widom test for significant eigenvalues")
        if (is(pco,"PcaGrid") || is(pco,"PcaHubert"))
            eigv <- getEigenvalues(pco)
        else
            eigv <- pco$eigenval[!is.na(pco$eigenval)]
        
        twpcList <- twTest(eigv,pval)
        ei <- twpcList$index
        
        if (is(pco,"PcaGrid") || is(pco,"PcaHubert")) {
            thePcs <- getLoadings(pco)
            thePcs <- thePcs[,ei,drop=FALSE]
        }
        else {
            thePcs <- as.data.frame(pco$eigenvect[,ei,drop=FALSE])
            rownames(thePcs) <- as.character(pco$sample.id)
            colnames(thePcs) <- paste0("PC",ei)
        }
        disp("    selected ",length(ei)," PCs")
    }
    else {
        disp("  with hard-cutoff at ",npc," PCs")
        if (is(pco,"PcaGrid") || is(pco,"PcaHubert"))
            thePcs <- getLoadings(pco)
        else {
            thePcs <- as.data.frame(pco$eigenvect)
            rownames(thePcs) <- as.character(pco$sample.id)
            colnames(thePcs) <- paste0("PC",seq_len(ncol(thePcs)))
        }
        
        if (npc > ncol(thePcs)) {
            warning("The number of desired covariate PCs is larger than the ",
                "number of available eigenvectors!\nWill use the available ",
                "ones (",ncol(thePcs),")...",immediate.=TRUE)
            npc <- ncol(thePcs)
        }

        thePcs <- thePcs[,seq_len(npc),drop=FALSE]
    }
    
    pcaCovariates(obj) <- thePcs
    return(obj)
}

.guessPcaParamsFromObject <- function(obj,ld=0.2) {
    m <- metadata(obj)
    if ("pcaCov" %in% names(m)) {
        if (is(m$pcaCov,"PcaGrid"))
            pseudof <- list(LD=ld,pcaRobust="grid")
        else if (is(m$pcaCov,"PcaGrid"))
            pseudof <- list(LD=ld,pcaRobust="hubert")
        else
            pseudof <- list(LD=ld)
    }
    else
        pseudof <- list(LD=ld)
    return(pseudof)
}

twTest <- function(eigv,p=0.05,tol=1e-8) {
    if (!requireNamespace("RMTstat"))
        stop("R package RMTstat is required!")
    
    .checkNumArgs("p-value cutoff",p,"numeric",c(0,1),"both")
    .checkNumArgs("Tolereance tol",tol,"numeric",0,"gte")
    
    eigv[eigv<=tol] <- tol
    
    # Calculation of the TW statistic
    L1 <- rev(cumsum(rev(eigv)))
    L2 <- rev(cumsum(rev(eigv^2)))
    N <- rev(seq_len(length(eigv)))
    S2 <- N^2*L2/(L1^2)
    v <- N*(N+2)/(S2-N) # Effective number of markers
    
    L <- N*eigv/L1

    vSt <- sqrt(v-1)
    Nst <- sqrt(N)

    mu  <- (vSt+Nst)^2/v
    sig <- (vSt+Nst)/v * (1/vSt+1/Nst)^(1/3)

    twstat <- (L-mu)/sig
    
    ps <- 1-ptw(twstat[!is.na(twstat)])
    d <- which(ps<p)
    if (length(d) == 0) {
        warning("No signficant PCs found! Will return the first one...")
        d <- 1
    }
    
    return(list(
        tw=twstat[!is.na(twstat)],
        pval=ps,
        index=d
    ))
}

# TODO: Plot histogram with overlayed normal curve + Q-Q plot
normalityCheck <- function(x,pval=0.05,tests=c("sw","ks","jb"),lower=30,
    #forceTest=TRUE,
    combine=c("fisher","simes","max","min"),nsample=1000) {
    # upper not implemented
    if ("jb" %in% tests && !requireNamespace("tseries"))
        stop("R package tseries is required for Jarque-Bera test!")
    
    combine <- combine[1]
    .checkTextArgs("Normality tests",tests,c("sw","ks","jb"),multiarg=TRUE)
    .checkTextArgs("p-value combination",combine,
        c("fisher","simes","max","min"),multiarg=TRUE)
    .checkNumArgs("Lower sample size",lower,"numeric",0,"gt")
    .checkNumArgs("Resampling iterations",nsample,"numeric",0,"gt")
    #if (!is.logical(forceTest))
    #    stop("Testing above upper sample size (forceTest) must be logical!")
    
    n <- length(x)
    p <- s <- numeric(length(tests))
    names(p) <- names(s) <- tests
    if (n <= lower) {
        disp("Sample size ",n," <= ",lower,". Will execute combined ",
            "normality checks using requested tests.")
        if ("sw" %in% tests) {
            res <- shapiro.test(x)
            p["sw"] <- res$p.value
            s["sw"] <- res$statistic
        }
        if ("ks" %in% tests) {
            res <- ks.test(x,y="pnorm",mean=mean(x),sd=sd(x))
            p["ks"] <- res$p.value
            s["ks"] <- res$statistic
        }
        if ("jb" %in% tests) {
            xx <- x[!is.na(x)]
            res <- jarque.bera.test(xx)
            p["jb"] <- res$p.value
            s["jb"] <- res$statistic
        }
        
        if (all(p < pval)) # Done
            pp <- max(p)
        else {
            if (length(tests) > 1) {
                switch(combine,
                    min = {
                        pp <- combineMinp(p)
                    },
                    max = {
                        pp <- combineMaxp(p)
                    },
                    fisher = {
                        tmp <- fisherMethod(t(as.matrix(p)),p.corr="none",
                            zeroSub=.Machine$double.xmin)
                        pp <- tmp$p.value
                    },
                    simes = {
                        pp <- combineSimes(p)
                    }
                )
            }
            else
                pp <- p
        }
        
        return(list(
            statistic=s,
            pvalues=p,
            pval=pp,
            normal=pp>=pval
        ))
    }
    else { # Bootstrap per 50 observations covering all x ranges
        disp("Sample size ",n," > ",lower,". Will execute multiple ",
            "normality checks using the Shapiro-Wilk tests in subsamples of ",
            "size ",lower,".\n  resampling...")
        ind <- createSplit(y=x,n=nsample,frac=round(50/n,digits=3))
        res <- lapply(ind,function(i,s) {
            tmp <- shapiro.test(s[i])
            return(list(p=tmp$p.value,s=tmp$statistic))
        },x)
        ps <- unlist(lapply(res,function(x) return(x$p)),use.names=FALSE)
        ss <- unlist(lapply(res,function(x) return(x$s)),use.names=FALSE)
        
        # Solution with quantiles!
        # If 90% (or 100*(1-2*pval)%) of p-value distribution is larger than
        # pval, then we conclude that the distribution is normal and we report
        # quantile(ps,1-pval) as p-value
        # Otherwise, we conclude that the distribution is not normal and we
        # report quantile(ps,pval) as p-value
        pctCut <- 1 - 2*pval
        pctAct <- length(which(ps>=pval))/nsample
        if (pctAct > pctCut)
            p <- quantile(ps,1-pval)
        else
            p <- quantile(ps,pval)
        
        return(list(
            statistic=ss,
            pvalues=ps,
            pval=p,
            normal=p>=pval
        ))
    }
}

rankTransform <- function(x) {
    x <- x[!is.na(x)]
    z <- (x-mean(x))/sd(x)
    out <- rank(z) - 0.5
    out[is.na(z)] <- NA
    mP <- 0.5/max(out,na.rm=TRUE)
    out <- out/(max(out,na.rm=TRUE)+0.5)
    return(qnorm(out))
}

combineSimes <- function(p,zerofix=NULL) {
    p <- .zeroFix(p,zerofix)
    m <- length(p)
    y <- sort(p)
    s <- min(m*(y/(seq_len(m))))
    return(min(c(s,1)))
}

combineWeight <- function(p,w=NULL,zerofix=NULL) {
    p <- .zeroFix(p,zerofix)
    if (is.null(w))
        w <- rep(1/length(p),length(p))
    return(prod(p^w))
}

combineHarmonic <- function(p,w=NULL,zerofix=NULL) {
    if (!requireNamespace("harmonicmeanp"))
        stop("R package harmonicmeanp is required!")
    p <- .zeroFix(p,zerofix)
    if (is.null(w))
        w <- rep(1/length(p),length(p))
    return(p.hmp(p,w,L=length(p),multilevel=FALSE))
}

combineMinp <- function(p) { 
    return(min(p))
}

combineMaxp <- function(p) { 
    return(max(p))
}

fisherMethod <- function(pvals,zerofix=NULL) {
    stopifnot(all(pvals>=0 & all(pvals<=1)))
    
    if (is.null(dim(pvals)))
        stop("pvals must have a dim attribute")
    
    zeroSub <- min(apply(pvals,1,.zeroFix,zerofix))
    pvals[pvals == 0] <- zeroSub
    
    fisher.sums <- data.frame(do.call(rbind,apply(pvals,1,fisherSum,
        zeroSub=zeroSub,na.rm=TRUE)))
    rownames(fisher.sums) <- rownames(pvals)
    fisher.sums$p.value <- 1-pchisq(fisher.sums$S,df=2*fisher.sums$num.p)
    fisher.sums$p.adj <- fisher.sums$p.value
    return(fisher.sums)
}

fisherSum <- function(p,zeroSub=0.00001,na.rm=FALSE) {
    if(any(p>1, na.rm=TRUE)||any(p<0, na.rm=TRUE))
        stop("You provided bad p-values")
    stopifnot(zeroSub>=0 & zeroSub<=1 || length(zeroSub)!=1)
    p[p==0] <- zeroSub
    if (na.rm)
        p <- p[!is.na(p)]
    S = -2*sum(log(p))
    res <- data.frame(S=S,num.p=length(p))
    return(res)
}

.zeroFix <- function(p,z=NULL) {
    if (!is.null(z) && !is.numeric(z))
        stop("zerofix must be NULL or a numeric greater than 0 and less than ",
            "1!")
    if (!is.null(z) && (z <= 0 || z >= 1))
        stop("When zerofix is not NULL it must be a numeric greater than 0 ",
            "and less than 1!")
    #ze <- which(p==0)
    ze <- which(p < 1e-300) # A very small value as 0 causes also problems...
    if (length(ze)>0) {
        if (!is.null(z))
            p[ze] <- z*min(p[-ze])
        else
            p[ze] <- 0.5*runif(length(ze))*min(p[-ze])
    }
    return(p)
}
