combineSimes <- function(p,zerofix=NULL) {
    p <- .zeroFix(p,zerofix)
    m <- length(p)
    y <- sort(p)
    s <- min(m*(y/(seq_len(m))))
    return(min(c(s,1)))
}

#combineWeight <- function(p,w,zerofix=NULL) {
#    p <- .zeroFix(p,zerofix)
#    return(prod(p^w))
#}
#

combineHarmonic <- function(p,w,zerofix=NULL) {
    if (!requireNamespace("harmonicmeanp"))
        stop("R package harmonicmeanp is required!")
    p <- .zeroFix(p,zerofix)
    return(p.hmp(p,w,L=length(p),multilevel=FALSE))
}

combineMinp <- function(p) { 
    return(min(p))
}

combineMaxp <- function(p) { 
    return(max(p))
}

fisherMethod <- function(pvals,method=c("fisher"),p.corr=c("bonferroni","BH",
    "none"),zeroSub=0.00001,na.rm=FALSE) {
    stopifnot(method %in% c("fisher"))
    stopifnot(p.corr %in% c("none","bonferroni","BH"))
    stopifnot(all(pvals>=0, na.rm=TRUE) & all(pvals<=1, na.rm=TRUE))
    stopifnot(zeroSub>=0 & zeroSub<=1 || length(zeroSub)!=1)
    if(is.null(dim(pvals)))
        stop("pvals must have a dim attribute")
    p.corr <- ifelse(length(p.corr)!=1, "BH", p.corr)
    ##substitute p-values of 0
    pvals[pvals == 0] <- zeroSub
    
    fisher.sums <- data.frame(do.call(rbind,apply(pvals,1,fisherSum,
        zeroSub=zeroSub,na.rm=na.rm)))
        
    rownames(fisher.sums) <- rownames(pvals)
    fisher.sums$p.value <- 1-pchisq(fisher.sums$S,df=2*fisher.sums$num.p)
    fisher.sums$p.adj <- switch(p.corr,
        bonferroni = p.adjust(fisher.sums$p.value,"bonferroni"),
        BH = p.adjust(fisher.sums$p.value,"BH"),
        none = fisher.sums$p.value
    )
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
