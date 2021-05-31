calcPcaCovar <- function(obj,ld=0.2,method=c("auto","snprel","grid","hubert"),
    npc=0,pval=0.05,rc=NULL) {
    if (!is(obj,"GWASExperiment"))
        stop("obj must be an object of class GWASExperiment!")
    
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
            m <- metadata(obj)
            # 1. Redo LD pruning with SNPRelate
            # 2. Robust PCA with the LD-pruned SNPs
            # We can use the .wrapPcaWithSnpStatsLd combined function               
            if ("pcaCov" %in% names(m)) {
                if (is(m$pcaCov,"PcaGrid"))
                    pseudof <- list(LD=ld,pcaRobust="grid")
                else if (is(m$pcaCov,"PcaGrid"))
                    pseudof <- list(LD=ld,pcaRobust="hubert")
                else
                    pseudof <- list(LD=ld)
            }
            else
                pseudof <- list(LD=ld,pcaRobust="grid")
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
    if (npc == 0) {
        disp("  with Tracy-Widom test for significant eigenvalues")
        m <- metadata(obj)
        pco <- m$pcaCov
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

    twstat <-(L-mu)/sig
    
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
