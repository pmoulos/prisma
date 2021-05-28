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
