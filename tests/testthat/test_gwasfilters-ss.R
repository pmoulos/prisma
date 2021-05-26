test_that("Basic snpStats filters work",{
    input <- .gimmeTestFiles()
    gwe <- importGWAS(input,backend="snpStats")
    
    filts <- list(
        snpCallRate=0.98,
        sampleCallRate=0.95,
        maf=0.05,
        hwe=1e-3,
        LD=0.2,
        heteroStat="median",
        heteroFac=3,
        heteroHard=NA,
        pcaOut=TRUE,
        IBD=0.2,
        inbreed=0.1
    )
    
    # Preparation
    y <- t(assay(gwe,1))
    snpSumY <- col.summary(y)
    sampleSumY <- row.summary(y)
    loc <- mean(sampleSumY$Heterozygosity,na.rm=TRUE)
    sca <- sd(sampleSumY$Heterozygosity,na.rm=TRUE)
    hetF <- .calcInbreedFromSnpMatrix(y,snpSumY,sampleSumY)
    
    # Test
    gweFilt <- .filterWithSnpStatsBasic(gwe,filts)
    expect_equal(nrow(gweFilt),5)
    expect_equal(ncol(gweFilt),37)
    
    x <- t(assay(gweFilt,1))
    snpSum <- col.summary(x)
    sampleSum <- row.summary(x)
    
    expect_true(all(snpSum$Call.rate > filts$snpCallRate))
    expect_true(all(snpSum$MAF > filts$maf))
    expect_true(all(abs(snpSum$z.HWE) < abs(qnorm(filts$hwe/2))))
    expect_true(all(sampleSum$Call.rate > filts$sampleCallRate))
    expect_true(all(sampleSum$Heterozygosity > loc - filts$heteroFac*sca
        & sampleSum$Heterozygosity < loc + filts$heteroFac*sca))
    expect_true(all(sampleSum$Call.rate > filts$sampleCallRate))
    expect_equal(length(which(hetF > filts$inbreed)),64)
})

test_that("IBD filter and SNPRelate PCA calculation work",{
    input <- .gimmeTestFiles()
    gwe <- importGWAS(input,backend="snpStats")
    
    filts <- getDefaults("filters")
    # Change this as number of SNPs in test make it really strict
    filts$IBD <- 0.45
    
    # Test
    gweFilt <- .filterWithSnpStatsIbd(gwe,filts,.testing=TRUE)
    expect_equal(nrow(gweFilt),20)
    expect_equal(ncol(gweFilt),7)
    
    m <- metadata(gweFilt)
    expect_true(!is.null(m$pcaOut))
    expect_true(is(m$pcaOut,"snpgdsPCAClass"))
})

test_that("Robust sample PCA filtering works",{
    input <- .gimmeTestFiles()
    gwe <- importGWAS(input,backend="snpStats")
    
    filts <- getDefaults("filters")
    # Do not IBD
    filts$IBD <- NA
    
    # Test PcaGrid
    gweFilt1 <- .filterWithSnpStatsRobustPca(gwe,filts)
    expect_equal(nrow(gweFilt1),20)
    expect_equal(ncol(gweFilt1),86)
    
    m <- metadata(gweFilt1)
    expect_true(!is.null(m$pcaRob))
    expect_true(is(m$pcaRob,"PcaGrid"))
    
    # Test PcaHubert
    filts$pcaRobust <- "hubert"
    gweFilt2 <- .filterWithSnpStatsRobustPca(gwe,filts)
    expect_equal(nrow(gweFilt2),20)
    expect_equal(ncol(gweFilt2),93)
    
    m <- metadata(gweFilt2)
    expect_true(!is.null(m$pcaRob))
    expect_true(is(m$pcaRob,"PcaHubert"))
})

test_that("Imputation with snpStats works",{
    input <- .gimmeTestFiles()
    gwe <- importGWAS(input,backend="snpStats")
    
    expect_true(any(is.na(assay(gwe,1))))
    
    # Small hack as the dataset is small and non-heterogeneous enough
    x <- assay(gwe,1)
    x["TSC0101718","430"] <- as.raw(1)
    assay(gwe,1) <- x

    o1 <- .internalImputeWithSnpStats(gwe)
    expect_false(any(is.na(assay(o1,1))))
})

test_that("kNN impute works",{
    input <- .gimmeTestFiles()
    gwe <- importGWAS(input,backend="snpStats")
    a <- as(assay(gwe,1),"numeric")
    
    expect_true(any(is.na(a)))
    b <- .internalImputeKnn(a)
    expect_false(any(is.na(b)))
})

test_that(".checkSelection works",{
    s1 <- NULL
    expect_silent(.checkSelection(s1))
    
    s2 <- list(alpha=2:4,beta=3:5)
    expect_error(.checkSelection(s2))
    
    s3 <- list(samples=2:4,beta=3:5)
    expect_error(.checkSelection(s3))
    
    s4 <- list(samples=2:4,snps=c("foo","bar"))
    expect_error(.checkSelection(s4))
    
    s5 <- list(samples=2:4,snps=3:5)
    expect_silent(.checkSelection(s5))
})

test_that(".checkFilters works",{
    f1 <- getDefaults("filters")
    expect_silent(.checkFilters(f1))
    
    f2 <- list(snpCallRate=0.98,maf=0.05,hwe=1e-6,IBD=0.1)
    expect_silent(.checkFilters(f2))
    f21 <- .checkFilters(f2)
    expect_true("sampleCallRate" %in% names(f21))
    expect_equal(f21$sampleCallRate,f1$sampleCallRate)
    
    f3 <- list(snpCallRat=0.98,maf=0.05,hwe=1e-6)
    expect_warning(.checkFilters(f3))
    
    f4 <- list(snpCallRat=0.98,ma=0.05,we=1e-6)
    expect_error(.checkFilters(f4))
    
    f5 <- list(snpCallRate=0.98,maf=0.05,heteroStat="median",heteroFac=NULL)
    expect_warning(.checkFilters(f5))
    
    f6 <- list(snpCallRate=0.98,maf=0.05,heteroStat=NULL,heteroFac=3)
    expect_warning(.checkFilters(f6))
    
    f7 <- list(snpCallRate=0.98,heteroStat="median",heteroFac=3,heteroHard=0.1)
    expect_warning(.checkFilters(f7))
})
