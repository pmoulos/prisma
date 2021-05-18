test_that("GWASExperiment from snpStats PLINK works",{
    input <- .gimmeTestFiles()
    sample <- snpStats::read.plink(input$bed,input$bim,input$fam)
    
    # Without phenotypes
    gwe <- importGWAS(input,backend="snpStats")
    expect_true(is(gwe,"GWASExperiment"))
    expect_equal(nrow(gwe),nrow(sample$map))
    expect_equal(ncol(gwe),nrow(sample$fam))
    expect_identical(rowData(gwe),DataFrame(sample$map))
    expect_identical(colData(gwe),DataFrame(sample$fam))
    expect_identical(phenotypes(gwe),NULL)
    expect_identical(metadata(gwe)$backend,"snpStats")
    
    # With phenotypes
    set.seed(42)
    pseudopheno <- data.frame(
        case_control=sample(c(0,1),nrow(sample$fam),replace=TRUE),
        other_pheno=sample(c("drug","nodrug"),nrow(sample$fam),replace=TRUE),
        yap=sample(c("normal","demipsek","psek"),nrow(sample$fam),replace=TRUE),
        cont=round(runif(nrow(sample$fam)),3),
        row.names=rownames(sample$fam)
    )
    gwe <- importGWAS(input,phenos=pseudopheno,backend="snpStats")
    expect_true(is(gwe,"GWASExperiment"))
    expect_identical(phenotypes(gwe),pseudopheno)
})

test_that("Basic snpStats filters work",{
    fam <- system.file("extdata/sample.fam",package="snpStats")
    bim <- system.file("extdata/sample.bim",package="snpStats")
    bed <- system.file("extdata/sample.bed",package="snpStats") 
    input <- list(bed=bed,bim=bim,fam=fam)
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
