# Preflight - 100 samples is the minimum for SNPTEST to run
dat <- makeSimData(nsnp=100,nsam=100,nphe=3,csnp=10)
gwe <- GWASExperiment(
    genotypes=dat$snp,
    features=dat$feature,
    samples=dat$sample,
    phenotypes=dat$pheno
)
map <- gfeatures(gwe)
gwe <- gwe[order(map$chromosome,map$position),]
gdsgwe <- writeGdsFile(gwe)
gdsfile(gwe) <- gdsgwe
gwe <- calcPcaCovar(gwe,method="hubert",npc=1)

test_that("gwa works",{
    gwe1 <- gwa(gwe,response="Trait_1",covariates=c("Trait_2","Trait_3"),
        pcs=TRUE,psig=0.05,methods=c("glm","statgen"))
    expect_true(all(c("glm","statgen") %in% colnames(pvalues(gwe1))))
})

test_that("gwaGlm works",{
    glmOut1 <- gwaGlm(gwe,response="Trait_1",covariates=c("Trait_2","Trait_3"),
        pcs=FALSE)
    glmOut2 <- gwaGlm(gwe,response="Trait_1",covariates=c("Trait_2","Trait_3"),
        pcs=TRUE)
    p1 <- glmOut1[,4] < 0.05
    p2 <- glmOut2[,4] < 0.05
    expect_false(identical(p1,p2))
})

test_that("gwaBlup works",{
    # Just test that run...
    expect_warning(bluOut1 <- gwaBlup(gwe,response="Trait_1",
        covariates=c("Trait_2","Trait_3"),usepc="rrint",npcs=2))
    expect_warning(bluOut3 <- gwaBlup(gwe,response="Trait_1",
        covariates=c("Trait_2","Trait_3"),usepc="none"))
})

test_that("gwaStatgen works",{
    stOut1 <- gwaStatgen(gwe,response="Trait_1",
        covariates=c("Trait_2","Trait_3"),pcs=FALSE)
    stOut2 <- gwaStatgen(gwe,response="Trait_1",
        covariates=c("Trait_2","Trait_3"),pcs=TRUE)
    p1 <- stOut1[,3] < 0.05
    p2 <- stOut2[,3] < 0.05
    expect_false(identical(p1,p2))
})

test_that("gwaSnptest works",{
    # For the time being, just that runs
    .findTool("snptest")
    skip_if(is.null(.EXTERNALS[["snptest"]]$exec))
    snpOut1 <- gwaSnptest(gwe,response="Trait_1",
        covariates=c("Trait_2","Trait_3"),pcs=FALSE)
    snpOut2 <- gwaSnptest(gwe,response="Trait_1",
        covariates=c("Trait_2","Trait_3"),pcs=TRUE)
    p1 <- snpOut1[,3] < 0.05
    p2 <- snpOut2[,3] < 0.05
    expect_false(identical(p1,p2))
})

test_that("gwaPlink works",{
    # For the time being, just that runs
    .findTool("plink")
    skip_if(is.null(.EXTERNALS[["plink"]]$exec))
    plOut1 <- gwaPlink(gwe,response="Trait_1",covariates=c("Trait_2","Trait_3"),
        pcs=FALSE)
    plOut2 <- gwaPlink(gwe,response="Trait_1",covariates=c("Trait_2","Trait_3"),
        pcs=TRUE)
    p1 <- plOut1[,3] < 0.05
    p2 <- plOut2[,3] < 0.05
    expect_false(identical(p1,p2))
})

test_that(".gwaGlmWorker works",{
    # Simple dummy data test
    expect_true(TRUE)
})

test_that(".canRunGwa works",{
    input <- .gimmeTestFiles()
    err <- importGWAS(input,backend="snpStats",writeGds=FALSE)
    expect_error(.canRunGwa(err))
})

test_that(".prepareGenotypesForBlup works",{
    # Maybe with a real test dataset or a test "thing"
    expect_true(TRUE)
})

test_that(".estimateNPCinBlup works",{
    # How?
    expect_true(TRUE)
})

test_that(".preparePlinkInputForSnptest works",{
    # Maybe with a real test dataset or a test "thing"
    expect_true(TRUE)
})

test_that(".preparePlinkInputForSnptest works",{
    # Maybe with a real test dataset or a test "thing"
    expect_true(TRUE)
})

test_that(".initSnptestSampleFirstRow works",{
    # Maybe with a real test dataset or a test "thing"
    expect_true(TRUE)
})

test_that(".validateResponseAndCovariates works",{
    df <- as.data.frame(matrix(runif(40),10,4))
    names(df) <- c("y","x1","x2","x3")
    
    res1 <- "y"
    cvs1 <- c("x1","x2","x3")
    v1 <- .validateResponseAndCovariates(df,res1,cvs1)
    expect_equal(v1$res,res1)
    expect_equal(v1$cvs,cvs1)
    
    res2 <- 1
    cvs2 <- c(2,3)
    v2 <- .validateResponseAndCovariates(df,res2,cvs2)
    expect_equal(v2$res,names(df)[res2])
    expect_equal(v2$cvs,names(df)[cvs2])
    
    res3 <- "y1"
    cvs3 <- c("x1","x2")
    expect_error(.validateResponseAndCovariates(df,res3,cvs3))
    
    res4 <- "y"
    cvs4 <- c("x1","x2","x4")
    expect_error(.validateResponseAndCovariates(df,res4,cvs4))
    
    res5 <- 8
    cvs5 <- c("x1","x2")
    expect_error(.validateResponseAndCovariates(df,res5,cvs5))
})

test_that(".validateBinaryForBinomial works",{
    set.seed(42)
    
    x1 <- runif(10)
    expect_error(x1 <- .validateBinaryForBinomial(x1))
    
    x2 <- sample(c(0,1),10,replace=TRUE)
    expect_silent(x2 <- .validateBinaryForBinomial(x2))
    
    x3 <- sample(c(1,2),10,replace=TRUE)
    expect_warning(x3 <- .validateBinaryForBinomial(x3))
    expect_true(is.factor(x3))
})
