test_that("twTest works",{
    set.seed(42)
    a <- matrix(runif(1000),50,20)
    b <- PcaHubert(a,k=10)
    e <- getEigenvalues(b)
    
    # OK - random dataset
    expect_warning(tw <- twTest(e))
    expect_equal(tw$index,1)
    expect_true(all(tw$pval > 0.05))
})

test_that("calcPcaCovar works",{
    input <- .gimmeTestFiles()
    gwe <- importGWAS(input,backend="snpStats")
    n1 <- ceiling(ncol(gwe)/3)
    n2 <- 2*floor(ncol(gwe)/3)
    ph <- data.frame(x=round(runif(ncol(gwe))),y=c(rep("A",n1),rep("B",n2)),
        row.names=colnames(gwe))
    phenotypes(gwe) <- ph
    
    expect_true(is(gwe,"GWASExperiment"))
    
    # Warning OK - small dataset
    expect_warning(o1 <- calcPcaCovar(obj=gwe,ld=NA,method="grid",npc=0))
    m <- metadata(o1)
    expect_true(is(m$pcaCov,"PcaGrid"))
    pc <- pcaCovariates(o1)
    expect_equal(ncol(pc),1)
    
    expect_warning(o2 <- calcPcaCovar(obj=gwe,ld=NA,method="hubert",npc=0))
    m <- metadata(o2)
    expect_true(is(m$pcaCov,"PcaHubert"))
    pc <- pcaCovariates(o2)
    expect_equal(ncol(pc),1)
    
    o3 <- calcPcaCovar(obj=gwe,ld=NA,method="grid",npc=2)
    m <- metadata(o3)
    expect_true(is(m$pcaCov,"PcaGrid"))
    pc <- pcaCovariates(o3)
    expect_equal(ncol(pc),2)
    
    o4 <- calcPcaCovar(obj=gwe,ld=NA,method="auto",npc=3)
    m <- metadata(o4)
    expect_true(is(m$pcaCov,"PcaGrid"))
    pc <- pcaCovariates(o4)
    expect_equal(ncol(pc),3)
})

test_that("normalityCheck works",{
    set.seed(42)

    # Simple test with default p-value combination (Fisher)
    x1 <- rnorm(20)
    
    o1 <- normalityCheck(x1)
    expect_true(o1$normal)
    
    # Simple test with another p-value combination (Simes)
    o2 <- normalityCheck(x1,combine="simes")
    expect_true(o2$normal)
    expect_true(o2$pval != o1$pval)
    
    # Simple test with another p-value combination (Union)
    o3 <- normalityCheck(x1,combine="min")
    expect_true(o3$normal)
    expect_true(o3$pval != o1$pval)
    expect_true(o3$pval != o2$pval)
    expect_equal(o3$pval,min(o3$pvalues))
    
    # Test with 2 instead of 3 tests
    o4 <- normalityCheck(x1,tests=c("sw","jb"),combine="simes")
    expect_true(o4$normal)
    expect_equal(length(o4$pvalues),2)
    expect_equal(length(o4$statistic),2)
    expect_true(o4$pval != o3$pval)
    
    # Test with 1 test
    o5 <- normalityCheck(x1,tests="sw")
    expect_true(o5$normal)
    expect_equal(length(o5$pvalues),1)
    expect_equal(length(o5$statistic),1)
    expect_equal(o5$pval,o5$pvalues)
    
    # Larger sample size
    x2 <- rnorm(100)
    o6 <- normalityCheck(x2)
    expect_true(o6$normal)
    expect_equal(length(o6$pvalues),1000)
    
    # Different lower parameter and also check resample
    x3 <- rnorm(50)
    o7 <- normalityCheck(x3,lower=50)
    expect_true(o7$normal)
    
    o8 <- normalityCheck(x2,lower=50,nsample=500)
    expect_true(o8$normal)
    expect_equal(length(o8$pvalues),500)
    
    # Not normally distributed
    x4 <- runif(50)
    o9 <- normalityCheck(x4,lower=50)
    expect_true(!o9$normal)
    
    # Not normally distributed - larger sample
    x5 <- runif(100)
    o10 <- normalityCheck(x5)
    expect_true(!o10$normal)
})
