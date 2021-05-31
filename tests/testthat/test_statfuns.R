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

test_that("partitionGWAS works",{
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
