test_that("partitionGWAS works",{
    input <- .gimmeTestFiles()    
    gwe <- importGWAS(input,backend="snpStats")
    n1 <- ceiling(ncol(gwe)/3)
    n2 <- 2*floor(ncol(gwe)/3)
    ph <- data.frame(x=round(runif(ncol(gwe))),y=c(rep("A",n1),rep("B",n2)),
        row.names=colnames(gwe))
    phenotypes(gwe) <- ph
    
    expect_true(is(gwe,"GWASExperiment"))
    
    expect_warning(o1 <- partitionGWAS(gwe,n=10))
    expect_true(is(o1,"list"))
    
    o2 <- partitionGWAS(gwe,by="y",n=10,frac=0.8)
    expect_true(is(o2,"list"))
    expect_true(all(lengths(o2) == 0.8*(n1+n2)))
    
    o3 <- partitionGWAS(gwe,by="y",n=1,out="train")
    expect_true(is(o3,"GWASExperiment"))
    expect_equal(ncol(o3),0.5*ncol(gwe))
    
    o4 <- partitionGWAS(gwe,by="y",n=1,out="ttboth")
    expect_true(is.list(o4))
    expect_equal(names(o4),c("train","test"))
    
    # Partition a non-GWAS object - expect error
    expect_error(partitionGWAS(ph,n=10))
    expect_error(partitionGWAS(gwe,by="y",n=10,frac=1.1))
})

test_that("createSplit works",{
    y <- runif(100)
    y[y<0.5] <- 0
    y[y>=0.5] <- 1
    
    o1 <- createSplit(y,n=10)
    expect_equal(length(o1),10)
    expect_true(all(lengths(o1)==50))
    
    o2 <- createSplit(y,n=10,frac=0.8)
    expect_true(all(lengths(o2)==80))
    
    o3 <- createSplit(y,n=5,out="binary")
    expect_true(all(unlist(o3,use.names=FALSE) %in% c(0,1)))
    
    #TODO: Write o4 to test replacing
    
    expect_error(createSplit(y,n=10,frac=1.1))
    expect_error(createSplit(y,n=10,object="brick"))
})

# Test suite for .checkEFOFormat
test_that(".checkEFOFormat works",{
    # Properly formatted scalar input
    i1 <- "EFO_0001234"
    expect_true(length(.checkEFOFormat(i1)) == 1)
    
    # Properly formatted vector input
    i2 <- c("EFO_0001234","EFO_0005678")
    expect_true(length(.checkEFOFormat(i2)) == 2)
    
    # Malformed entry in vector input
    i3 <- c("EFO_0005937","EFO_0005937","MalformedPrefix_000000")
    expect_warning(o <- .checkEFOFormat(i3))
    expect_true(length(o) == 2)
    
    # All malformed entries - exception raised
    i4 <- c("EFO0001234","_0001234","0001234","TotallyMalformedId")
    expect_error(.checkEFOFormat(i4))
})

# Test suite for .checkPGSFormat
test_that(".checkPGSFormat works",{
    # Properly formatted scalar input
    i1 <- "PGS000001"
    expect_true(length(.checkPGSFormat(i1)) == 1)
    
    # Properly formatted vector input
    i2 <- c("PGS000001","PGS000298")
    expect_true(length(.checkPGSFormat(i2)) == 2)
    
    # Malformed entry in vector input
    i3 <- c("PGS000298","PGS_000298","PGS298","Malformed123")
    expect_warning(length(.checkPGSFormat(i3)) == 1)
    expect_true(length(suppressWarnings(.checkPGSFormat(i3))) == 1)
    
    # All malformed entries - exception raised
    i4 <- c("PGS123","_0001234","PGS","TotallyMalformedId")
    expect_error(.checkPGSFormat(i4))
})

# Test suite for .checkPGSFormat
test_that(".checkPMIDFormat works",{
    # Properly formatted scalar input
    i1 <- "33407065"
    expect_true(length(.checkPMIDFormat(i1)) == 1)
    
    # Properly formatted vector input
    i2 <- c("33407065","2369")
    expect_true(length(.checkPMIDFormat(i2)) == 2)
    
    # Malformed entry in vector input
    i3 <- c("33407065","2369_","PGS298","Malformed123")
    expect_warning(length(.checkPMIDFormat(i3)) == 1)
    expect_true(length(suppressWarnings(.checkPMIDFormat(i3))) == 1)
    
    # All malformed entries - exception raised
    i4 <- c("PGS123","_0001234","EFOPGS","TotallyMalformedId")
    expect_error(.checkPMIDFormat(i4))
})

test_that(".checkTextArgs works",{
    expect_silent(.checkTextArgs("test","test","test",FALSE))
    expect_silent(.checkTextArgs("test","test1",c("test1","test2"),TRUE))
    expect_error(.checkTextArgs("test","tes","test",FALSE))
})

test_that(".checkNumArgs works",{
    expect_silent(.checkNumArgs("test",4L,"integer",3L,"gt"))
    expect_silent(.checkNumArgs("test",4,"numeric",5,"lt"))
    expect_silent(.checkNumArgs("test",4,"numeric",c(3,4 ),"both"))
    expect_error(.checkNumArgs("test",3L,"integer",3L,"gt"))
})

test_that("Getting test files works",{
    input <- .gimmeTestFiles()
    expect_equal(length(input),3)
}

test_that(".isEmpty works",{
    x1 <- NULL
    expect_true(.isEmpty(x1))
    
    x2 <- NA
    expect_true(.isEmpty(x2))
    
    x3 <- ""
    expect_true(.isEmpty(x3))
    
    x4 <- integer(0)
    expect_true(.isEmpty(x4))
    
    x5 <- list()
    expect_true(.isEmpty(x5))
    
    x6 <- 42
    expect_false(.isEmpty(x6))
    
    x7 <- list(a=42)
    expect_false(.isEmpty(x7))
}

test_that(".validateBinaryForBinomial works",{
    x1 <- c(rep("A",5),rep("B",5))
    expect_warning(.validateBinaryForBinomial(x1))
    
    x2 <- c(rep(1,5),rep(2,5))
    expect_warning(.validateBinaryForBinomial(x2))
    
    x3 <- c(rep(1,5),rep(2,5),rep(3,5))
    expect_error(.validateBinaryForBinomial(x3))
    
    x4 <- c(rep(0,5),rep(1,5))
    expect_silent(x <- .validateBinaryForBinomial(x4))
    expect_identical(x,x4)
})

test_that("prismaVerbosity works",{
    expect_error(prismaVerbosity(1))
    expect_error(prismaVerbosity("invalid"))
    expect_silent(prismaVerbosity())
    expect_silent(prismaVerbosity("silent"))
    expect_equal(prismaVerbosity(),"silent")
})
