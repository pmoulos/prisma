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
    expect_true(length(.checkEFOFormat(i3)) == 2)
    
    # All malformed entries - exception raised
    i4 <- c("EFO0001234","_0001234","0001234","TotallyMalformedId")
    expect_true(.checkEFOFormat(i4))
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

# Test suite for .emptyVariantsDf
test_that(".emptyVariantsDf works",{
    expect_true(nrow(.emptyVariantsDf()) == 0)
})

test_that("Get/Set API base works",{
    # Set API base
    mock <- "https://www.example.com/api"
    setAPIBase(mock)
    expect_true(getAPIBase() == mock)    
    
    # Get API base
    setAPIBase() # revert first
    base <- getAPIBase()
    expect_true(base == .defaultUrlBase())
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
