test_that(".gwaGlmWorker works",{
    # Simple dummy data test
})

test_that("gwaGlm works",{
    # Maybe with a real test dataset or a test "thing"
})

test_that(".prepareGenotypesForBlup works",{
    # Maybe with a real test dataset or a test "thing"
})

test_that(".estimateNPCinBlup works",{
    # How?
})

test_that("gwaBlup works",{
    # Maybe with a real test dataset or a test "thing"
})

test_that("gwaStatgen works",{
    # Maybe with a real test dataset or a test "thing"
})

test_that(".validateResponseAndCovariates",{
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
