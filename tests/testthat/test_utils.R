test_that("Get/Set API base",{
    # Set API base
    mock <- "https://www.example.com/api"
    setAPIBase(mock)
    expect_true(getAPIBase() == mock)    
    
    # Get API base
    setAPIBase() # revert first
    base <- getAPIBase()
    expect_true(base == .defaultUrlBase())
})
