test_that("downloadSnptest works",{
    # Simple dummy data test
})

test_that(".tryNumerize works",{
    expect_true(is.numeric(.tryNumerize("2.3")))
    expect_true(is.numeric(.tryNumerize("2.3.4")))
    expect_true(is.numeric(.tryNumerize("v2.3.4-rc123@45")))
    expect_true(is.numeric(.tryNumerize("1980-04-25_14:45:56")))
    expect_true(is.null(.tryNumerize("rcABS")))
    expect_true(is.null(.tryNumerize(NULL)))
})
