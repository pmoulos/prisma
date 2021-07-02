# Test suite for getGWASVariants
test_that("getGWASVariants works",{
    efoId <- c("EFO_0004566","EFO_0005937")
    efoTrait <- c("body weight gain","longitudinal BMI measurement") 
    
    test1 <- getGWASVariants(efoId=efoId,removeUnknownRisk=FALSE)
    expect_true(is.data.frame(test1))
    expect_true(any(is.na(test1$risk_allele)))
    
    test2 <- getGWASVariants(efoTrait=efoTrait,removeUnknownRisk=FALSE)
    expect_true(is.data.frame(test2))
    expect_true(any(is.na(test2$risk_allele)))
    
    test3 <- getGWASVariants(efoId=efoId,efoTrait=efoTrait,
        removeUnknownRisk=FALSE)
    expect_true(is.data.frame(test3))
    expect_true(any(is.na(test3$risk_allele)))
    
    test4 <- getGWASVariants(efoId=efoId)
    expect_true(is.data.frame(test4))
    expect_true(!any(is.na(test4$risk_allele)))
    expect_true(nrow(test4) <= nrow(test1))
})

# Test suite for .iterGWASVariantsAPICall
test_that(".iterGWASVariantsAPICall works",{
    # In test testing functions
    # Basic expected result informatin
    testBasic <- function(x,n) {
        expect_true(is.list(x))
        expect_true(length(x) == 2)
        expect_true(all(names(x) %in% n))
    }
    
    # Test expected data back
    testFound <- function(x) {
        expect_true(nrow(x[[1]]) > 0)
        expect_true(nrow(x[[2]]) > 0)
    }
    
    # Test expected data not found
    testNotFound <- function(x) {
        expect_true(nrow(x[[1]]) == 0)
        expect_true(nrow(x[[2]]) == 0)
    }
    
    inputId <- c("EFO_0004566","EFO_0005937")
    test <- .iterGWASVariantsAPICall(input=inputId,type="id")
    testBasic(test,inputId)
    testFound(test)
    
    inputTrait <- c("body weight gain","longitudinal BMI measurement")
    test <- .iterGWASVariantsAPICall(input=inputTrait,type="trait")
    testBasic(test,inputTrait)
    testFound(test)
    
    inputIdF <- c("EFO_0001234","EFO_0000000")
    test <- .iterGWASVariantsAPICall(input=inputIdF,type="id")
    testBasic(test,inputIdF)
    testNotFound(test)
    
    inputTraitF <- c("BODY gain","shortitudinal BMI measurement")
    test <- .iterGWASVariantsAPICall(input=inputTraitF,type="trait")
    testBasic(test,inputTraitF)
    testNotFound(test)
})

# Test suite for .GWASVariantsWorker
test_that(".GWASVariantsWorker works",{
    # Should return data rows if API call succesful
    inputId <- "EFO_0004566"
    test <- .GWASVariantsWorker(input=inputId,type="id")
    if (test$success)
        expect_true(nrow(test$data) > 0)
    else
        expect_true(test$data$input == "id")
    
    # Should return data rows if API call succesful
    inputTrait <- "body weight gain"
    test <- .GWASVariantsWorker(input=inputTrait,type="trait")
    if (test$success)
        expect_true(nrow(test$data) > 0)
    else
        expect_true(test$data$input == "trait")
    
    # Should return data rows if API call succesful - non-existent id
    inputIdF <- "EFO_0001234"
    test <- .GWASVariantsWorker(input=inputIdF,type="id")
    if (test$success)
        expect_true(nrow(test$data) == 0)
    else
        expect_true(test$data$input == "id")
    
    # Should return data rows if API call succesful - non-existent trait
    # as traits are case-sensitive
    inputTraitF <- "Body Weight Gain"
    test <- .GWASVariantsWorker(input=inputTraitF,type="trait")
    if (test$success)
        expect_true(nrow(test$data) == 0)
    else
        expect_true(test$data$input == "trait")
})

# Test suite for getGWASVariants
test_that("getPGSVariants works",{
    pgsId <- c("PGS000034","PGS000299")
    efoId <- c("EFO_0004340","EFO_0007788")
    pubmedId <- c("29212779","32527150") 
    
    test1 <- getPGSVariants(pgsId=pgsId)
    expect_true(is.data.frame(test1))
    #expect_true(any(is.na(test1$risk_allele)))
    
    test2 <- getPGSVariants(efoId=efoId)
    expect_true(is.data.frame(test2))
    #expect_true(any(is.na(test2$risk_allele)))
    
    test3 <- getPGSVariants(pgsId=pgsId,efoId=efoId)
    expect_true(is.data.frame(test3))
    #expect_true(any(is.na(test3$risk_allele)))
    
    test4 <- getPGSVariants(pubmedId=pubmedId)
    expect_true(is.data.frame(test4))
    #expect_true(!any(is.na(test4$risk_allele)))
    expect_true(nrow(test4) <= nrow(test1))
})

# Test suite for .iterGWASVariantsAPICall
test_that(".iterPGSVariantsAPICall works",{
    # In test testing functions
    # Basic expected result informatin
    testBasic <- function(x,n) {
        expect_true(is.list(x))
        expect_true(length(x) == 2)
        expect_true(all(names(x) %in% n))
    }
    
    # Test expected data back
    testFound <- function(x) {
        expect_true(nrow(x[[1]]) > 0)
        expect_true(nrow(x[[2]]) > 0)
    }
    
    # Test expected data not found
    testNotFound <- function(x) {
        expect_true(nrow(x[[1]]) == 0)
        expect_true(nrow(x[[2]]) == 0)
    }
    
    inputId <- c("EFO_0004566","EFO_0005937")
    test <- .iterGWASVariantsAPICall(input=inputId,type="id")
    testBasic(test,inputId)
    testFound(test)
    
    inputTrait <- c("body weight gain","longitudinal BMI measurement")
    test <- .iterGWASVariantsAPICall(input=inputTrait,type="trait")
    testBasic(test,inputTrait)
    testFound(test)
    
    inputIdF <- c("EFO_0001234","EFO_0000000")
    test <- .iterGWASVariantsAPICall(input=inputIdF,type="id")
    testBasic(test,inputIdF)
    testNotFound(test)
    
    inputTraitF <- c("BODY gain","shortitudinal BMI measurement")
    test <- .iterGWASVariantsAPICall(input=inputTraitF,type="trait")
    testBasic(test,inputTraitF)
    testNotFound(test)
})

# Test suite for .GWASVariantsWorker
test_that(".PGSVariantsWorker works",{
    # Should return data rows if API call succesful
    inputId <- "EFO_0004566"
    test <- .GWASVariantsWorker(input=inputId,type="id")
    if (test$success)
        expect_true(nrow(test$data) > 0)
    else
        expect_true(test$data$input == "id")
    
    # Should return data rows if API call succesful
    inputTrait <- "body weight gain"
    test <- .GWASVariantsWorker(input=inputTrait,type="trait")
    if (test$success)
        expect_true(nrow(test$data) > 0)
    else
        expect_true(test$data$input == "trait")
    
    # Should return data rows if API call succesful - non-existent id
    inputIdF <- "EFO_0001234"
    test <- .GWASVariantsWorker(input=inputIdF,type="id")
    if (test$success)
        expect_true(nrow(test$data) == 0)
    else
        expect_true(test$data$input == "id")
    
    # Should return data rows if API call succesful - non-existent trait
    # as traits are case-sensitive
    inputTraitF <- "Body Weight Gain"
    test <- .GWASVariantsWorker(input=inputTraitF,type="trait")
    if (test$success)
        expect_true(nrow(test$data) == 0)
    else
        expect_true(test$data$input == "trait")
})

# Test suite for .emptyVariantsDf
test_that(".emptyVariantsDf works",{
    expect_true(nrow(.emptyVariantsDf()) == 0)
})
