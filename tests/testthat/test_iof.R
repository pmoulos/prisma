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
