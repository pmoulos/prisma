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

test_that("writePlink from GWAS experiment works",{
    dat <- makeSimData(nsnp=100,nsam=10,nphe=1,csnp=10)
    dat$pheno$test <- sample(c(1L,2L),10,replace=TRUE)

    gwe <- GWASExperiment(
        genotypes=dat$snp,
        features=dat$feature,
        samples=dat$sample,
        phenotypes=dat$pheno
    )
    
    # Tests general, reverse, with phenotype, per chr
    
    # General write
    outBase1 <- tempfile()
    writePlink(gwe,outBase=outBase1)
    expect_true(file.exists(paste0(outBase1,".bed")))
    expect_true(file.exists(paste0(outBase1,".bim")))
    expect_true(file.exists(paste0(outBase1,".fam")))
    
    # With reverse alleles
    outBase2 <- tempfile()
    writePlink(gwe,outBase=outBase2,reverse=TRUE)
    gwe1 <- importGWAS(outBase2,backend="snpStats",writeGds=FALSE)
    g <- as(genotypes(gwe1),"numeric") + as(genotypes(gwe),"numeric")
    expect_true(all(g==2))
    
    # With phenotype
    outBase3 <- tempfile()
    writePlink(gwe,outBase=outBase3,pheno="test")
    gwe2 <- importGWAS(outBase3,backend="snpStats",writeGds=FALSE)
    g <- gsamples(gwe2)
    expect_identical(dat$pheno$test,g$affected)
    
    # Per chromosome
    outBase4 <- tempfile()
    writePlink(gwe,outBase=outBase4,pheno="test",perChr=TRUE)
    cc <- unique(dat$feature$chromosome)
    expect_true(file.exists(paste0(outBase4,"_chr",cc[1],".bed")))
    expect_true(file.exists(paste0(outBase4,"_chr",cc[2],".bed")))
})

test_that("writeGdsFile from GWAS experiment works",{
    dat <- makeSimData(nsnp=100,nsam=10,nphe=1,csnp=10)
    gwe <- GWASExperiment(
        genotypes=dat$snp,
        features=dat$feature,
        samples=dat$sample,
        phenotypes=dat$pheno
    )
    gdsgwe <- writeGdsFile(gwe)
    expect_true(file.exists(gdsgwe))
})
