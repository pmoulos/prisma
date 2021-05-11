# Pseudodata generating function
.makeThingData <- function(error=c("none","snp_class","feature_size",
    "feature_name","sample_size","sample_name","pheno_size","pheno_name")) {
    if (missing(error))
        error <- "none"
    
    # FAM table
    proband <- c(1,4,7,10,13,16)
    father <- mother <- affected <- sex <- rep(NA,18)
    father[proband] <- paste0("sample_",proband+1)
    mother[proband] <- paste0("sample_",proband+2)
    affected[proband] <- 2
    sex[proband] <- sample(c(1,2),6,replace=TRUE)
    sex[proband+1] <- 1
    sex[proband+2] <- 2
    fam <- data.frame(
        pedigree=rep(seq_len(6),each=3),
        member=paste0("sample_",seq_len(18)),
        father=father,
        mother=mother,
        sex=sex,
        affected=affected,
        row.names=paste0("sample_",seq_len(18))
    )
    
    # MAP table
    nuc <- c("A","C","G","T")
    allele_a <- sample(nuc,100,replace=TRUE)
    allele_b <- sample(nuc,100,replace=TRUE)
    while(any(allele_a==allele_b)) {
        s <- which(allele_a==allele_b)
        replace=ifelse(length(s)<4,TRUE,FALSE)
        allele_a[s] <- sample(nuc,length(s),replace=TRUE)
    }
    map <- data.frame(
        chromosome=sample(seq_len(18),100,replace=TRUE),
        snp_name=paste0("snp_",seq_len(100)),
        position=sample(1000:1000000,100),
        allele_a=allele_a,
        allele_b=allele_b,
        row.names=paste0("snp_",seq_len(100))
    )

    # Phenotypes
    pseudopheno <- data.frame(
        case_control=sample(c(0,1),nrow(fam),replace=TRUE),
        other_pheno=sample(c("drug","nodrug"),nrow(fam),replace=TRUE),
        yap=sample(c("normal","demipsek","psek"),nrow(fam),replace=TRUE),
        cont=round(runif(nrow(fam)),3),
        row.names=rownames(fam)
    )
    
    # Genotypes
    snp <- matrix(as.raw(sample(c(0,1,2),1800,replace=TRUE)),18,100)
    rownames(snp) <- rownames(fam)
    colnames(snp) <- rownames(map)
    
    if ("none" %in% error)
        return(list(
            snp=SnpMatrix(t(snp)),
            sample=fam,
            feature=map,
            pheno=pseudopheno
        ))
    else {
        if ("snp_class" %in% error)
            snp <- t(snp)
        else
            snp <- SnpMatrix(t(snp))
        
        if ("feature_size" %in% error)
            map <- map[-sample(nrow(map),1),,drop=FALSE]
        if ("feature_name" %in% error)
            rownames(map)[sample(nrow(map),1)] <- "bad_name"
        if ("sample_size" %in% error)
            fam <- fam[-sample(nrow(fam),1),,drop=FALSE]
        if ("sample_name" %in% error)
            rownames(fam)[sample(nrow(fam),1)] <- "bad_name"
        if ("pheno_size" %in% error)
            pseudopheno <- pseudopheno[-sample(nrow(pseudopheno),1),,drop=FALSE]
        if ("pheno_name" %in% error)
            rownames(pseudopheno)[sample(nrow(pseudopheno),1)] <- "bad_name"
        
        return(list(
            snp=snp,
            sample=fam,
            feature=map,
            pheno=pseudopheno
        ))
    }
}


test_that("GWASExperiment definition works",{
    expect_error(validObject(.GWASExperiment())) # should not be used
    expect_silent(GWASExperiment())
    expect_true(validObject(GWASExperiment()))
})

test_that("GWASExperiment validation works - warnings",{
    set.seed(42)
    
    thing <- .makeThingData()
    
    # All should be ok
    expect_silent(GWASthing <- GWASExperiment(
        genotypes=thing$snp,
        features=thing$feature,
        samples=thing$sample,
        phenotypes=thing$pheno,
        metadata=list(genome=NA_character_,backend=NA_character_,
            filters=.initFilterInfo()
        )
    ))
    expect_true(validObject(GWASthing))
    
    # Warnings should be issued
    expect_warning(GWASExperiment(
        genotypes=thing$snp,
        features=thing$feature,
        samples=thing$sample,
        phenotypes=thing$pheno,
        metadata=list()
    ))
    expect_warning(GWASExperiment(
        genotypes=thing$snp,
        features=thing$feature,
        samples=thing$sample,
        phenotypes=thing$pheno,
        metadata=list(genome="hg19")
    ))
    expect_warning(GWASExperiment(
        genotypes=thing$snp,
        features=thing$feature,
        samples=thing$sample,
        phenotypes=thing$pheno,
        metadata=list(genome="hg19")
    ))    
})

test_that("GWASExperiment validation works - single errors",{
    set.seed(42)
    
    # Filter error
    expect_error(GWASExperiment(
        metadata=list(genome=NA_character_,backend=NA_character_,
            filters=data.frame(a=NA,b=NA))
    ))
    
    # Single errors
    thing <- .makeThingData(error="snp_class")
    expect_error(GWASExperiment(
        genotypes=thing$snp,
        features=thing$feature,
        samples=thing$sample,
        phenotypes=thing$pheno
    ))
    expect_error(GWASExperiment(
        genotypes=thing$snp,
        features=thing$feature,
        samples=thing$sample,
        phenotypes=thing$pheno
    ))
    
    thing <- .makeThingData(error="feature_size")
    expect_error(GWASExperiment(
        genotypes=thing$snp,
        features=thing$feature,
        samples=thing$sample,
        phenotypes=thing$pheno
    ))
    
    # Should not throw error as withDimnames=TRUE and assay rownames are always
    # replacing rowData rownames
    thing <- .makeThingData(error="feature_name")
    expect_silent(GWASExperiment(
        genotypes=thing$snp,
        features=thing$feature,
        samples=thing$sample,
        phenotypes=thing$pheno
    ))
    
    thing <- .makeThingData(error="sample_size")
    expect_error(GWASExperiment(
        genotypes=thing$snp,
        features=thing$feature,
        samples=thing$sample,
        phenotypes=thing$pheno
    ))
    
    # Should not throw error as withDimnames=TRUE and assay colnames are always
    # replacing colData rownames
    thing <- .makeThingData(error="sample_name")
    expect_silent(GWASExperiment(
        genotypes=thing$snp,
        features=thing$feature,
        samples=thing$sample,
        phenotypes=thing$pheno
    ))
    
    thing <- .makeThingData(error="pheno_size")
    expect_error(GWASExperiment(
        genotypes=thing$snp,
        features=thing$feature,
        samples=thing$sample,
        phenotypes=thing$pheno
    ))
    
    # Should not throw error as withDimnames=TRUE and assay colnames are always
    # replacing phenotypes rownames
    thing <- .makeThingData(error="pheno_name")
    expect_silent(GWASExperiment(
        genotypes=thing$snp,
        features=thing$feature,
        samples=thing$sample,
        phenotypes=thing$pheno
    ))
})

test_that("GWASExperiment validation works - multiple errors",{
    set.seed(42)
    
    # Multiple errors
    thing <- .makeThingData(error=c("feature_size","sample_size"))
    expect_error(GWASExperiment(
        genotypes=thing$snp,
        features=thing$feature,
        samples=thing$sample,
        phenotypes=thing$pheno
    ))
    
    # Should not throw error as withDimnames=TRUE and assay dimnames are always
    # replacing annotation data dimnames
    thing <- .makeThingData(error=c("feature_name","sample_name"))
    expect_silent(GWASExperiment(
        genotypes=thing$snp,
        features=thing$feature,
        samples=thing$sample,
        phenotypes=thing$pheno
    ))
})

test_that("GWASExperiment getters work",{
    set.seed(42)
   
    thing <- .makeThingData()
    GWASthing <- GWASExperiment(
        genotypes=thing$snp,
        features=thing$feature,
        samples=thing$sample,
        phenotypes=thing$pheno
    )
    
    expect_identical(rownames(assays(GWASthing)[[1]]),rownames(thing$snp))
    expect_identical(colnames(assays(GWASthing)[[1]]),colnames(thing$snp))
    
    expect_identical(rownames(rowData(GWASthing)),rownames(thing$feature))
    expect_identical(colnames(rowData(GWASthing)),colnames(thing$feature))
    expect_identical(rowData(GWASthing),DataFrame(thing$feature))
    
    expect_identical(rownames(colData(GWASthing)),rownames(thing$sample))
    expect_identical(colnames(colData(GWASthing)),colnames(thing$sample))
    expect_identical(colData(GWASthing),DataFrame(thing$sample))
    
    expect_identical(rownames(phenotypes(GWASthing)),rownames(thing$pheno))
    expect_identical(colnames(phenotypes(GWASthing)),colnames(thing$pheno))
    expect_identical(phenotypes(GWASthing),thing$pheno)
    
    i <- seq_len(10)
    j <- seq_len(5)
    Gslice <- GWASthing[i,j]
    expect_true(nrow(Gslice)==length(i))
    expect_true(ncol(Gslice)==length(j))
    expect_identical(rowData(Gslice),DataFrame(thing$feature[i,,drop=FALSE]))
    expect_identical(colData(Gslice),DataFrame(thing$sample[j,,drop=FALSE]))
    expect_identical(phenotypes(Gslice),thing$pheno[j,,drop=FALSE])
})

test_that("GWASExperiment setters work",{
    set.seed(42)
    thingOld <- .makeThingData()
    GWASthing <- GWASExperiment(
        genotypes=thingOld$snp,
        features=thingOld$feature,
        samples=thingOld$sample,
        phenotypes=thingOld$pheno
    )
    
    set.seed(43)
    thingNew <- .makeThingData()
    
    gfeatures(GWASthing) <- thingNew$feature
    expect_identical(rowData(GWASthing),DataFrame(thingNew$feature))
    
    gsamples(GWASthing) <- thingNew$sample
    expect_identical(colData(GWASthing),DataFrame(thingNew$sample))
    
    phenotypes(GWASthing) <- thingNew$pheno
    expect_identical(phenotypes(GWASthing),thingNew$pheno)
    
    genome(GWASthing) <- "hg19"
    expect_identical(genome(GWASthing),"hg19")
})

test_that("GWASExperiment setters work",{
    set.seed(42)
    thingOld <- .makeThingData()
    GWASthing <- GWASExperiment(
        genotypes=thingOld$snp,
        features=thingOld$feature,
        samples=thingOld$sample,
        phenotypes=thingOld$pheno
    )
    
    set.seed(43)
    thingNew <- .makeThingData()
    
    gfeatures(GWASthing) <- thingNew$feature
    expect_identical(rowData(GWASthing),DataFrame(thingNew$feature))
    
    gsamples(GWASthing) <- thingNew$sample
    expect_identical(colData(GWASthing),DataFrame(thingNew$sample))
    
    phenotypes(GWASthing) <- thingNew$pheno
    expect_identical(phenotypes(GWASthing),thingNew$pheno)
    
    genome(GWASthing) <- "hg19"
    expect_identical(genome(GWASthing),"hg19")
})

test_that("GWASExperiment cbind works",{
    set.seed(42)
    thing1 <- .makeThingData()
    GWASthing1 <- GWASExperiment(
        genotypes=thing1$snp,
        features=thing1$feature,
        samples=thing1$sample,
        phenotypes=thing1$pheno
    )
    
    set.seed(43)
    thing2 <- .makeThingData()
    colnames(thing2$snp) <- rownames(thing2$sample) <- 
        rownames(thing2$pheno) <- paste0(colnames(thing1$snp),"_1")
    GWASthing2 <- GWASExperiment(
        genotypes=thing2$snp,
        features=thing1$feature,
        samples=thing2$sample,
        phenotypes=thing2$pheno
    )
    GWASthingF <- GWASExperiment(
        genotypes=thing2$snp,
        features=thing2$feature,
        samples=thing2$sample,
        phenotypes=thing2$pheno
    )
    
    expect_identical(rowData(GWASthing1),rowData(GWASthing2))
    expect_false(identical(colData(GWASthing1),colData(GWASthing2)))
    expect_false(identical(phenotypes(GWASthing1),phenotypes(GWASthing2)))
    
    expect_error(GWASthing <- cbind(GWASthing1,GWASthingF))
    expect_silent(GWASthing <- cbind(GWASthing1,GWASthing2))
    expect_equal(length(GWASthing),length(GWASthing1))
    expect_equal(ncol(GWASthing),ncol(GWASthing1) + ncol(GWASthing2))
})

test_that("GWASExperiment rbind works",{
    set.seed(42)
    thing1 <- .makeThingData()
    GWASthing1 <- GWASExperiment(
        genotypes=thing1$snp,
        features=thing1$feature,
        samples=thing1$sample,
        phenotypes=thing1$pheno
    )
    
    set.seed(43)
    thing2 <- .makeThingData()
    rownames(thing2$snp) <- rownames(thing2$feature) <- 
        paste0(rownames(thing1$snp),"_1")
    GWASthing2 <- GWASExperiment(
        genotypes=thing2$snp,
        features=thing2$feature,
        samples=thing1$sample,
        phenotypes=thing1$pheno
    )
    GWASthingF <- GWASExperiment(
        genotypes=thing2$snp,
        features=thing2$feature,
        samples=thing2$sample,
        phenotypes=thing2$pheno
    )
    
    expect_identical(colData(GWASthing1),colData(GWASthing2))
    expect_false(identical(rowData(GWASthing1),rowData(GWASthing2)))
    expect_identical(phenotypes(GWASthing1),phenotypes(GWASthing2))
    
    expect_error(GWASthing <- rbind(GWASthing1,GWASthingF))
    expect_silent(GWASthing <- rbind(GWASthing1,GWASthing2))
    expect_equal(length(GWASthing),length(GWASthing1) + length(GWASthing2))
    expect_equal(ncol(GWASthing),ncol(GWASthing1))
})
