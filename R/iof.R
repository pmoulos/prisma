# input may be a prefix of bim, bed, fam or a list with each
importGWAS <- function(input,phenos=NULL,backend=c("snpStats","bigsnpr"),
    selection=NULL,genome=NA_character_,
    gdsfile=ifelse(backend=="snpStats",tempfile(),NA)) {
    backend <- backend[1]
    
    if (backend == "snpStats") {
        if (!requireNamespace("snpStats"))
            stop("Bioconductor package snpStats is required!")
        if (!requireNamespace("SNPRelate"))
            stop("Bioconductor package SNPRelate is required!")
        if (is.na(gdsfile) || !is.character(gdsfile)) {
            warning("A valid path to store GDS file for filtering is required!",
                " Assuming default (tempfile())...")
            gdsfile <- tempfile()
        }
    }
    if (backend == "bigsnpr" && !requireNamespace("snpStats"))
        stop("Bioconductor package bigsnpr is required!")
    
    if (!is.null(phenos)) {
        if (!is.character(phenos) && !is.data.frame(phenos))
            stop("The phenotypes argument (pheno) must be a data.frame or an ",
                "existing file!")
        if (is.character(phenos) && !file.exists(phenos))
            stop("When character, the phenotypes argument (pheno) must be an ",
                "existing file!")
        if (file.exists(phenos))
            phenos <- read.delim(phenos)
    }
    
    selection <- .checkSelection(selection)
    
    if (is.character(input)) {
        if (!dir.exists(dirname(input)))
            stop("The input directory containing PLINK files does not exist!")
        input <- list(bed=paste0(input,".bed"),bim=paste0(input,".bim"),
            fam=paste0(input,".fam"))
    }
    plExist <- vapply(input,file.exists,logical(1))
    if (!all(plExist))
        stop("Input PLINK file(s) ",paste(input[!plExist],collapse=", "),
            " do not exist!")
    
    if (backend == "snpStats") {
        disp("Reading PLINK files with snpStats framework")
        snpObj <- read.plink(input$bed,input$bim,input$fam,
            select.subjects=selection$samples,select.snps=selection$snps)
        disp("Reading PLINK files with SNPRelate framework and storing ",
            "output to ",gdsfile)
        snpgdsBED2GDS(input$bed,input$fam,input$bim,gdsfile,family=TRUE)
        # We also need to read a GDS file for LD and IBD filtering
        return(GWASExperiment(
            genotypes=t(snpObj$genotypes),
            features=snpObj$map,
            samples=snpObj$fam,
            phenotypes=phenos,
            metadata=list(
                genome=genome,
                backend=backend,
                filters=.initFilterInfo(),
                gdsfile=gdsfile
            )
        ))
    }
    else if (backend == "bigsnpr") {
        # TODO: Add bigsnpr options
        #snpObj <- ... snp_readSth
    }
        
}

filterGWAS <- function(obj,filters=getDefaults("filters"),imputeMissing=TRUE,
    rc=NULL) {
    filters <- .checkFilters(filters)
    
    # Later input type may be more native, e.g. read from CSV files
    if (is(assay(obj,1),"SnpMatrix"))
        return(.filterWithSnpStats(obj,filters,imputeMissing,rc))
    else if (is(assay(obj,1),"bigsnp"))
        return(.filterWithBigSnpr(x,filters,imputeMissing))
}

imputeGWAS <- function(obj,mode=c("single","split"),rc=NULL) {
    mode <- mode[1]
    .checkTextArgs("Imputation mode (mode)",mode,c("single","split"),
        multiarg=FALSE)
    
    m <- metadata(obj)
    if (m$backend == "snpStats") {
        disp("\nImputing missing values with snpStats rules and scrime kNN")
        return(.internalImputeWithSnpStats(obj,rc=rc))
    }
    else
        return(obj) # For now
}

.checkSelection <- function(s) {
    if (is.null(s))
        return(s)
    
    if (!is.list(s))
        stop("The selection argument must be a list!")
    if (!(all(names(s) %in% c("samples","snps"))))
        stop("The selection argument must be a named list with members ",
            "'samples' and 'snps'!")
            
    problemFound <- vapply(s,function(x) {
        if (!is.null(x)) {
            if (is.numeric(x) && all(x>0))
                return(FALSE)
            else
                return(TRUE)
        }
        else
            return(FALSE)
    },logical(1))
            
    if (any(problemFound))
        stop("All selection list members must be non-negative integers ",
            "denoting SNP matrix rows\nand/or columns! ",
            paste(names(s)[problemFound],collapse=", ")," are not...")
    
    return(s)
}

writePlink <- function(obj,outBase,salvage=FALSE) {
    if (!is(obj,"GWASExperiment"))
        stop("The input object (obj) must be a GWASExperiment object!")
    if (missing(outBase))
        stop("A path and filename base for the resulting PLINK files must be ",
            "provided!")
    
    map <- gfeatures(obj)
    fam <- gsamples(obj)
    gen <- t(genotypes(obj))
    
    if (salvage) {
        if (is.na(genome(obj)))
            warning("Cannot salvage SNP locations if genome version/build is ",
                "not specified! Skipping...")
        else {
            noChrRsInd <- which(is.na(map$chromosome) & 
                grepl("^rs",rownames(map),perl=TRUE))
            noChrRs <- rownames(map)[noChrRsInd]
            names(noChrRsInd) <- noChrRs
            rsInfo <- rsLocsFromEnsembl(noChrRs,genome(obj),canonical=TRUE)
        }
            
        map$chromosome[noChrRsInd[rsInfo$refsnp_id]] <- rsInfo$chr_name
        map$position[noChrRsInd[rsInfo$refsnp_id]] <- rsInfo$chrom_start
    }
    
    # Preflight...
    map$chromosome <- as.character(map$chromosome)
    if (any(map$chromosome %in% c("23","24","25","26"))) {
        map$chromosome[map$chromosome=="23"] <- "X"
        map$chromosome[map$chromosome=="24"] <- "Y"
        chromosome[map$chromosome=="25"] <- "XY"
        map$chromosome[map$chromosome=="26"] <- "MT"
    }
    
    # Filter out those with remaining NA chromosome, position (controls?)
    na <- which(is.na(map$chromosome))

    # Write the PLINK files that we will work with!
    write.plink(
        file.base=outBase,
        snps=gen[,-na],
        pedigree=fam$pedigree,
        id=rownames(fam),
        father=fam$father,
        mother=fam$mother,
        sex=fam$sex,
        phenotype=fam$affected,
        chromosome=map$chromosome[-na],
        position=map$position[-na],
        allele.1=map$allele.1[-na],
        allele.2=map$allele.2[-na]
    )
}
