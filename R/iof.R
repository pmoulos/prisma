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
    
    selection <- .checkSelection(selection)
    
    if (is.character(input)) {
        if (!dir.exists(basename(input)))
            stop("The input directory containing PLINK files does not exist!")
        input <- list(bed=paste0(input,".bed"),bim=paste0(input,".bim"),
            bed=paste0(input,".fam"))
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

