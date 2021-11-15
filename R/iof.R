# Update:
# I checked manually the below... It is not true and genotypes are not
# reversed... There may be differences on how each software counts...
# Until we figure what is going on, we reverse the manual PRS
#
# !!! SOS !!! snpStats read.plink **REVERSES** genotypes!
# From the official PLINK documentation
# https://www.cog-genomics.org/plink/1.9/formats#bed
#
# Allele 1 (usually minor), 'X' if absent
# Allele 2 (usually major), 'X' if absent
#
# 00    Homozygous for first allele in .bim file
# 01    Missing genotype
# 10    Heterozygous
# 11    Homozygous for second allele in .bim file
#
# Therefore it should be (with as.numeric conversion)
# 0 = homozygous minor (risk/effect)
# 1 = heterozygous 
# 2 = homozygous major (non-risk)
#
# In map file read by snpStats we have (since it's like that in default PLINK)
# allele.1 - minor
# allele.2 - major
#
# On the contrary it is
#
# 0 = homozygous major (non-risk)
# 1 = heterozygous 
# 2 = homozygous minor (risk/effect)
#
# The same notation is followed by the gaston R package and possibly also by
# bigsnpr (SOS)
#

# input may be a prefix of bim, bed, fam or a list with each
importGWAS <- function(input,phenos=NULL,backend=c("snpStats","bigsnpr"),
    selection=NULL,genome=NA_character_,alleleOrder=c("plink","reverse"),
    gdsfile=ifelse(backend=="snpStats",tempfile(),NA),gdsOverwrite=TRUE) {
    
    backend <- backend[1]
    alleleOrder <- alleleOrder[1]
    if (!(backend %in% c("snpStats","bigsnpr")))
        stop("backend should be one of \"snpStats\" or \"bigsnpr\"")
    .checkTextArgs("Allelic order in input (alleleOrder)",alleleOrder,
        c("plink","reverse"),multiarg=FALSE)
    
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
        if (is.character(phenos) && !file.exists(phenos))
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
        if (file.exists(gdsfile) && !gdsOverwrite)
            disp("Skipping GDS file creation as ",gdsfile," already exists...")
        else {
            disp("Reading PLINK files with SNPRelate framework and storing ",
                "output to ",gdsfile)
            # We also need to read a GDS file for LD and IBD filtering
            snpgdsBED2GDS(input$bed,input$fam,input$bim,gdsfile,family=TRUE)
        }
        ## ***SWITCH ALLELES*** if PLINK default
        #if (alleleOrder == "plink") {
        #    disp("Switching SnpMatrix alleles to comply with default PLINK ",
        #        "behaviour")
        #    snpObj$genotypes <- switch.alleles(snpObj$genotypes,
        #        seq_len(ncol(snpObj$genotypes)))
        #}
        
        return(GWASExperiment(
            genotypes=t(snpObj$genotypes),
            features=snpObj$map,
            samples=snpObj$fam,
            phenotypes=phenos,
            metadata=list(
                genome=genome,
                backend=backend,
                filters=.initFilterInfo(),
                gdsfile=gdsfile,
                alleleOrder=alleleOrder
            )
        ))
    }
    else if (backend == "bigsnpr") {
        # TODO: Add bigsnpr options
        #snpObj <- ... snp_readSth
    }
        
}

filterGWAS <- function(obj,filters=getDefaults("filters"),imputeMissing=TRUE,
    imputeMode=c("single","split"),rc=NULL) {
    filters <- .checkFilters(filters)
    imputeMode <- imputeMode[1]
    
    # Later input type may be more native, e.g. read from CSV files
    if (is(assay(obj,1),"SnpMatrix"))
        return(.filterWithSnpStats(obj,filters,imputeMissing,imputeMode,rc))
    else if (is(assay(obj,1),"bigsnp"))
        return(.filterWithBigSnpr(x,filters,imputeMissing))
}

imputeGWAS <- function(obj,mode=c("single","split"),failed=c("scrime","none"),
    rc=NULL) {
    mode <- mode[1]
    failed <- failed[1]
    .checkTextArgs("Imputation mode (mode)",mode,c("single","split"),
        multiarg=FALSE)
    .checkTextArgs("Imputation failure action (failed)",failed,
        c("scrime","none"),multiarg=FALSE)
    
    m <- metadata(obj)
    if (m$backend == "snpStats") {
        ex <- ifelse(failed=="scrime"," and scrime kNN","")
        disp("\nImputing missing values with snpStats rules",ex)
        return(.internalImputeWithSnpStats(obj,mode,failed,rc=rc))
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

writePlink <- function(obj,pheno=NULL,outBase=NULL,salvage=FALSE,
    reverse=FALSE) {
    if (!is(obj,"GWASExperiment"))
        stop("The input object (obj) must be a GWASExperiment object!")
    if (is.null(outBase))
        stop("A path and filename base for the resulting PLINK files must be ",
            "provided!")
    
    map <- gfeatures(obj)
    fam <- gsamples(obj)
    gen <- t(genotypes(obj))
    
    # If the input object comes with default PLINK allele order, then they were
    # reversed upon importing and also have the alleleOrder metadata field.
    if (reverse)
        gen <- switch.alleles(gen,seq_len(ncol(gen)))
    
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
    na <- is.na(map$chromosome)
    
    # Finally, attach the desired phenotype (if provided, otherwise default)
    ph <- phenotypes(obj)
    if (is.null(ph) || is.null(pheno) || (is(ph,"DataFrame") && nrow(ph) == 0)
        || (is.data.frame(ph) && nrow(ph) == 0))
        phenoVec <- fam$affected
    else {
        if (is.numeric(pheno)) {
            if (pheno > ncol(ph)) {
                warning("The main phenotype index ",pheno," is larger than ",
                    "the object's available phenotypes! Exporting default...")
                phenoVec <- fam$affected
            }
            else {
                pheno <- colnames(ph)[pheno]
                phenoVec <- ph[,pheno]
            }
        }
        else if (is.character(pheno)) {
            if (!(pheno %in% colnames(ph))) {
                warning("The main phenotype ",pheno," was not found in the ",
                    "object's phenotypes! Exporting default...")
                phenoVec <- fam$affected
            }
            else
                phenoVec <- ph[,pheno]
        }
    }

    # Write the PLINK files that we will work with!
    write.plink(
        file.base=outBase,
        snps=gen[,!na],
        pedigree=fam$pedigree,
        id=rownames(fam),
        father=fam$father,
        mother=fam$mother,
        sex=fam$sex,
        phenotype=phenoVec,
        chromosome=map$chromosome[!na],
        position=map$position[!na],
        allele.1=map$allele.1[!na],
        allele.2=map$allele.2[!na]
    )
}
