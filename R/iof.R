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
    writeGds=TRUE,gdsfile=ifelse(backend=="snpStats",tempfile(),NA),
    gdsOverwrite=TRUE) {
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
        if (writeGds && (is.null(gdsfile) || !is.character(gdsfile))) {
            warning("A valid path to store GDS file for filtering is required!",
                " Assuming default...",immediate.=TRUE)
            gdsfile <- file.path(dirname(input[[1]]),
                sub("(.bed|.bim|.fam)$","",basename(input[[1]]),
                    ignore.case=TRUE))
            gdsfile <- paste0(gdsfile,".gds")
        }
    }
    if (backend == "bigsnpr" && !requireNamespace("bigsnpr"))
        stop("Bioconductor package bigsnpr is required!")
    
    if (!is.null(phenos)) {
        if (!is.character(phenos) && !is.data.frame(phenos)
            && !is(phenos,"DataFrame"))
            stop("The phenotypes argument (pheno) must be a data.frame or a ",
                "DataFrame or an existing file!")
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
        if (writeGds) {
            if (file.exists(gdsfile) && !gdsOverwrite)
                disp("Skipping GDS file creation as ",gdsfile,
                    " already exists...")
            else
                .writeGdsFile(input,snpObj,gdsfile,selection)
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

writeGdsFile <- function(obj,gdsfile=NULL) {
    if (!is(obj,"GWASExperiment"))
        stop("Input object must be a GWASExperiment")
    if (is.null(gdsfile))
        gdsfile <- tempfile()
    
    map <- as.data.frame(gfeatures(obj))
    fam <- as.data.frame(gsamples(obj))
    geno <- t(genotypes(obj))
    tmpBase <- tempfile()
    
    write.plink(
        file.base=tmpBase,
        snps=geno,
        pedigree=fam$pedigree,
        id=rownames(fam),
        father=fam$father,
        mother=fam$mother,
        sex=fam$sex,
        phenotype=fam$affected,
        chromosome=map$chromosome,
        position=map$position,
        allele.1=map$allele.1,
        allele.2=map$allele.2
    )
    
    snpgdsBED2GDS(paste0(tmpBase,".bed"),paste0(tmpBase,".fam"),
        paste0(tmpBase,".bim"),gdsfile,family=TRUE)
        
    unlink(paste0(tmpBase,".bed"),recursive=TRUE,force=TRUE)
    unlink(paste0(tmpBase,".bim"),recursive=TRUE,force=TRUE)
    unlink(paste0(tmpBase,".fam"),recursive=TRUE,force=TRUE)
    
    return(gdsfile)
}

.writeGdsFile <- function(input,snpObj,gdsfile,selection) {
    disp("Reading PLINK files with SNPRelate framework and storing output tp ",
        gdsfile)
    # We also need to read a GDS file for LD and IBD filtering
    # If selection is not NULL, then the GDS file will not be 
    # correct... We have to write temporary BED/BIM/FAMs...
    if (!is.null(selection)) {
        tmpBase <- tempfile()
        disp("  SNP/Sample selection was made! Writing temporary PLINK files ",
            "for correct GDS creation...")
        write.plink(
            file.base=tmpBase,
            snps=snpObj$genotypes,
            pedigree=snpObj$fam$pedigree,
            id=rownames(snpObj$fam),
            father=snpObj$fam$father,
            mother=snpObj$fam$mother,
            sex=snpObj$fam$sex,
            phenotype=snpObj$fam$affected,
            chromosome=snpObj$map$chromosome,
            position=snpObj$map$position,
            allele.1=snpObj$map$allele.1,
            allele.2=snpObj$map$allele.2
        )
        snpgdsBED2GDS(paste0(tmpBase,".bed"),paste0(tmpBase,".fam"),
            paste0(tmpBase,".bim"),gdsfile,family=TRUE)
        unlink(paste0(tmpBase,".bed"))
        unlink(paste0(tmpBase,".bim"))
        unlink(paste0(tmpBase,".fam"))
    }
    else
        snpgdsBED2GDS(input$bed,input$fam,input$bim,gdsfile,
            family=TRUE)
}

filterGWAS <- function(obj,filters=getDefaults("filters"),imputeMissing=TRUE,
    imputeMode=c("single","split"),rc=NULL) {
    filters <- .checkFilters(filters)
    imputeMode <- imputeMode[1]
    
    # Later input type may be more native, e.g. read from CSV files
    if (is(assay(obj,1),"SnpMatrix"))
        return(.filterWithSnpStats(obj,filters,imputeMissing,imputeMode,rc))
    else if (is(assay(obj,1),"bigsnp"))
        return(.filterWithBigSnpr(obj,filters,imputeMissing))
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
    reverse=FALSE,perChr=FALSE,overwrite=TRUE) {
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
    if (reverse) {
        disp("Switching alleles...")
        gen <- switch.alleles(gen,seq_len(ncol(gen)))
    }
    
    if (salvage) {
        disp("Trying to salvage SNPs with missing locations")
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
    disp("Preflight...")
    map$chromosome <- as.character(map$chromosome)
    if (any(map$chromosome %in% c("23","24","25","26"))) {
        map$chromosome[map$chromosome=="23"] <- "X"
        map$chromosome[map$chromosome=="24"] <- "Y"
        map$chromosome[map$chromosome=="25"] <- "XY"
        map$chromosome[map$chromosome=="26"] <- "MT"
    }
    
    # Has info score from imputation? If yes must be written to separate file
    hasInfo <- ifelse("info" %in% names(map),TRUE,FALSE)
    
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
    if (perChr) {
        gen <- gen[,!na]
        map <- map[!na,,drop=FALSE]
        S <- split(seq_len(nrow(map)),map$chromosome)
        if (!any(grepl("chr",names(S))))
            names(S) <- paste("chr",names(S),sep="")
        lapply(names(S),function(n,R) {
            disp("===== Writing files for chromosome ",
                gsub("chr","",n,ignore.case=TRUE))
            if (file.exists(paste0(outBase,"_",n,".bed")) && !overwrite)
                disp("  fileset ",paste0(outBase,"_",n)," exists! Skipping...")
            else {
                ii <- R[[n]]
                write.plink(
                    file.base=paste0(outBase,"_",n),
                    snps=gen[,ii],
                    pedigree=fam$pedigree,
                    id=rownames(fam),
                    father=fam$father,
                    mother=fam$mother,
                    sex=fam$sex,
                    phenotype=phenoVec,
                    chromosome=map$chromosome[ii],
                    position=map$position[ii],
                    allele.1=map$allele.1[ii],
                    allele.2=map$allele.2[ii]
                )
            }
            
            if (hasInfo) {
                if (file.exists(paste0(outBase,"_",n,".impinfo")) && !overwrite)
                    disp("  imputation info file ",paste0(outBase,"_",n,
                        ".impinfo")," exists! Skipping...")
                else {
                    impinfo <- as.data.frame(map[ii,c("snp.name","info")])
                    names(impinfo)[1] <- "snp_name"
                    write.table(impinfo,paste0(outBase,"_",n,".impinfo"),
                        sep="\t",quote=FALSE,row.names=FALSE)
                }
            }
        },S)
    }
    else {
        if (file.exists(paste0(outBase,".bed")) && !overwrite)
            disp("  fileset ",outBase," exists! Skipping...")
        else {
            write.plink(
                file.base=outBase,
                snps=gen[,!na],
                #pedigree=fam$pedigree,
                pedigree=rownames(fam),
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
        
        if (hasInfo) {
            if (file.exists(paste0(outBase,".impinfo")) && !overwrite)
                disp("  imputation info file ",paste0(outBase,".impinfo"),
                    " exists! Skipping...")
            else {
                impinfo <- as.data.frame(map[!na,c("snp.name","info")])
                names(impinfo)[1] <- "snp_name"
                write.table(impinfo,paste0(outBase,".impinfo"),sep="\t",
                    quote=FALSE,row.names=FALSE)
            }
        }
    }
    
    # For some reason the process takes a lot of memory...
    gc(verbose=FALSE)
    invisible(return(NULL))
}
