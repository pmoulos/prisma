# Convert GWASExperiment to gData (statgenGWAS)
GWASExperiment2gData <- function(obj,covariates=NULL,pcs=FALSE,reverse=FALSE) {
    if (!requireNamespace("statgenGWAS"))
        stop("R package statgenGWAS is required!")
    
    # Phenotypes - drop those which are not numeric for now... and check with
    # covariates... If not everything satisfied, we do not proceed.
    pheno <- phenotypes(obj)
    chn <- vapply(names(pheno),function(x,p) {
        return(is.character(p[[x]]))
    },logical(1),pheno)
    if (any(chn)) {
        warning("Some phenotypes are not numeric (",
            paste0(names(pheno)[which(chn)],collapse=", "),")!\nThese will ",
            "not be passed to the resulting gData object.",immediate.=TRUE)
        pheno <- pheno[,!chn,drop=FALSE]
    }
    
    if (pcs) { # Include precalculated PCs if any
        if (.hasPcaCovariates(obj)) {
            pcov <- pcaCovariates(obj)
            pheno <- cbind(pheno,pcov)
            covariates <- c(covariates,colnames(pcov))
        }
        else
            warning("PC covariates reqeusted in the model, but no calculated ",
                "PC covariates found! Ignoring...",immediate.=TRUE)
    }
    
    covar <- NULL
    if (!is.null(covariates)) {
        if (!(all(covariates %in% names(pheno)))) {
            no <- which(!(covariates %in% names(pheno)))
            stop("Covariate(s) ",paste0(covariates[no],collapse=", "),
                " not found in gData phenotypes!\nHave they been dropped ",
                "because they were not numeric?")
        }
        
        covar <- pheno[,covariates,drop=FALSE]
        pheno <- pheno[,-match(covariates,names(pheno)),drop=FALSE]
    }
    pheno <- data.frame(genotype=rownames(pheno),pheno)
    
    
    # Genotypes
    geno <- t(genotypes(obj))
    if (reverse)
        geno <- switch.alleles(geno,seq_len(ncol(geno)))
    geno <- as(geno,"numeric")
    
    # Map
    map <- as.data.frame(gfeatures(obj))
    map <- map[,c("chromosome","position")]
    names(map) <- c("chr","pos")
    
    return(createGData(geno=geno,map=map,pheno=pheno,covar=covar))
}

GWASExperiment2GDS <- function(obj) {
    if (is.null(gdsfile(obj)))
        stop("A valid GDS file location must be provided! Use gdsfile()")
    
    # Write temporary PLINK
    tmp <- tempfile()
    writePlink(obj,outBase=tmp)
    
    # Write GDS file
    bed <- paste0(tmp,".bed")
    fam <- paste0(tmp,".fam")
    bim <- paste0(tmp,".bim")
    snpgdsBED2GDS(bed,fam,bim,gdsfile(obj),family=TRUE)
}

GWASExperimentLiftOver <- function(obj,from,to) {
    if (!requireNamespace("liftOver"))
        stop("Bioconductor package liftOver is required!")
    if (!is(obj,"GWASExperiment"))
        stop("Input object must be a GWASExperiment object!")
    
    from <- tolower(from[1])
    to <- tolower(to[1])
    .checkTextArgs("Source genome (from)",from,c("hg19","hg38"),
        multiarg=FALSE)
    .checkTextArgs("Destination genome (to)",to,c("hg19","hg38"),
        multiarg=FALSE)
    
    if (from == to) # Nothing to do
        return(obj)
    
    # 1. Convert gfeatures to GPos
    f <- gfeatures(obj)
    gf <- data.frame(
        chromosome=f$chromosome,
        start=f$position,
        end=f$position
    )
    gfp <- GPos(gf,stitched=FALSE)
    names(gfp) <- rownames(f)
    
    # 2. Lift over
    disp("Lifting over features")
    lifted <- .liftOverSNPs(gfp,from,to)
    
    # 3. Common reference
    common <- intersect(rownames(f),names(lifted))
    if (length(common) < nrow(f)) {
        disp(nrow(f) - length(common)," SNPs did not map to ",to)
        disp("These are:",level="full")
        disp("    ",paste(setdiff(rownames(f),common),collapse="\n    "),
            level="full")
    }
    
    # 4. Drop annotation and genotypes for SNPs that do not match and reorder
    obj <- obj[common,]
    f <- f[common,]
    lifted <- lifted[common]
    f$chromosome <- as.character(seqnames(lifted))
    f$position <- start(lifted)
    gfeatures(obj) <- f
    genome(obj) <- to
    
    return(obj)
}

guessHumanGenomeVersion <- function(o) {
    # Sample 5 random SNPs... If not with rs, then no possibility
    s <- sample(rownames(o),5)
    s <- s[grepl("^rs",s)]
    if (length(s) == 0) {
        warning("Cannot determine genome version without rs ids!",
            immediate.=TRUE)
        return(NA)
    }
        
    # rsnps fetch from hg38
    hits <- suppressWarnings(rsnps::ncbi_snp_query(s))
    pos <- gfeatures(o)[s,"position"]
    hpo <- hits$bp
    if (!all(pos %in% hpo))
        return("hg19")
    else
        return("hg38")
}
