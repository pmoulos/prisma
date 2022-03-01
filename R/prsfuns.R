# More arguments to come, like PRS formula like in PRSice etc.
PRS <- function(gwe,snpSelection,type=c("avg","sum","std")) {
    # Some validation here
    if (!is(gwe,"GWASExperiment"))
        stop("Input object must be a GWASExperiment object!")
    if (!any(c("effect","effect_weight","OR") %in% colnames(snpSelection)))
        stop("snpSelection does not seem to be an output from prisma or ",
            "aggregatePrsMarkers functions! Please check!")
    
    type <- type[1]
    .checkTextArgs("PRS calculation type (type)",type,c("avg","sum","std"),
        multiarg=FALSE)
    
    ii <- grep("effect",colnames(snpSelection))
    X <- t(as(genotypes(gwe[rownames(snpSelection),,drop=FALSE]),"numeric"))
    return(.prs(X,snpSelection[,ii],type))
}

.prs <- function(g,e,s=c("avg","sum","std")) {
    s <- s[1]
    switch(s,
        avg = {
            return((g %*% e)/length(e))
        },
        sum = {
            return(g %*% e)
        },
        std = {
            p <- g %*% e
            return((p - mean(p))/sd(p))
        }
    )
}

runPRS <- function(base,target=base,response,covariates=NULL,pcs=FALSE,
    methods=c("lassosum","prsice"),lassosumOpts=getDefaults("lassosum"),
    prsiceOpts=getDefaults("prsice"),wspace=NULL,rc=NULL) {
    # Validate arguments
    .checkTextArgs("PRS algorithm(s) (methods)",methods,c("lassosum","prsice"),
        multiarg=TRUE)
    lassosumOpts <- .checkPrsArgs(lassosumOpts,"lassosum")
    prsiceOpts <- .checkPrsArgs(prsiceOpts,"prsice")
    
    prsResults <- vector("list",length(methods))
    names(prsResults) <- methods
    disp("\nStarting PRS analysis with ",paste0(methods,
        collapse=", "))
    for (m in methods) {
        disp("\n========== Running PRS analysis with ",m)
        switch(m,
            lassosum = {
                prsResults[[m]] <- lassosumPRS(base,target,response,covariates,
                    pcs,wspace,lassosumOpts$anc,lassosumOpts$valid)
            },
            prsice = {
                runtype <- "calculate"
                prsResults[[m]] <- prsicePRS(base,target,response,covariates,
                    pcs,runtype,wspace,prsiceOpts$clump_kb,prsiceOpts$clump_r2,
                    prsiceOpts$clump_p,prsiceOpts$score,prsiceOpts$perm,
                    prsiceOpts$seed,rc)
            }
        )
    }
    disp("\nRunning PRS algorithms finished! Processing the results...")
    
    theBetas <- matrix(0,nrow(target),length(methods))
    rownames(theBetas) <- rownames(target)
    colnames(theBetas) <- methods
    for (m in methods) {
        # There are extreme cases where there are 1-2 mismatches
        keep <- intersect(prsResults[[m]]$variant_id,rownames(target))
        theBetas[keep,m] <- prsResults[[m]][keep,"effect_weight"]
    }
    
    npcs <- ifelse(pcs && .hasPcaCovariates(base),ncol(pcaCovariates(base)),0)
    prsbetas(target,response,covariates,npcs) <- theBetas
    
    disp("Done!")
    return(target)
}

lassosumPRS <- function(base,target=base,response,covariates=NULL,pcs=FALSE,
    wspace=NULL,anc=c("eur","asn","afr"),valid=c("auto","stack","split")) {
    # Can run PRS? Has necessary fields?
    .canRunPrs(base,response)
    .canRunPrs(target,response,isBase=FALSE)
    
     # Workspace valid?
    wspace <- .validateWorkspacePath(wspace,"lassosum")
    
    # Ancestry for LD blocks based on 1000G provided by lassosum
    anc <- anc[1]
    valid <- valid[1]
    .checkTextArgs("Ancestry (anc)",anc,c("eur","asn","afr"),multiarg=FALSE)
    .checkTextArgs("Validation type (valid)",valid,c("auto","stack","split"),
        multiarg=FALSE)
    
    if (valid == "auto")
        valid <- ifelse(ncol(base) >= 1000,"split","stack")
    
    # Prepare the phenotypes for later prediction/validation
    p <- phenotypes(target)
    chResCov <- .validateResponseAndCovariates(p,response,covariates)
    response <- chResCov$res
    covariates <- chResCov$cvs
    pheno <- p[,c("FID","IID",response)]
    covar <- p[,c("FID","IID",covariates)]
    # Make sure FID and IID matches, otherwise cannot be properly handled by
    # snpStats::write.plint
    if (!identical(pheno[,"FID"],pheno[,"IID"])) {
        pheno[,"IID"] <- pheno[,"FID"]
        covar[,"IID"] <- covar[,"FID"]
    }
    if (pcs) { # Include robust PCs in the model
        if (.hasPcaCovariates(target)) {
            pcov <- pcaCovariates(target)
            covar <- cbind(covar,pcov)
            npcs <- ncol(pcov)
        }
        else {
            disp("PC covariates requested in the model, but no calculated PC ",
                "covariates found in the target!\nWill try to guess from the ",
                "provided base...")
            if (.hasPcaCovariates(base)) {
                f <- .guessPcaParamsFromObject(base)
                target <- calcPcaCovar(target,f$LD,rc=rc)
                pcov <- pcaCovariates(target)
                covar <- cbind(covar,pcov)
                npcs <- ncol(pcov)
            }
            else {
                 #warning("PC covariates requested in the model, but not ",
                 #   "calculated in target or base objects! Ignoring...",
                 #   immediate.=TRUE)
                 #npcs <- 0
                 stop("PC covariates requested in the model, but not ",
                    "calculated in target or base objects!\nInput objects are ",
                    "probably broken!")
             }
        }
    }
    
    # If all ok, prepare the run
    prepList <- .prepareLassosumRun(base,target,response,covariates,npcs,wspace,
        toupper(anc))
    
    # Transform the P-values into correlation
    ss <- prepList$ss
    cor <- p2cor(p=ss$pvalue,n=ncol(base),sign=ss$effect)
    
    # We need to attach marker ids to ss for later extraction
    marks <- rownames(ss)
    names(marks) <- paste0(ss$chromosome,"_",ss$position,"_",ss$allele.1,"_",
        ss$allele.2)
    
    # Run the lassosum pipeline
    disp("Running lassosum PRS algorithm")
    cwd <- getwd()
    setwd(wspace)
    out <- lassosum.pipeline(
        cor=cor,
        chr=ss$chromosome,
        pos=ss$position,
        A1=ss$allele.1,
        A2=ss$allele.2,
        #ref.bfile=prepList$base,
        ref.bfile=prepList$target,
        test.bfile=prepList$target,
        LDblocks=paste0(toupper(anc),".",prepList$gb)
    )
    # Adjust the output - actual location of temporary PLINK file
    out$test.bfile <- file.path(wspace,prepList$target)
    # Marker names
    tmp <- out$sumstats
    rownames(tmp) <- paste0(tmp$chr,"_",tmp$pos,"_",tmp$A1,"_",tmp$A2)
    # There are some extreme cases of not agreeing marker names...
    keep <- intersect(rownames(tmp),names(marks))
    tmp <- tmp[keep,,drop=FALSE]
    rownames(tmp) <- marks[keep]
    out$sumstats <- tmp
    setwd(cwd)

    # PRS? If sample size is >1000 then splitvalidate else validate
    disp("Running lassosum validation and PRS calculation")
    if (valid == "stack")
        val <- lassosum::validate(out,pheno=pheno,covar=covar,plot=FALSE)
    else if (valid == "split")
        val <- lassosum::splitvalidate(out,pheno=pheno,covar=covar,plot=FALSE)
    
    disp("Done! R^2 is: ",max(val$validation.table$value)^2)
    return(.extractLassosumPrsComponents(out,val,genome(base)))
}

# --no-clump, --extract to be used in ensembl snps
# Probably more arguments related to output will be added
# A mode arg will be added, controlling new PRS or with provided SNPs that will
# be extracted from the obj along with the parameters above
prsicePRS <- function(base,target=base,response,covariates=NULL,pcs=FALSE,
    mode=c("calculate","apply"),wspace=NULL,clump_kb=250,clump_r2=0.1,clump_p=1,
    score=c("avg","sum","std","con-std"),perm=10000,seed=42,rc=NULL) {
    # Can run PRS? Has necessary fields?
    .canRunPrs(base,response)
    .canRunPrs(target,response,isBase=FALSE)
    
    # Check further parameters
    score <- tolower(score[1])
    mode <- tolower(mode[1])
    .checkNumArgs("Clumping distance (clump_kb)",clump_kb,"numeric",0,"gt")
    .checkNumArgs("Clumping R2 (clump_r2)",clump_r2,"numeric",c(0,1),"botheq")
    .checkNumArgs("Clumping p-value (clump_p)",clump_p,"numeric",c(0,1),"both")
    .checkNumArgs("Number of permutations (perm)",perm,"numeric",0,"gte")
    .checkNumArgs("Number of permutations (perm)",perm,"numeric",0,"gte")
    .checkTextArgs("PRS scoring scheme (score)",score,c("avg","sum",
        "std","con-std"),multiarg=FALSE)
    .checkTextArgs("PRSice run mode (mode)",mode,c("calculate","apply"),
        multiarg=FALSE)
    if (!is.numeric(seed))
        stop("The PRSice permutation seed must be a number!")
    
    # Workspace valid?
    wspace <- .validateWorkspacePath(wspace,"prsice")
    
    # Prepare the phenotypes for later prediction/validation
    p <- phenotypes(target)
    chResCov <- .validateResponseAndCovariates(p,response,covariates)
    response <- chResCov$res
    covariates <- chResCov$cvs
    pheno <- p[,c("FID","IID",response)]
    covar <- p[,c("FID","IID",covariates)]
    # Make sure FID and IID matches, otherwise cannot be properly handled by
    # snpStats::write.plink
    if (!identical(pheno[,"FID"],pheno[,"IID"])) {
        pheno[,"IID"] <- pheno[,"FID"]
        covar[,"IID"] <- covar[,"FID"]
    }
    if (pcs) { # Include robust PCs in the model
        if (.hasPcaCovariates(target)) {
            pcov <- pcaCovariates(target)
            covar <- cbind(covar,pcov)
            npcs <- ncol(pcov)
        }
        else {
            disp("PC covariates requested in the model, but not calculated ",
                "PC covariates found in the target!\nWill try to guess from ",
                "the provided base...")
            if (.hasPcaCovariates(base)) {
                f <- .guessPcaParamsFromObject(base)
                #target <- .wrapPcaWithSnpStatsLd(target,f,rc)
                target <- calcPcaCovar(target,f$LD,rc=rc)
                pcov <- pcaCovariates(target)
                covar <- cbind(covar,pcov)
                npcs <- ncol(pcov)
            }
            else {
                #warning("PC covariates requested in the model, but not ",
                #    "calculated in target or base objects! Ignoring...",
                #    immediate.=TRUE)
                #npcs <- 0
                stop("PC covariates requested in the model, but not ",
                    "calculated in target or base objects!\nInput objects are ",
                    "probably broken!")
            }
        }
    }
    
    # If all ok, prepare the run
    prepList <- .preparePrsiceRun(base,target,pheno,covar,response,covariates,
        npcs,wspace)
    
    # Run the PRSice pipeline
    disp("Running PRSice PRS algorithm in ",mode," mode")
    prsice <- .getToolPath("prsice")
    prsiceScript <- file.path(dirname(prsice),"PRSice.R")
    binaryTarget <- ifelse(.maybeBinaryForBinomial(pheno[,response]),"T","F")
    threads <- .coresFrac(rc)
    
    if (mode == "calculate")
        command <- paste(
            paste0("Rscript ",prsiceScript," \\"),
            paste0("  --prsice ",prsice," \\"),
            paste0("  --base ",prepList$base," \\"),
            "  --beta \\",
            paste0("  --binary-target ",binaryTarget," \\"),
            paste0("  --pheno ",prepList$pheno," \\"),
            paste0("  --target ",prepList$target," \\"),
            paste0("  --cov ",prepList$covar," \\"),
            paste0("  --clump-kb ",clump_kb," \\"),
            paste0("  --clump-r2 ",clump_r2," \\"),
            paste0("  --clump-p ",clump_p," \\"),
            paste0("  --score ",score," \\"),
            paste0("  --seed ",seed," \\"),
            paste0("  --out ",prepList$out," \\"),
            "  --print-snp",
            sep="\n"
        )
    else if (mode == "apply") {
        command <- paste(
            paste0("Rscript ",prsiceScript," \\"),
            paste0("  --prsice ",prsice," \\"),
            paste0("  --base ",prepList$base," \\"),
            "  --beta \\",
            paste0("  --binary-target ",binaryTarget," \\"),
            paste0("  --pheno ",prepList$pheno," \\"),
            paste0("  --target ",prepList$target," \\"),
            paste0("  --cov ",prepList$covar," \\"),
            paste0("  --no-clump \\"),
            paste0("  --score ",score," \\"),
            paste0("  --seed ",seed," \\"),
            paste0("  --out ",prepList$out," \\"),
            "  --print-snp \\",
            "  --bar-levels 1 \\",
            "  --fastscore",
            sep="\n"
        )
    }
    
    if (threads > 1) {
        addThreads <- paste0("  --thread ",threads)
        command <- paste0(command," \\\n",addThreads)
    }
    if (perm > 0 && mode == "calculate") {
        addPerm <- paste0("  --perm ",perm)
        command <- paste0(command," \\\n",addPerm)
    }
    ## Not working - future work: create LD from 1000 genomes
    #if (ncol(base) < 500) {
    #   addLd <- " --ld ",prepList$ldb
    #   command <- paste0(command,"\\\n",addLd)
    #}
    
    message("\nExecuting:\n",command)
    cwd <- getwd()
    setwd(wspace)
    out <- tryCatch({
        system(command,ignore.stdout=TRUE,ignore.stderr=TRUE)
        if (prismaVerbosity() == "full") {
            logfile <- file.path(wspace,paste0(prepList$out,".log"))
            log <- .formatSystemOutputForDisp(readLines(logfile))
            disp("\nPRSice output is:\n")
            disp(paste(log,collapse="\n"))
        }
    },error=function(e) {
        message("Caught error: ",e$message)
        return(1L)
    },interrupt=function() {
        setwd(cwd)
    },finally="")
    setwd(cwd)
    
    if (!.isEmpty(out) && out == 1L) {
        message("PRSice run failed! Will return empty result...")
        return(.emptyVariantsDf("pgs"))
    }
    
    # Read the summary file so as to extract p-value threshold and R2
    sumFile <- file.path(wspace,paste0("prsice_out_",prepList$runid,".summary"))
    sumData <- read.delim(sumFile)
    
    disp("Done!\n")
    disp("Adjusted PRSice R^2 is: ",sumData$PRS.R2)
    disp("Full model R^2 is: ",sumData$Full.R2)
    disp("Null model R^2 is: ",sumData$Null.R2)
    
    return(.extractPrsicePrsComponents(base,prepList$runid,sumData$Threshold,
        response,covariates,npcs,wspace))
}

.extractLassosumPrsComponents <- function(out,val,gb) {
    # Construct initial data frame (we have matched marker names with locations)
    df <- out$sumstats[,c("chr","pos","A1","A2")]
    colnames(df) <- c("chromosome","position","risk_allele","reference_allele")
    df$variant_id <- rownames(df)
    
    # Assuming that betas are aligned with sumstats (they are)
    nonZero <- val$best.beta!=0
    df <- df[nonZero,,drop=FALSE]
    # Until we get an answer from lassosum authors, we reverse betas so as to be
    # concordant with the classical X %*% betas PRS definition
    ### It was WRONG as snpStats was reversing genotypes
    #df$effect_weight <- -val$best.beta[nonZero]
    ###
    df$effect_weight <- val$best.beta[nonZero]
    # and conversion to odds ratios
    df$OR <- exp(df$effect_weight)
    df$locus_name <- rep(NA,nrow(df))
    
    # Final alignment with the external API fetch outcomes from PGS catalog
    df <- df[,c("chromosome","position","variant_id","risk_allele",
        "reference_allele","locus_name","effect_weight","OR")]
    if (is.null(gb))
        gb <- "nr"
    df$asm <- rep(gb,nrow(df))
    
    return(df)
}

.extractPrsicePrsComponents <- function(obj,runid,pcut,res,cvs,npcs,wspace) {
    # Read in detailed SNP output
    snpFile <- file.path(wspace,paste0("prsice_out_",runid,".snp"))
    snpData <- read.delim(snpFile)
    snps <- as.character(snpData$SNP[snpData$P<=pcut])
    
    # Construct the output
    obj <- obj[snps,]
    df <- as.data.frame(gfeatures(obj))
    # Basic info
    out <- df[,c("chromosome","position","snp.name","allele.1","allele.2")]
    # Attach effects
    e <- effects(obj,res,cvs,npcs)
    pri <- .getGwaLinArgPrior()
    if (any(pri %in% colnames(e)))
        out$effect_weight <- e[,pri[which(pri %in% colnames(e))[1]]]
    out$OR <- exp(out$effect_weight)
    out$locus_name <- rep(NA,nrow(out))
    
    # Final alignment with the external API fetch outcomes from PGS catalog
    names(out)[c(3,4,5)] <- c("variant_id","risk_allele","reference_allele")
    out <- out[,c("chromosome","position","variant_id","risk_allele",
        "reference_allele","locus_name","effect_weight","OR")]
    gb <- genome(obj)
    if (is.null(gb))
        gb <- "nr"
    out$asm <- rep(gb,nrow(out))
    
    return(out)
}

.prepareLassosumRun <- function(base,target,pheno,cvs,npcs,wspace,anc) {
    disp("Preparing lassosum run at workspace directory: ",wspace)
    
    # Prepare the sum stat
    disp("  preparing summary statistics...")
    statsCoords <- gfeatures(base)
    statsCoords <- as.data.frame(statsCoords[,c("chromosome","position",
        "allele.1","allele.2")])
    preP <- pvalues(base,pheno,cvs,npcs)
    p <- preP[,ncol(preP)]
    e <- effects(base,pheno,cvs,npcs)
    
    # Until we find a better way of summarizing effects across methods,
    # prioritize effect selection according to algorithm popularity
    pri <- .getGwaLinArgPrior()
    if (any(pri %in% colnames(e)))
        oe <- e[,pri[which(pri %in% colnames(e))[1]]]
    sumStat <- data.frame(statsCoords,effect=oe,pvalue=p)

    # Prepare PLINK files
    runId <- .randomString()
    #disp("  preparing PLINK temporary files for reference...")
    #bfileBase <- paste0("plink_lassosum_base_",runId)
    #writePlink(base,pheno,outBase=file.path(wspace,bfileBase))
    disp("  preparing PLINK temporary files for target...")
    #if (identical(colnames(base),colnames(target)))
    #    bfileTarget <- bfileBase
    #else {
    #    bfileTarget <- paste0("plink_lassosum_target_",runId)
    #    writePlink(target,pheno,outBase=file.path(wspace,bfileTarget))
    #}
    bfileTarget <- paste0("plink_lassosum_target_",runId)
    writePlink(target,pheno,outBase=file.path(wspace,bfileTarget))
    
    # Copy LD block files to workspace directory
    disp("  copying linkage disequilibrium block files...")
    ldHome <- system.file("data",package="lassosum")
    gb <- ifelse(is.null(genome(base)),"hg19",genome(base))
    fromLd <- file.path(ldHome,paste0("Berisa.",anc,".",gb,".bed"))
    toLd <- file.path(wspace,paste0("Berisa.",anc,".",gb,".bed"))
    file.copy(from=fromLd,to=toLd,overwrite=TRUE)
    
    disp("Preparation done!")
    #return(list(ss=sumStat,base=bfileBase,target=bfileTarget,ldb=toLd,gb=gb))\
    return(list(ss=sumStat,base=bfileTarget,target=bfileTarget,ldb=toLd,gb=gb))
}

.preparePrsiceRun <- function(base,target,pheno,covars,res,cvs,npcs,wspace) {
    disp("Preparing PRSice run at workspace directory: ",wspace)
    
    # Prepare the sum stat
    disp("  preparing summary statistics for the base file...")
    statsInfo <- gfeatures(base)
    statsInfo <- as.data.frame(statsInfo[,c("chromosome","position",
        "snp.name","allele.1","allele.2")])
    #if (identical(colnames(base),colnames(target)))
    #    p <- rep(1,nrow(target))
    #else {
        preP <- pvalues(base,res,cvs,npcs)
        p <- preP[,ncol(preP)]
    #}
    e <- effects(base,res,cvs,npcs)
    
    # Until we find a better way of summarizing effects across methods,
    # prioritize effect selection according to algorithm popularity
    pri <- .getGwaLinArgPrior()
    if (any(pri %in% colnames(e)))
        oe <- e[,pri[which(pri %in% colnames(e))[1]]]
    sumStat <- data.frame(statsInfo,BETA=oe,P=p)
    # PRSice needs(?) MAF and INFO, we fake them as input data are supposed to 
    # be already filtered
    nr <- nrow(sumStat)
    sumStat$MAF <- rep(0.5,nr)
    sumStat$INFO <- rep(1,nr)
    # Also sample size, since it's included in their examples
    sumStat$N <- rep(ncol(base),nr)
    # Let's give the default expected names by PRSice
    colnames(sumStat)[seq_len(5)] <- c("CHR","BP","SNP","A1","A2")
    
    # Random run id
    runId <- .randomString()
    
    # A file must be written - whitespace separator, not tab!
    baseFile <- file.path(wspace,paste0("base_",runId))
    write.table(sumStat,file=baseFile,quote=FALSE,row.names=FALSE)
    
    # Prepare the covariates - a file must be written with space delimited
    disp("  preparing the covariates file...")
    covarFile <- file.path(wspace,paste0("covar_",runId))
    # ALL the covariates must be numeric - there may be factors which are 
    # character - R is deceiving about this as it handles
    covars <- .checkAndCorrectFactorsFormat(covars)
    write.table(covars,file=covarFile,quote=FALSE,row.names=FALSE)
    
    # Prepare the phenotype - a file must be written with space delimited
    disp("  preparing the phenotype file...")
    phenoFile <- file.path(wspace,paste0("pheno_",runId))
    pheno <- .checkAndCorrectFactorsFormat(pheno)
    write.table(pheno,file=phenoFile,quote=FALSE,row.names=FALSE)

    # Prepare PLINK files
    disp("  preparing PLINK temporary files - target...")
    resp <- names(pheno)[3]
    targetBase <- file.path(wspace,paste0("plink_prsice_",runId))
    #writePlink(target,resp,outBase=targetBase,reverse=TRUE)
    writePlink(target,resp,outBase=targetBase)
    
    ## Let's see what will happen with lassosum LD blocks...
    #anc <- "EUR" # Hardcode for now
    #disp("  copying linkage disequilibrium block files...")
    #ldHome <- system.file("data",package="lassosum")
    gb <- ifelse(is.null(genome(base)),"hg19",genome(base))
    #fromLd <- file.path(ldHome,paste0("Berisa.",anc,".",gb,".bed"))
    #toLd <- file.path(wspace,paste0("Berisa.",anc,".",gb,".bed"))
    #file.copy(from=fromLd,to=toLd,overwrite=TRUE)
    
    disp("Preparation done!")
    return(list(base=baseFile,target=targetBase,pheno=phenoFile,covar=covarFile,
        ldb=NULL,out=paste0("prsice_out_",runId),runid=runId,gb=gb))
}

.getGwaLinArgPrior <- function() {
    return(c("snptest","statgen","plink","glm","rrblup"))
}

.canRunPrs <- function(o,p,isBase=FALSE) {
    # Object class?
    if (!is(o,"GWASExperiment"))
        stop("PRS input must be an object of class GWASExperiment!")
    # Has an association test been performed?
    if (isBase && (is.null(pvalues(o)) || is.null(effects(o))))
        stop("At least one association test must be performed prior to ",
            "building PRS!")
    # The given phenotype exists?
    ph <- phenotypes(o)
    if (is.null(ph))
        stop("PRS cannot be performed without phenotypes!")
    else {
        if (is.numeric(p)) {
            if (p > ncol(ph))
                stop("The phenotype index ",p," is larger than the object's ",
                    "available phenotypes!")
        }
        else if (is.character(p)) {
            if (!(p %in% colnames(ph)))
                stop("The phenotype ",p," cannot be found in the object's ",
                    "available phenotypes (",paste(colnames(ph),collapse=", "),
                    ")!")
        }
    }
}

.checkPrsArgs <- function(args,prs=c("lassosum","prsice")) {
    prs <- prs[1]
    
    # Allowed and given values
    defaults <- getDefaults(prs)
    allowed <- names(defaults)
    given <- names(args)
    
    # Check if illegal filter names have been provided
    check <- given %in% allowed
    if (!any(check))
        stop("No valid ",prs," parameter name found!")
    if (!all(check)) {
        warning("The following ",prs," parameters names are invalid and will ",
            "be ignored:\n",paste(given[!check],collapse=", "))
        args <- args[given[check]]
    }
    
    # Proceed with per algorithm check, will stop here if something wrong
    switch(prs,
        lassosum = {
            .checkLassosumArgs(args)
        },
        prsice = {
            .checkPrsiceArgs(args)
        }
    )
    
    # If all OK
    return(.setArg(defaults,args))
}

.checkLassosumArgs <- function(a) {
    if (!.isEmpty(a$anc))
        .checkTextArgs("Ancestry (anc)",a$anc,c("eur","asn","afr"),
            multiarg=FALSE)
    if (!.isEmpty(a$valid))
        .checkTextArgs("Validation type (valid)",a$valid,c("auto","stack",
            "split"),multiarg=FALSE)
}

.checkPrsiceArgs <- function(a) {
    if (!.isEmpty(a$clump_kb))
        .checkNumArgs("Clumping distance (clump_kb)",a$clump_kb,"numeric",0,
            "gt")
    if (!.isEmpty(a$clump_r2))
        .checkNumArgs("Clumping R2 (clump_r2)",a$clump_r2,"numeric",c(0,1),
            "botheq")
    if (!.isEmpty(a$clump_p))
        .checkNumArgs("Clumping p-value (clump_p)",a$clump_p,"numeric",c(0,1),
            "both")
    if (!.isEmpty(a$perm))
        .checkNumArgs("Number of permutations (perm)",a$perm,"numeric",0,"gte")
    if (!.isEmpty(a$score))
        .checkTextArgs("PRS scoring scheme (score)",a$score,c("avg","sum","std",
            "con-std"),multiarg=FALSE)
    if (!.isEmpty(a$seed) && !is.numeric(a$seed))
        stop("The PRSice permutation seed must be a number!")
}
