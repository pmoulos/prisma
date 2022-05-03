# Make intervals for IMPUTE2
# The logic is
# 1. Find min max SNP positions from gfeatures per chromosome
# 2. Create intervals (see recoup function)
# 3. Export gwe per chromosome
# 4. Convert to PED... PLINK needed...
# 5. Convert to IMPUTE gen... GTOOL needed...
# 6. Run per chromosome
# We are going to need some kind of automation for reference panels

# The output should ideally be an imputed GWASExperiment with some additional
# INFO score in the gfeatures
# Ideally, some kind of break and continue mechanism should be in place like
# in prs pipelines
extendGWAS <- function(obj,intSize=1e+6,wspace=NULL,refSpace=NULL,
    continue=FALSE,cleanup=c("none","intermediate","all"),runId=NULL,
    convTool=c("gtool","qctool"),rc=NULL) {
    if (!.toolAvailable("impute"))
        stop("IMPUTE2 program not found in the system!")
    
    cleanup <- cleanup[1]
    convTool <- convTool[1]
    .checkTextArgs("Cleanup option (cleanup)",cleanup,
        c("none","intermediate","all"),multiarg=FALSE)
    .checkTextArgs("Format conversion tool (convTool)",convTool,
        c("gtool","qctool"),multiarg=FALSE)
    
    if (convTool == "gtool" && !.toolAvailable("gtool"))
        stop("GTOOL program not found in the system!")
    if (convTool == "qctool" && !.toolAvailable("qctool"))
        stop("QCTOOL program not found in the system!")
    
    # Validate workspace
    wspace <- .validateWorkspacePath(wspace,"impute")
    
    # Check refSpace contains 1000GP files
    # ...
    
    # Check chromosomes are there
    map <- gfeatures(obj)
    if (!("chromosome" %in% names(map))) # Trying to impute with what?
        stop("Cannot extend the dataset without location info!")
    
    # We need a runid for potentially working in the same space
    if (is.null(runId))
        runId <- .randomString()
    
    # Better to create a runId subdir to do the processing instead of attaching
    # to filename and move the final beds outside this in the end
    wspaceBak <- wspace
    wspace <- file.path(wspace,runId)
    if (!dir.exists(wspace))
        dir.create(wspace,recursive=TRUE,mode="0755",showWarnings=FALSE)
    
    # Prepare the main input files to impute per chromosome
    if (convTool == "gtool")
        gensam <- .prepareInputFilesForImpute2(obj,wspace,rc)
    else if (convTool == "qctool")
        gensam <- .prepareInputFilesForImpute2_1(obj,wspace,rc)
    genFiles <- gensam$gen
    sampleFiles <- gensam$sam
    # These must be named BASE_chrZ.gen, BASE_chrZ.sam
    
    # Create the imputation intervals
    chunkList <- .makeImputationIntervals(map,intSize)
    
    # Impute2 path
    itool <- .getToolPath("impute")
    
    chrs <- as.character(unique(map$chromosome))
    disp("\n",.symbolBar("#",64))
    disp("External imputation on ",length(chrs)," chromosomes")
    disp(.symbolBar("#",64))
    
    # Pre create the chromosome directories
    chrDirs <- file.path(wspace,"chromosomes",chrs)
    toCreate <- chrDirs[!dir.exists(chrDirs)]
    for (d in toCreate)
        dir.create(d,recursive=TRUE,mode="0755",showWarnings=FALSE)
    #chrs <- chrs[!(chrs %in% c("1","10","11","12","13","14","15","16"))]
    genFiles <- genFiles[chrs]
    sampleFiles <- sampleFiles[chrs]
    null <- lapply(chrs,function(x) {
        disp("\n",.symbolBar("=",64))
        disp("Imputing on chromosome ",x)
        disp(.symbolBar("=",64))
        
        # Input and intervals
        g <- genFiles[[x]]
        ints <- chunkList[[x]]
        
        # Actual imputation on intervals in parallel
        .runImpute2(g,x,ints,refSpace,wspace,itool,rc)
    })
    
    # After finish, harvest interval results for each chromosome
    disp("\nImputation finished, re-merging imputation intervals")
    impList <- .summarizeImpute2Run(wspace) # Contains named lists
    impFiles <- impList$gen
    infoFiles <- impList$info
    
    # And convert back to PLINK
    .postProcessImpute2Files(impFiles,sampleFiles,rc=rc)
    
    # TODO: Urgent, need to add a step of parsing/cutting GEN files to get
    # the alleles as they are returned by gtool as 1, 2 instead of ACGT if
    # there are indels in the reference
    # Seems that we didn't have to implement this... Alleles are in the a0 and
    # a1 columns of INFO files...
    #disp("\nExtracting alleles from GEN files for later restoration")
    #alleleFiles <- .extractAllelesFromGen(impFiles)
    
    # At this point we must be having BASE_chrZ.imputed file triplets
    # These must be moved somewhere? Covnerted to GWASExperiment?
    # Read back to GWASExperiment
    disp("\nInterval merging finished, reading back to GWASExperiment")
    ibedFiles <- dir(file.path(wspace,"output"),pattern=".imputed.bed$",
        full.names=TRUE)
    pheno <- phenotypes(obj)
    objList <- cmclapply(ibedFiles,function(x) {
        disp("Prosessing ",x)
        x <- sub("\\.bed$","",x)
        y <- file.path(wspace,"chromosomes",paste0(basename(x),".info"))
        
        # It's possible that duplicate SNPs are created through the imputation
        # process... We remove them through SNP selection while reading PLINK
        selection <- list(samples=NULL,snps=NULL)
        bim <- read.delim(paste0(x,".bim"),header=FALSE)
        nodup <- which(!duplicated(bim[,seq_len(4)]))
        if (length(nodup) != nrow(bim))
            selection$snps <- nodup
        
        # GDS files will not be needed further though
        gwe <- importGWAS(
            input=x,
            phenos=pheno,
            backend=metadata(obj)$backend,
            genome=genome(obj),
            selection=selection,
            gdsfile=file.path(dirname(x),paste0(basename(x),".gds")),
            gdsOverwrite=FALSE
        )
        
        # INFO score should also be attached
        infoData <- read.delim(y)
        if (!is.null(selection$snps))
            infoData <- infoData[!duplicated(infoData[,c(2,3)]),,drop=FALSE]
        rownames(infoData) <- infoData$rs_id
        infoData <- infoData[rownames(gwe),]
        
        g <- gfeatures(gwe)
        # Correct alleles in A, C, G, T format and add info
        g$allele.1 <- infoData$a0
        g$allele.2 <- infoData$a1
        g$info <- infoData$info
        gfeatures(gwe) <- g
        
        return(gwe)
    },rc=rc)
    
    # Now we have to cbind and order
    disp("\nFinalizing...")
    impObj <- suppressWarnings(do.call("rbind",objList))
    tmp <- gfeatures(impObj)
    theOrder <- order(tmp$chromosome,tmp$position)
    # There is a warning about the size of SnpMatrix objects... OK
    impObj <- suppressWarnings(impObj[theOrder,])
    
    # We need to restore original workspace, without runId as it may have to be
    # deleted overall
    wspace <- wspaceBak
    switch(cleanup,
        none = {
            disp("All imputation pipeline output can be found at ",
                file.path(wspace,runId))
        },
        intermediate = {
            disp("Cleaning up temporary interval imputation files and other ",
                "intermediate files from ",wspace)
            .imputePartialCleanup(wspace,runId)
        },
        all = {
            disp("Cleaning up all imputation pipeline files!")
            unlink(wspace,recursive=TRUE,force=TRUE)
        }
    )
    
    disp("Done! The run id was ",runId,"\n")
    
    return(impObj)
}

#~ .extractAllelesFromGen <- function(gen) {
#~     ..cutCommand <- function(x) {
#~         o <- paste0(gsub("\\.gen$","",x),".als")
#~         return(paste0("cut -f2,4,5 -d' ' ",x," > ",o))
#~     }
    
#~     # Is cut command available? We may be in Windows
#~     useCut <- ifelse(Sys.which("cut")=="",FALSE,TRUE)
    
#~     if (useCut) {
#~         alsOut <- cmclapply(gen,function(x) {
#~             cmd <- ..cutCommand(x)
#~             disp("  extracting from ",x)
#~             disp("\nExecuting:\n",cmd,level="full")
#~             out <- tryCatch({
#~                 log <- capture.output({
#~                     system(cmd,intern=TRUE)
#~                 })
#~                 FALSE
#~             },error=function(e) {
#~                 message("Caught error: ",e$message)
#~                 return(TRUE)
#~             },finally="")
#~             return(out)
#~         },rc=rc)
#~     }
#~     else { # Vanilla R, will take a bit longer
#~         # Read a few lines from the first gen to get columns to remove
#~         tmp <- read.table(gen[1],nrow=10)
#~         colClasses <- c("NULL","character","NULL",rep("character",2),
#~             rep("NULL",ncol(tmp)-5))
#~         alsOut <- cmclapply(gen,function(x) {
#~             disp("  reading from ",x)
#~             out <- tryCatch({
#~                 als <- read.table(x,colClasses=colClasses,quote="")
#~                 write.table(als,file=paste0(gsub("\\.gen$","",x),".als"),
#~                     quote=FALSE,row.names=FALSE,col.names=FALSE)
#~                 FALSE
#~             },error=function(e) {
#~                 message("Caught error: ",e$message)
#~                 return(TRUE)
#~             },finally="")
#~             return(out)
#~         },rc=rc)
#~     }
    
#~     if (any(alsOut))
#~         stop("A problem occured during allele extraction from GEN files ",
#~             paste(gen[alsOut],collapse=", "),". Please check!")
    
#~     return(paste(gsub("\\.gen$","",gen),".als",sep=""))
#~ }

.imputePartialCleanup <- function(wspace,rid) {
    ..intMoveUp <- function(x,r) {
        if (length(x) > 0) {
            for (y in x) {
                to <- sub(paste0(r,"_"),"",y)
                to <- file.path(dirname(to),"..",basename(to))
                file.rename(from=y,to=to)
            }
        }
    }
    
    wspace <- file.path(wspace,rid)
    unlink(file.path(wspace,"input"),recursive=TRUE,force=TRUE)
    unlink(file.path(wspace,"chromosomes"),recursive=TRUE,force=TRUE)
    
    beds <- dir(file.path(wspace,"output"),
        #pattern=paste0("^",rid,".*.imputed.bed$"),full.names=TRUE)
        pattern=".imputed.bed$",full.names=TRUE)
    bims <- dir(file.path(wspace,"output"),
        pattern=".imputed.bim$",full.names=TRUE)
    fams <- dir(file.path(wspace,"output"),
        pattern=".imputed.fam$",full.names=TRUE)
    ..intMoveUp(beds,rid)
    ..intMoveUp(bims,rid)
    ..intMoveUp(fams,rid)
    unlink(file.path(wspace,"output"),recursive=TRUE,force=TRUE)
}

.prepareInputFilesForImpute2 <- function(obj,wspace,rc=NULL) {
    ..plinkToPedCommand <- function(x,p) {
        x <- gsub("\\.bed$","",x)
        return(paste(
           paste0(p," \\"),
           paste0("  --bfile ",x," \\"),
           "  --recode \\",
           paste0("  --out ",x),
           sep="\n"
        ))
    }
    
    ..gtoolToGenCommand <- function(x,g) {
        x <- gsub("\\.ped$","",x)
        return(paste(
           paste0(g," -P \\"),
           paste0("  --ped ",paste0(x,".ped")," \\"),
           paste0("  --map ",paste0(x,".map")," \\"),
           paste0("  --og ",paste0(x,".gen")," \\"),
           paste0("  --os ",paste0(x,".sample")),
           sep="\n"
        ))
    }
    
    gtool <- .getToolPath("gtool")
    plink <- .getToolPath("plink")
    
    # Export GWASExperiment as PLINK per chromosome
    disp("")
    inputDir <- file.path(wspace,"input")
    if (!dir.exists(inputDir))
        dir.create(inputDir,recursive=TRUE,mode="0755",showWarnings=FALSE)
    writePlink(obj,outBase=file.path(inputDir,"plink_impute2"),perChr=TRUE,
        overwrite=FALSE)
    
    # Convert to PED
    disp("\nConverting BED files to PED")
    bedFiles <- dir(inputDir,pattern=".bed$",full.names=TRUE)
    pedOut <- unlist(cmclapply(bedFiles,function(x,p) {
        dest <- sub(".bed$",".ped",x)
        if (!file.exists(dest)) {
            cmd <- ..plinkToPedCommand(x,p)
            disp("  converting ",x)
            disp("\nExecuting:\n",cmd,level="full")
            out <- tryCatch({
                log <- .formatSystemOutputForDisp(capture.output({
                    system(cmd,intern=TRUE)
                }))
                disp("\nPLINK output is:\n",level="full")
                disp(paste(log,collapse="\n"),"\n",level="full")
                FALSE
            },error=function(e) {
                message("Caught error: ",e$message)
                return(TRUE)
            },finally="")
        }
        else {
            disp("  file ",dest," already exists! Skipping...")
            return(FALSE)
        }
    },plink,rc=rc))
    
    # Conversion should be sucesfull...
    if (any(pedOut))
        stop("A problem occured during the conversion of BED files ",
            paste(bedFiles[pedOut],collapse=", "),". Please check!")
    
    pedFiles <- dir(inputDir,pattern="\\.ped$",full.names=TRUE)
    #mapFiles <- dir(wspace,pattern="\.map$")
    disp("\nConverting PED files to GEN")
    genOut <- unlist(cmclapply(pedFiles,function(x,g) {
        dest <- sub(".ped$",".gen",x)
        if (!file.exists(dest)) {
            cmd <- ..gtoolToGenCommand(x,g)
            disp("  converting ",x)
            disp("\nExecuting:\n",cmd,level="full")
            out <- tryCatch({
                log <- .formatSystemOutputForDisp(capture.output({
                    system(cmd,intern=TRUE)
                }))
                disp("\nPLINK output is:\n",level="full")
                disp(paste(log,collapse="\n"),"\n",level="full")
                FALSE
            },error=function(e) {
                message("Caught error: ",e$message)
                return(TRUE)
            },finally="")
        }
        else {
            disp("  file ",dest," already exists! Skipping...")
            return(FALSE)
        }
    },gtool,rc=rc))
    
    if (any(genOut))
        stop("A problem occured during the generation of GEN files ",
            paste(pedFiles[genOut],collapse=", "),". Please check!")
    
    genFiles <- dir(inputDir,pattern="\\.gen$")
    samFiles <- dir(inputDir,pattern="\\.sample$")
    
    # dir-ing is unstable... We must name these vectors with chromosomes
    # extracted from their names...
    names(genFiles) <- unlist(lapply(strsplit(genFiles,"_"),function(s) {
        sub(".gen","",sub("chr","",s[3]))
    }))
    names(samFiles) <- unlist(lapply(strsplit(samFiles,"_"),function(s) {
        sub(".sample","",sub("chr","",s[3]))
    }))
    
    return(list(gen=genFiles,sam=samFiles))
}

# Convert back to PLINK triplets with GTOOL and PLINK
# Or just input the workspace? As we need the chromosomes too
.postProcessImpute2Files <- function(impFiles,samFiles,rc=NULL) {
    ..plinkToBedCommand <- function(x,p) {
        x <- gsub("\\.ped$","",x)
        return(paste(
           paste0(p," \\"),
           paste0("  --file ",x," \\"),
           "  --make-bed \\",
           paste0("  --out ",x),
           sep="\n"
        ))
    }
    
    ..gtoolToPedCommand <- function(x,y,g) {
        x <- gsub("\\.imputed\\.gen$","",x)
        xs <- file.path(dirname(x),"..","input",basename(x))
        xp <- file.path(dirname(x),"..","output",basename(x))
        return(paste(
           paste0(g," -G \\"),
           paste0("  --g ",paste0(x,".imputed.gen")," \\"),
           paste0("  --s ",paste0(xs,".sample")," \\"),
           paste0("  --ped ",paste0(xp,".imputed.ped")," \\"),
           paste0("  --map ",paste0(xp,".imputed.map")," \\"),
           paste0("  --chr ",y),
           #paste0("  --chr ",y," \\"),
           #"  --snp ",
           sep="\n"
        ))
    }
    
    gtool <- .getToolPath("gtool")
    plink <- .getToolPath("plink")
    
    outputDir <- file.path(dirname(impFiles[[1]]),"..","output")
    if (!dir.exists(outputDir))
        dir.create(outputDir,recursive=TRUE,mode="0755",showWarnings=FALSE)
    
    # Convert GEN impute files to PED
    disp("\nConverting GEN files to PED")
    pedOut <- unlist(cmclapply(names(impFiles),function(n,g,D) {
        x <- D[[n]]
        dest <- file.path(dirname(x),"..","output",paste0(sub(".gen$",".ped",
            basename(x))))
        if (!file.exists(dest)) {
            cmd <- ..gtoolToPedCommand(x,n,g)
            disp("  converting ",x)
            disp("\nExecuting:\n",cmd,level="full")
            out <- tryCatch({
                log <- .formatSystemOutputForDisp(capture.output({
                    system(cmd,intern=TRUE)
                }))
                disp("\nPLINK output is:\n",level="full")
                disp(paste(log,collapse="\n"),"\n",level="full")
                FALSE
            },error=function(e) {
                message("Caught error: ",e$message)
                return(TRUE)
            },finally="")
        }
        else {
            disp("  file ",dest," already exists! Skipping...")
            return(FALSE)
        }
    },gtool,impFiles,rc=rc))
    
    if (any(pedOut))
        stop("A problem occured during the generation of GEN files ",
            paste(impFiles[pedOut],collapse=", "),". Please check!")
    
    # Convert PED to BED
    disp("\nConverting PED files to BED")
    pedFiles <- dir(file.path(dirname(impFiles[[1]]),"..","output"),
        pattern="\\.imputed\\.ped$",full.names=TRUE)
    bedOut <- unlist(cmclapply(pedFiles,function(x,p) {
        dest <- sub(".ped$",".bed",x)
        if (!file.exists(dest)) {
            cmd <- ..plinkToBedCommand(x,p)
            disp("  converting ",x)
            disp("\nExecuting:\n",cmd,level="full")
            out <- tryCatch({
                log <- .formatSystemOutputForDisp(capture.output({
                    system(cmd,intern=TRUE)
                }))
                disp("\nPLINK output is:\n",level="full")
                disp(paste(log,collapse="\n"),"\n",level="full")
                FALSE
            },error=function(e) {
                message("Caught error: ",e$message)
                return(TRUE)
            },finally="")
        }
        else {
            disp("  file ",dest," already exists! Skipping...")
            return(FALSE)
        }
    },plink,rc=rc))
    
    # Conversion should be sucesfull...
    if (any(pedOut))
        stop("A problem occured during the conversion of PED files ",
            paste(pedFiles[bedOut],collapse=", "),". Please check!")
}

# This function expects to find a directory with one subdir for each chromosome
# and in that subdir, results for each interval. We simply cat files.
.summarizeImpute2Run <- function(wspace,rc=NULL) {
    # We expect only to find a dir per chromosome, nothing else!
    chrDirs <- list.dirs(file.path(wspace,"chromosomes"),full.names=TRUE,
        recursive=FALSE)
    # The outout filename must be the part before ___interval___
    # We need to summarize both genotype and info files for INFO score
    sumGens <- cmclapply(chrDirs,function(x) {
        disp("Summarizing imputation files for chromosome ",basename(x))
        filesToCat <- dir(x,pattern=".gen$",full.names=TRUE)
        mainName <- strsplit(basename(filesToCat[1]),"___")[[1]][1]
        mainFile <- file.path(wspace,"chromosomes",
            paste0(mainName,".imputed.gen"))
        if (!file.exists(mainFile))
            file.append(mainFile,filesToCat)
        else
            disp("  file ",mainFile," exists! Skipping...")
        return(mainFile)
    },rc=rc)
    
    sumInfos <- lapply(chrDirs,function(x) {
        disp("Summarizing info files for chromosome ",basename(x))
        filesToCat <- dir(x,pattern=".gen_info$",full.names=TRUE)
        mainName <- strsplit(basename(filesToCat[1]),"___")[[1]][1]
        mainFile <- file.path(wspace,"chromosomes",
            paste0(mainName,".imputed.info"))
        if (!file.exists(mainFile)) {
            infoList <- cmclapply(filesToCat,function(z) {
                return(read.table(z,header=TRUE,quote=""))  
            },rc=rc)
            toWrite <- do.call("rbind",infoList)
            write.table(toWrite,file=mainFile,sep="\t",quote=FALSE,
                row.names=FALSE)
        }
        else
            disp("  file ",mainFile," exists! Skipping...")
        return(mainFile)
    })
    
    names(sumGens) <- names(sumInfos) <- basename(chrDirs)
    return(list(gen=sumGens,info=sumInfos))
}

# Impute expects the following patterns in a specific folder:
# -m genetic_map_$CHR_combined_b37.txt
# -h 1000GP_Phase3_$CHR.hap.gz 
# -l 1000GP_Phase3_$CHR.legend.gz
# chr is always present, therefore it's better to make the BED files using
# writePlink
# impute2 produces
# thiseas_nafld_common_phased_chr18.impute2
# thiseas_nafld_common_phased_chr18.impute2_info
# thiseas_nafld_common_phased_chr18.impute2_info_by_sample
# From these, thiseas_nafld_common_phased_chr18.impute2 with the previous 
# sample file is needed to go back to plink
# thiseas_nafld_common_phased_chr18.impute2_info_by_sample
.runImpute2 <- function(g,chr,ints,refSpace,wspace,exec,rc) {
    # Expected files
    mFile <- file.path(refSpace,paste0("genetic_map_chr",chr,
        "_combined_b37.txt"))
    hFile <- file.path(refSpace,paste0("1000GP_Phase3_chr",chr,".hap.gz"))
    lFile <- file.path(refSpace,paste0("1000GP_Phase3_chr",chr,".legend.gz"))
    
    # Output file base and dir per chromosome
    chrDir <- file.path(wspace,"chromosomes",chr)
    #if (!dir.exists(chrDir))
    #    dir.create(chrDir,recursive=TRUE,mode="0755",showWarnings=FALSE)
    oFileBase <- file.path(chrDir,sub("(.*)\\..*$","\\1",basename(g)))
    g <- file.path(wspace,"input",basename(g))
    
    # Create a list of intervals for parallel
    ints <- as.list(as.data.frame(t(ints)))
    disp("---------- ",length(ints)," intervals")
    
    # Parallely impute
    fails <- unlist(cmclapply(ints,function(x) {
        intext <- paste0(format(x[1],scientific=FALSE),"-",
            format(x[2],scientific=FALSE))
        oFile <- paste0(oFileBase,"___",intext,"___.gen")
        
        disp("Imputing for interval ",intext," output at ",oFile)
        
        #if (!file.exists(oFile)) { # For crash restart support later
        if (!.intervalComplete(intext,chrDir)) {
            cmd <- paste(
               paste0(exec," \\"),
               paste0("  -m ",mFile," \\"),
               paste0("  -h ",hFile," \\"),
               paste0("  -l ",lFile," \\"),
               paste0("  -g ",g," \\"),
               paste0("  -int ",x[1]," ",x[2]," \\"),
               paste0("  -Ne 20000 \\"),
               paste0("  -align_by_maf_g \\"),
               paste0("  -filt_rules_l 'EUR==0' \\"),
               paste0("  -o ",oFile),
               sep="\n"
            )
            
            disp("\nExecuting:\n",cmd,level="full")
            
            # Init JSON for later completion checking here
            .recordIntervalProgress(intext,chrDir)
            
            out <- tryCatch({
                log <- .formatSystemOutputForDisp(capture.output({
                    system(cmd,intern=TRUE)
                }))
                disp("\nIMPUTE2 output is:\n",level="full")
                disp(paste(log,collapse="\n"),"\n",level="full")
                .recordIntervalProgress(intext,chrDir,done=TRUE)
                FALSE
            },error=function(e) {
                message("Caught error: ",e$message)
                return(TRUE)
            },finally="")
        }
        else {
            disp("  result for interval ",intext," exists! Skipping...")
            out <- FALSE
        }
        
        return(out)
    },rc=rc))

    if (any(fails)) {
        fa <- which(fails)
        faInt <- ints[fa]
        faTxt <- unlist(lapply(faInt,function(x) {
            return(paste0(format(x[1],scientific=FALSE),"-",
                format(x[2],scientific=FALSE)))
        }))
        warning("A problem occured during imputation for interval(s) ",
            paste(faTxt,collapse=", "),". Please check!",immediate.=TRUE)
    }
}

# x is a gwe or gfeatures
# There must be a minimum number of SNPs in each interval...
.makeImputationIntervals <- function(x,size=1e+6,ms=5) {
    # Check x
    if (!is(x,"GWASExperiment") && !is(x,"DataFrame") && !is.data.frame(x))
        stop("Input object must be a GWASExperiment or its gfeatures!")
    if (is(x,"GWASExperiment"))
        x <- gfeatures(x)
        
    # Is x position sorted? Check with a small sample to avoid full re-ordering
    xx <- x[sort(sample(nrow(x),100)),]
    yy <- xx[order(xx$chromosome,xx$position),]
    if (!identical(xx,yy))
        x <- x[order(xx$chromosome,xx$position),]
    
    # Split positions per chromosome
    x <- as.data.frame(x[,c("chromosome","position")])
    #S <- split(x$position,x$chromosome)
    S <- split(x,x$chromosome)
    
    # Firstly create a GRanges with smallest and largest SNP position
    tmp <- GRanges(do.call("rbind",lapply(names(S),function(n,Y) {
        y <- Y[[n]]$position
        return(data.frame(chromosome=n,start=min(y),end=max(y)))        
    },S)))
    
    # Then calculate number of chunks of length ~size according to size
    nc <- vapply(S,function(x) {
        return(round((max(x$position,na.rm=TRUE) - 
            min(x$position,na.rm=TRUE))/size))
    },numeric(1))
    
    # Now tile
    ints <- tile(tmp,nc)
    
    ## We should now overlap with SNP locations... if <2 merge
    #ints <- .correctTheIntervals(ints,x)
    
    return(lapply(ints,function(x) {
        return(as.data.frame(unname(x))[,c("start","end")])
    }))
    
}

.correctTheIntervals <- function(ints,x,ms) {
    # Make GPos object from 
    xg <- x
    xg$end <- x$position
    names(xg)[2] <- "start"
    xg <- GPos(xg)
    
    # Split it so as to apply find overlaps
    xgl <- split(xg,seqnames(xg))
    
    # Find overlaps. Chromosome existence should not be a problem since ints
    # is constructed from gfeatures.
    hits <- lapply(names(ints),function(n,a,b) {
        findOverlaps(a[[n]],b[[n]])
    },ints,xgl)
    names(hits) <- names(ints)
    
    # Locate empty intervals and drop them
    tables <- lapply(hits,function(y) {
        return(table(queryHits(y)))
    })
    indsToKeep <- lapply(tables,function(y) {
        return(as.numeric(names(y)))
    })
    ints <- lapply(names(ints),function(n,a,b) {
        return(a[[n]][b[[n]]])
    },ints,indsToKeep)
    names(ints) <- names(hits)
    
    # Now we should merge intervals with number of snps < ms
    lapply(names(ints),function(n,a,b) {
        # Any interval has < ms SNPs? If not happy
        if (!any(b[[n]] < ms))
            return(a[[n]])
        
        # Else reduce... merge with the neighboring interval with the less SNPs
        small <- which(b[[n]] < ms)
        #if (any(small == 1)) # merge with 2
        #if (any(small == length(b[[n]]))) # merge with n-1
        # etc.
        
    },ints,tables)
}

.recordIntervalProgress <- function(int,chrDir,done=FALSE) {
    progFile <- file.path(chrDir,paste0(".",int,".json"))
    if (!file.exists(progFile)) { # Initial write
        prog <- list(done=FALSE)
        write_json(prog,path=progFile,auto_unbox=TRUE,pretty=TRUE)
    }
    else { # Exists, update
        curr <- fromJSON(progFile)
        curr$done <- done
        write_json(curr,path=progFile,auto_unbox=TRUE,pretty=TRUE)
    }
}

.intervalComplete <- function(int,chrDir) {
    progFile <- file.path(chrDir,paste0(".",int,".json"))
    if (!file.exists(progFile))
        return(FALSE)
    else {
        curr <- fromJSON(progFile)
        return(curr$done)
    }
}

#~ impute2 \
#~   -m /impute/1000GP_Phase3/genetic_map_chr1_combined_b37.txt \
#~   -h /impute/1000GP_Phase3/1000GP_Phase3_chr1.hap.gz \
#~   -l /1000GP_Phase3/1000GP_Phase3_chr1.legend.gz \
#~   -g ../thiseas_nafld_common_1.gen \
#~   -align_by_maf_g \
#~   -int 1 100000 \
#~   -Ne 20000 \
#~   -filt_rules_l 'EUR==0' \
#~   -o ./thiseas_nafld_common_phased_chr1.impute2

# The path is the place where the file will be unziped so a directory called
# 1000GP_Phase3 will be created there. It will not be stored in system's
# variables like software, it must be provided during imputation.
download1000GP3 <- function(path=NULL) {
    # Address from Oxford
    remote <- "https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.tgz"
    
    # Manage path
    if (!dir.exists(path)) {
        disp("1000GP imputation files home directory created at ",path)
        dir.create(path,recursive=TRUE,mode="0755",showWarnings=FALSE)
    }
    else
        disp("1000GP imputation files home directory found at ",path)
    
    # Download and uncompress
    disp("Downlading 1000GP Phase 3 imputation files...")
    des <- file.path(path,"1000GP_Phase3.tgz")
    download.file(remote,des)
    
    # Unzip the archive
    disp("Uncompressing...")
    untar(des,exdir=path)
    
    # Delete compressed
    disp("Deleting downloaded archive...")
    unlink(des,force=TRUE)
    
    # Test whether all properly fetched and installed
    disp("Verifying...")
    if (dir.exists(file.path(path,"1000GP_Phase3"))
        && file.exists(file.path(path,"1000GP_Phase3",
        "1000GP_Phase3_chr1.hap.gz"))) {
            disp("1000GP Phase3 imputation files seem to be in place!")
            return(file.path(path,"1000GP_Phase3"))
        }
    else {
        disp("1000GP Phase3 imputation files were not downloaded! Please ",
            "check manually")
        return(NULL)
    }
}

#~ .imputePartialCleanup <- function(wspace,ri) {
#~     chrdirs <- list.dirs(wspace,full.names=TRUE,recursive=FALSE)
#~     unlink(chrdirs,recursive=TRUE,force=TRUE)
    
#~     # Rename final .ped, .bim, .fam to remove run id
#~     beds <- dir(wspace,pattern=paste0("^",ri,".*.imputed.bed$"),full.names=TRUE)
#~     bims <- dir(wspace,pattern=paste0("^",ri,".*.imputed.bim$"),full.names=TRUE)
#~     fams <- dir(wspace,pattern=paste0("^",ri,".*.imputed.fam$"),full.names=TRUE)
#~     if (length(beds) > 0) {
#~         for (be in beds)
#~             file.rename(from=be,to=sub(paste0(ri,"_"),"",be))
#~     }
#~     if (length(bims) > 0) {
#~         for (bi in bims)
#~             file.rename(from=bi,to=sub(paste0(ri,"_"),"",bi))
#~     }
#~     if (length(fams) > 0) {
#~         for (fa in fams)
#~             file.rename(from=fa,to=sub(paste0(ri,"_"),"",fa))
#~     }
    
#~     tmps <- dir(wspace,pattern=paste0("^",ri),full.names=TRUE)
#~     if (length(tmps) > 0)
#~         unlink(tmps,recursive=TRUE,force=TRUE)   
#~ }

.prepareInputFilesForImpute2_1 <- function(obj,wspace,rc=NULL) {
    ..qctoolToGenCommand <- function(x,g) {
        x <- gsub("\\.bed$","",x)
        return(paste(
           paste0(g," \\"),
           paste0("  -g ",paste0(x,".bed")," \\"),
           paste0("  -og ",paste0(x,".tmp.gen")," \\"),
           paste0("  -os ",paste0(x,".sample")),
           sep="\n"
        ))
    }
    
    qctool <- .getToolPath("qctool")
    
    # Export GWASExperiment as PLINK per chromosome
    disp("")
    inputDir <- file.path(wspace,"input")
    if (!dir.exists(inputDir))
        dir.create(inputDir,recursive=TRUE,mode="0755",showWarnings=FALSE)
    writePlink(obj,outBase=file.path(inputDir,"plink_impute2"),perChr=TRUE,
        overwrite=FALSE)
    
    bedFiles <- dir(inputDir,pattern=".bed$",full.names=TRUE)
    disp("\nConverting BED files to GEN")
    genOut <- unlist(cmclapply(bedFiles,function(x,g) {
        dest <- sub(".bed$",".gen",x)
        if (!file.exists(dest)) {
            cmd <- ..qctoolToGenCommand(x,g)
            disp("  converting ",x)
            disp("\nExecuting:\n",cmd,level="full")
            out <- tryCatch({
                log <- .formatSystemOutputForDisp(capture.output({
                    system(cmd,intern=TRUE,ignore.stdout=TRUE,
                        ignore.stderr=TRUE)
                }))
                disp("\nQCTOOL output is:\n",level="full")
                disp(paste(log,collapse="\n"),"\n",level="full")
                
                xx <- sub(".bed$",".tmp.gen",x)
                y <- read.table(xx)
                xx <- gsub("\\.tmp\\.gen$","",xx)
                y <- y[,-3,drop=FALSE]
                xxx <- paste0(xx,".gen")
                write.table(y,file=xxx,row.names=FALSE,col.names=FALSE,
                    quote=FALSE)
                unlink(paste0(xx,".tmp.gen"),recursive=TRUE,force=TRUE)
                FALSE
            },error=function(e) {
                message("Caught error: ",e$message)
                return(TRUE)
            },finally="")
        }
        else {
            disp("  file ",dest," already exists! Skipping...")
            return(FALSE)
        }
    },qctool,rc=rc))
    
    if (any(genOut))
        stop("A problem occured during the generation of GEN files ",
            paste(bedFiles[genOut],collapse=", "),". Please check!")
    
    genFiles <- dir(inputDir,pattern="\\.gen$")
    #genFiles <- dir(inputDir,pattern="\\.gen$",full.names=TRUE)
    samFiles <- dir(inputDir,pattern="\\.sample$")
    
    ## Remove 3rd column from genFiles
    #genFiles <- basename(unlist(cmclapply(genFiles,function(x) {
    #    disp("  correcting format for ",x)
    #    y <- read.table(x)
    #    x <- gsub("\\.gentmp$","",x)
    #    y <- y[,-3,drop=FALSE]
    #    x <- paste0(x,".gen")
    #    write.table(y,file=x,row.names=FALSE,col.names=FALSE)
    #    return(x)
    #},rc=rc)))
    
    # dir-ing is unstable... We must name these vectors with chromosomes
    # extracted from their names...
    names(genFiles) <- unlist(lapply(strsplit(genFiles,"_"),function(s) {
        sub(".gen","",sub("chr","",s[3]))
    }))
    names(samFiles) <- unlist(lapply(strsplit(samFiles,"_"),function(s) {
        sub(".sample","",sub("chr","",s[3]))
    }))
    
    return(list(gen=genFiles,sam=samFiles))
}
