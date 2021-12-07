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
    continue=FALSE,cleanup=c("none","intermediate","all"),rc=NULL) {
    if (!.toolAvailable("gtool"))
        stop("GTOOL program not found in the system!")
    if (!.toolAvailable("impute"))
        stop("IMPUTE2 program not found in the system!")
    
    cleanup <- cleanup[1]
    .checkTextArgs("Cleanup option (cleanup)",cleanup,
        c("none","intermediate","all"),multiarg=FALSE)
    
    # Validate workspace
    wspace <- .validateWorkspacePath(wspace,"impute")
    
    # Check refSpace contains 1000GP files
    # ...
    
    # Check chromosomes are there
    map <- gfeatures(obj)
    if (!("chromosome" %in% names(map))) # Trying to impute with what?
        stop("Cannot extend the dataset without location info!")
    
    # We need a runid for potentially working in the same space
    runid <- .randomString()
    
    # Prepare the main input files to impute per chromosome
    gensam <- .prepareInputFilesForImpute2(obj,runid,wspace,rc)
    genFiles <- gensam$gen
    sampleFiles <- gensam$sam
    # These must be named BASE_chrZ.gen, BASE_chrZ.sam
    
    # writePlink writes with split map chromosome order
    S <- split(map,map$chromosome)
    chrs <- names(genFiles) <- names(sampleFiles) <- names(S)
    
    # Create the imputation intervals
    chunkList <- .makeImputationIntervals(map,intSize)
    
    # Impute2 path
    itool <- .getToolPath("impute")
    
    disp("\n",.symbolBar("#",64))
    disp("External imputation on ",length(chrs)," chromosomes")
    disp(.symbolBar("#",64))
    devNull <- lapply(chrs,function(x) {
        disp(.symbolBar("=",64))
        disp("\nImputing on chromosome ",x)
        disp(.symbolBar("=",64))
        
        # Input and intervals
        g <- genFiles[[x]]
        ints <- chunkList[[x]]
        
        # Actual imputation on intervals in parallel
        .runImpute2(g,x,ints,refSpace,wspace,itool,rc)
    })
    
    # After finish, harvest interval results for each chromosome
    disp("\nImputation finished, re-merging imputation intervals")
    impFiles <- .summarizeImpute2Run(wspace) # It's a named list!
    
    # And convert back to PLINK
    .postProcessImpute2Files(impFiles,sampleFiles,rc=NULL)
    
    # At this point we must be having BASE_chrZ.imputed file triplets
    # These must be moved somewhere? Covnerted to GWASExperiment?
    # Read back to GWASExperiment
    disp("\nInterval merging finished, reading back to GWASExperiment")
    ibedFiles <- dir(wspace,pattern="\\.imputed\\.bed$",full.names=TRUE)
    pheno <- phenotypes(obj)
    objList <- cmclapply(ibedFiles,function(x) {
        x <- sub("\\.bed$","",x)
        # GDS files will not be needed further though
        return(importGWAS(
            input=x,
            phenos=pheno,
            backend=metadata(obj)$backend,
            genome=genome(obj),
            gdsfile=file.path(dirname(x),paste0(basename(x),".gds")),
            gdsOverwrite=FALSE
        ))
    },rc=rc)
    
    # Now we have to cbind and order
    disp("\nFinalizing...")
    impObj <- do.call("rbind",objList)
    tmp <- gfeatures(impObj)
    theOrder <- order(tmp$chromosome,tmp$position)
    impObj <- impObj[theOrder,]
    
    switch(cleanup,
        none = {
            disp("All imputation pipeline output can be found at ",wspace)
        },
        intermediate = {
            disp("Cleaning up temporary interval imputation files and other ",
                "intermediate files from ",wspace)
            .imputePartialCleanup(wspace,runid)
        },
        all = {
            disp("Cleaning up all imputation pipeline files!")
            unlink(wspace,recursive=TRUE,force=TRUE)
        }
    )
    
    disp("Done!\n")
    
    return(impObj)
}

.imputePartialCleanup <- function(wspace,ri) {
    chrdirs <- list.dirs(wspace,full.names=TRUE,recursive=FALSE)
    unlink(chrdirs,recursive=TRUE,force=TRUE)
    
    # Rename final .ped, .bim, .fam to remove run id
    beds <- dir(wspace,pattern=paste0("^",ri,".*.imputed.bed$"),full.names=TRUE)
    bims <- dir(wspace,pattern=paste0("^",ri,".*.imputed.bim$"),full.names=TRUE)
    fams <- dir(wspace,pattern=paste0("^",ri,".*.imputed.fam$"),full.names=TRUE)
    if (length(beds) > 0) {
        for (be in beds)
            file.rename(from=be,to=sub(paste0(ri,"_"),"",be))
    }
    if (length(bims) > 0) {
        for (bi in bims)
            file.rename(from=bi,to=sub(paste0(ri,"_"),"",bi))
    }
    if (length(fams) > 0) {
        for (fa in fams)
            file.rename(from=fa,to=sub(paste0(ri,"_"),"",fa))
    }
    
    tmps <- dir(wspace,pattern=paste0("^",ri),full.names=TRUE)
    if (length(tmps) > 0)
        unlink(tmps,recursive=TRUE,force=TRUE)
}

.prepareInputFilesForImpute2 <- function(obj,runid,wspace,rc=NULL) {
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
    writePlink(obj,outBase=file.path(wspace,paste0(runid,"_plink_impute2")),
        perChr=TRUE)
    
    # Convert to PED
    disp("\nConverting BED files to PED")
    bedFiles <- dir(wspace,pattern="\\.bed$",full.names=TRUE)
    pedOut <- unlist(cmclapply(bedFiles,function(x,p) {
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
    },plink,rc=rc))
    
    # Conversion should be sucesfull...
    if (any(pedOut))
        stop("A problem occured during the conversion of BED files ",
            paste(bedFiles[pedOut],collapse=", "),". Please check!")
    
    pedFiles <- dir(wspace,pattern="\\.ped$",full.names=TRUE)
    #mapFiles <- dir(wspace,pattern="\.map$")
    disp("\nConverting PED files to GEN")
    genOut <- unlist(cmclapply(pedFiles,function(x,g) {
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
    },gtool,rc=rc))
    
    if (any(genOut))
        stop("A problem occured during the generation of GEN files ",
            paste(pedFiles[genOut],collapse=", "),". Please check!")
    
    genFiles <- dir(wspace,pattern="\\.gen$")
    samFiles <- dir(wspace,pattern="\\.sample$")
    
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
        return(paste(
           paste0(g," -G \\"),
           paste0("  --g ",paste0(x,".imputed.gen")," \\"),
           paste0("  --s ",paste0(x,".sample")," \\"),
           paste0("  --ped ",paste0(x,".imputed.ped")," \\"),
           paste0("  --map ",paste0(x,".imputed.map")," \\"),
           paste0("  --chr ",y),
           sep="\n"
        ))
    }
    
    gtool <- .getToolPath("gtool")
    plink <- .getToolPath("plink")
    
    # Convert GEN impute files to PED
    disp("\nConverting GEN files to PED")
    pedOut <- unlist(cmclapply(names(impFiles),function(n,g,D) {
        x <- D[[n]]
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
    },gtool,impFiles,rc=rc))
    
    if (any(pedOut))
        stop("A problem occured during the generation of GEN files ",
            paste(impFiles[pedOut],collapse=", "),". Please check!")
    
    # Convert PED to BED
    disp("\nConverting PED files to BED")
    pedFiles <- dir(dirname(impFiles[[1]]),pattern="\\.imputed\\.ped$",
        full.names=TRUE)
    bedOut <- unlist(cmclapply(pedFiles,function(x,p) {
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
    },plink,rc=rc))
    
    # Conversion should be sucesfull...
    if (any(pedOut))
        stop("A problem occured during the conversion of PED files ",
            paste(pedFiles[bedOut],collapse=", "),". Please check!")
}

# This function expects to find a directory with one subdir for each chromosome
# and in that subdir, results for each interval. We simply cat files.
.summarizeImpute2Run <- function(wspace) {
    # We expect only to find a dir per chromosome, nothing else!
    chrDirs <- list.dirs(wspace,full.names=TRUE,recursive=FALSE)
    # The outout filename must be the part before ___interval___
    sumGens <- cmclapply(chrDirs,function(x) {
        disp("Summarizing imputation files for chromosome ",x)
        filesToCat <- dir(x,pattern=".gen$",full.names=TRUE)
        mainName <- strsplit(basename(filesToCat[1]),"___")[[1]][1]
        mainFile <- file.path(wspace,paste0(mainName,".imputed.gen"))
        file.append(mainFile,filesToCat)
        return(mainFile)
    },rc=rc)
    names(sumGens) <- basename(chrDirs)
    return(sumGens) # Or may later just dir them to proceed
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
    chrDir <- file.path(wspace,chr)
    if (!dir.exists(chrDir))
        dir.create(chrDir,recursive=TRUE,mode="0755",showWarnings=FALSE)
    oFileBase <- file.path(chrDir,sub("(.*)\\..*$","\\1",basename(g)))
    g <- file.path(wspace,basename(g))
    
    # Create a list of intervals for parallel
    ints <- as.list(as.data.frame(t(ints)))
    disp("---------- ",length(ints)," intervals")
    
    # Parallely impute
    fails <- unlist(cmclapply(ints,function(x) {
        intext <- paste0(format(x[1],scientific=FALSE),"-",
            format(x[2],scientific=FALSE))
        oFile <- paste0(oFileBase,"___",intext,"___.gen")
        
        disp("Imputing for interval ",intext," output at ",oFile)
        
        if (!file.exists(oFile)) { # For crash restart support later
            cmd <- paste(
               paste0(exec," -m \\"),
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
    },rc=rc))

    if (any(fails)) {
        fa <- which(fails)
        faInt <- ints[fa]
        faTxt <- unlist(lapply(faInt,function(x) {
            return(paste0(format(x[1],scientific=FALSE),"-",
                format(x[2],scientific=FALSE)))
        }))
        warning("A problem occured during imputation for interval(s) ",
            past(faTxt,collapse=", "),". Please check!",immediate.=TRUE)
    }
}

# x is a gwe or gfeatures
.makeImputationIntervals <- function(x,size=1e+6) {
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
    S <- split(x$position,x$chromosome)
    
    # Firstly create a GRanges with smallest and largest SNP position
    tmp <- GRanges(do.call("rbind",lapply(names(S),function(n,Y) {
        y <- Y[[n]]
        return(data.frame(chromosome=n,start=min(y),end=max(y)))        
    },S)))
    
    # Then calculate number of chunks of length ~size according to size
    nc <- vapply(S,function(x) {
        return(round((max(x,na.rm=TRUE)-min(x,na.rm=TRUE))/size))
    },numeric(1))
    
    # Now tile
    ints <- tile(tmp,nc)
    return(lapply(ints,function(x) {
        return(as.data.frame(unname(x))[,c("start","end")])
    }))
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
