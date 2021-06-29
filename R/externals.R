# This will download in default paths... If another path is desired, use 
# individual download functions directly
downloadExternalTools <- function(tools=c("snptest","plink","prsice"),
    ver=getDefaults("externalTools"),path=NULL) {
    .checkTextArgs("tools",tools,c("snptest","plink","prsice"),multiarg=TRUE)
    
    if ("snptest" %in% tools)
        downloadSnptest(ver=ver$snptest,path=path)
    if ("plink" %in% tools)
        downloadPlink(ver=ver$plink,path=path)
    if ("prsice" %in% tools)
        downloadPrsice(ver=ver$prsice,path=path)
}

downloadSnptest <- function(ver=c("v2.5.6","v2.5.4","v2.5.2"),path=NULL) {
    # This downloads the latest version by default... Might not work for older
    # Linux operating systems
    ver <- ver[1]
    .checkTextArgs("SNPTEST version (ver)",ver,c("v2.5.6","v2.5.4","v2.5.2"),
        multiarg=FALSE)
    
    base <- "http://www.well.ox.ac.uk/~gav/resources/"
    if (ver == "v2.5.6") {
        ext <- "snptest_v2.5.6_CentOS_Linux7.8-x86_64_dynamic.tgz"
        exec <- "snptest_v2.5.6"
    }
    else if (ver == "v2.5.4") {
        ext <- "snptest_v2.5.4-beta3_linux_x86_64_static.tgz"
        exec <- "snptest_v2.5.4-beta3"
    }
    else if (ver == "v2.5.2") {
        ext <- "snptest_v2.5.2_linux_x86_64_static.tgz"
        exec <- "snptest_v2.5.2"
    }
        
    src <- paste0(base,ext)
    unext <- sub(".tgz","",ext)
    sver <- substring(ver,2,nchar(ver))
    defPath <- .createToolPath("snptest",sver,"SNPTEST",path)
    
    # Download and uncompress
    disp("Preparing to download SNPTEST...")
    des <- file.path(tempdir(),ext)
    download.file(src,des)
    
    # Move and rename the executable to the common style
    disp("Installing...")
    untar(des,exdir=defPath)
    tomove <- dir(file.path(defPath,unext))
    for (f in tomove)
        file.rename(from=file.path(defPath,unext,f),to=file.path(defPath,f))
    unlink(file.path(defPath,unext),recursive=TRUE)
    file.rename(from=file.path(defPath,exec),to=file.path(defPath,"snptest"))
    Sys.chmod(file.path(defPath,"snptest"),"0755")
    
    # Test whether all properly fetched and installed
    disp("Verifying...")
    cmd <- paste0(file.path(defPath,"snptest")," -help")
    disp("Trying the command: ",cmd)
    v <- system(cmd,ignore.stdout=TRUE,ignore.stderr=TRUE)
    if (v == 0) {
        disp("SNPTEST seems to be working! Updating environment...")
        .EXTERNALS[["snptest"]]$dir <- defPath
        .EXTERNALS[["snptest"]]$exec <- "snptest"
        .EXTERNALS[["snptest"]]$version <- .tryNumerize(sver)
    }
    else
        disp("SNPTEST does NOT seem to be working! Please install manually ",
            "and make it accessible through the system's path or try another ",
            "version!")
}

downloadPlink <- function(ver=c("v1.90"),path=NULL) {
    # This downloads the latest version by default... Might not work for older
    # Linux operating systems
    ver <- ver[1]
    .checkTextArgs("PLINK version (ver)",ver,c("v1.90"),multiarg=FALSE)
    
    base <- "https://s3.amazonaws.com/plink1-assets/"
    if (ver == "v1.90") {
        ext <- "plink_linux_x86_64_20210606.zip"
        exec <- "plink"
    }
        
    src <- paste0(base,ext)
    sver <- substring(ver,2,nchar(ver))
    defPath <- .createToolPath("plink",sver,"PLINK",path)
    
    # Download and uncompress
    disp("Preparing to download PLINK...")
    des <- file.path(tempdir(),ext)
    download.file(src,des)
    
    # Move and rename the executable to the common style
    disp("Installing...")
    unzip(des,exdir=defPath)
    Sys.chmod(file.path(defPath,"plink"),"0755")
    
    # Test whether all properly fetched and installed
    disp("Verifying...")
    cmd <- paste0(file.path(defPath,"plink")," --help")
    disp("Trying the command: ",cmd)
    v <- system(cmd,ignore.stdout=TRUE,ignore.stderr=TRUE)
    if (v == 0) {
        disp("PLINK seems to be working! Updating environment...")
        .EXTERNALS[["plink"]]$dir <- defPath
        .EXTERNALS[["plink"]]$exec <- "plink"
        .EXTERNALS[["plink"]]$version <- .tryNumerize(sver)
    }
    else
        disp("PLINK does NOT seem to be working! Please install manually ",
            "and make it accessible through the system's path or try another ",
            "version!")
}

downloadPrsice <- function(ver=c("v2.3.3"),path=NULL) {
    # This downloads the latest version by default... Might not work for older
    # Linux operating systems
    ver <- ver[1]
    .checkTextArgs("PRSice version (ver)",ver,c("v2.3.3"),multiarg=FALSE)
    
    base <- "https://github.com/choishingwan/PRSice/releases/download/"
    if (ver == "v2.3.3") {
        prext <- "2.3.3/"
        ext <- "PRSice_linux.zip"
        exec <- "PRSice_linux"
    }
        
    src <- paste0(base,prext,ext)
    sver <- substring(ver,2,nchar(ver))
    defPath <- .createToolPath("prsice",sver,"PRSice",path)
    
    # Download and uncompress
    disp("Preparing to download PRSice...")
    des <- file.path(tempdir(),ext)
    download.file(src,des)
    
    # Move and rename the executable to the common style
    disp("Installing...")
    unzip(des,exdir=defPath)
    Sys.chmod(file.path(defPath,"PRSice_linux"),"0755")
    
    # Test whether all properly fetched and installed
    disp("Verifying...")
    cmd <- paste0(file.path(defPath,"PRSice_linux")," --help")
    disp("Trying the command: ",cmd)
    v <- system(cmd,ignore.stdout=TRUE,ignore.stderr=TRUE)
    if (v == 0) {
        disp("PRSice seems to be working! Updating environment...")
        .EXTERNALS[["prsice"]]$dir <- defPath
        .EXTERNALS[["prsice"]]$exec <- "plink"
        .EXTERNALS[["prsice"]]$version <- .tryNumerize(sver)
    }
    else
        disp("PRSice does NOT seem to be working! Please install manually ",
            "and make it accessible through the system's path or try another ",
            "version!")
}

.getToolPath <- function(tool=.EXTERNAL_TOOLS_LIST,version=NULL,error=TRUE) {
    if (.toolAvailable(tool[1],version,error)) {
        components <- as.list(.EXTERNALS[[tool]])
        return(file.path(components$dir,components$exec))
    }
}

# Inspired from rmarkdown package
.toolAvailable <- function(tool=.EXTERNAL_TOOLS_LIST,version=NULL,error=FALSE) {
    # Ensure we've scanned for the tool required
    tool <- tool[1]
    .findTool(tool)
    # Check availability
    version <- .tryNumerize(version)
    found <- !is.null(.EXTERNALS[[tool]]$dir) && 
        (is.null(version) || .EXTERNALS[[tool]]$version >= version)
    if (error && !found)
        stop(tool,if (!is.null(version)) c("version",version,"or higher"),
            " is required to run GWAS with ",tool," but was not found!\nIf ",
            "not already in the system PATH, please check the download",
            .capFirst(tool),"() function for possible solutions.",call.=FALSE)
    return(found)
}

.findTool <- function(tool=.EXTERNAL_TOOLS_LIST,cache=TRUE) {
    tool <- tool[1]
    .checkTextArgs("tool",tool,.EXTERNAL_TOOLS_LIST,multiarg=FALSE)
    if (!is.null(.EXTERNALS[[tool]]$dir) && cache) 
        return(invisible(as.list(.EXTERNALS[[tool]])))
    
    # Define potential sources, strategy:
    # 1. Search PATH
    # 2. Search predefined locations
    switch(tool,
        snptest = {
            execs <- c("snptest_v2.5.6","snptest_v2.5.4-beta3",
                "snptest_v2.5.2","snptest")
        },
        plink = {
            execs <- "plink"
        },
        prsice = {
            execs <- "PRSice_linux"
        }
    )
    
    wh <- Sys.which(execs)
    check <- vapply(wh,.isEmpty,logical(1))
    if (all(check)) { # Not found, try automatically defined path
        contents <- dir(dirname(.getDefaultToolPath(tool)))
        ii <- grep(tool,contents)
        if (length(ii) > 0) { # Found dir, now check executables
            toolDir <- dir(dirname(.getDefaultToolPath(tool)),
                full.names=TRUE)[ii]
            execCands <- match(execs,dir(toolDir))
            if (all(is.na(execCands))) # Shame!
                warning("Unable to locate ",tool," executable! Please ",
                    "install manually or check the download",.capFirst(tool),
                    "() function for possible solutions.")
            else { # Finally!
                .EXTERNALS[[tool]]$dir <- toolDir
                .EXTERNALS[[tool]]$exec <- 
                    execs[execCands[which(!is.na(execCands))[1]]]
                .EXTERNALS[[tool]]$version <- 
                    .tryNumerize(strsplit(toolDir,"-")[[1]][2])
            }
        }
    }
    else {
        toolExec <- wh[which(!check)]
        .EXTERNALS[[tool]]$dir <- dirname(toolExec)
        .EXTERNALS[[tool]]$exec <- basename(toolExec)
        .EXTERNALS[[tool]]$version <- NULL
    }
    return(invisible(as.list(.EXTERNALS[[tool]])))
}

.createToolPath <- function(tool,ver=NULL,name,path=NULL) {
    if (!is.null(path) && is.character(path)) {
        d <- tool
        if (!is.null(ver))
            d <- paste0(d,"-",ver)
        defPath <- file.path(path,d)
        if (!dir.exists(defPath))
            dir.create(defPath,recursive=TRUE)
        disp("\n",name," will be downloaded and installed in user specified ",
            "path: ",defPath,"\nPlease add this to the system PATH")
    }
    else {
        defPath <- .getDefaultToolPath(tool,ver)
        if (!interactive())
            .checkToolPath(defPath)
        else {
            if (!dir.exists(defPath)) {
                confirm <- menu(c("Yes","No"),
                    title=sprintf(paste0(name," is going to be installed in ",
                        "%s. Do you agree?"),defPath))
                if (confirm == 1)
                    .checkToolPath(defPath)
                else
                    stop(paste0(name," cannot be installed without user ",
                        "agreement.\nPlease agree or use the 'path' argument ",
                        "in the 'downloadSnptest' function."))
            }
        }
    }
    return(defPath)
}

.checkToolPath <- function(tp=NULL,to=.EXTERNAL_TOOLS_LIST,v=NULL) {
    to <- to[1]
    if (is.null(tp))
        tp <- .getDefaultToolPath(to,v)
    if (!dir.exists(tp)) {
        disp(to," home directory created at ",tp)
        dir.create(tp,recursive=TRUE,mode="0755",showWarnings=FALSE)
    }
    else
        disp(tp," home directory found at ",tp)
}

.getDefaultToolPath <- function(tool,version=NULL) {
    d <- tool
    if (!is.null(version))
        d <- paste0(d,"-",version)
    return(file.path(tools::R_user_dir(d)))
}

.tryNumerize <- function(x) {
    if (is.null(x))
        return(NULL)
    else if (is.numeric(x))
        return(x)
    else if (!is.na(suppressWarnings(as.numeric(x))))
        return(as.numeric(x))
    else { # Try to strip letters and non-numericals and convert
        rx <- "[A-Za-z]|[\\\\/\\;\\:\\,\\.\\~\\!\\@\\#\\%\\^\\&\\*\\(\\)\\_\\-]"
        x <- gsub(rx,"",x,perl=TRUE)
        if (!.isEmpty(x)) # Some numerics remain
            return(as.numeric(x))
        else # NULL if nothing can be done
            return(NULL)
    }
}
