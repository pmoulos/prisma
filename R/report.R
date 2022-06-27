prismaReport <- function(gwe,prismaOut,cvMetricsOut,lookupOut,path=NULL) {
    if (is.null(path)) { # prsWorkspace from prismaOut
        path <- paste0("prisma_report_",format(Sys.time(),"%Y%m%d%H%M%S"))
        if (!is.null(prismaOut$params$args$prsWorkspace))
            path <- file.path(prismaOut$params$args$prsWorkspace,path)
    }
    if (!dir.exists(path))
        dir.create(path,showWarnings=FALSE,recursive=TRUE)
    
    # Now create the respective files containing the PRS candidates (an Excel)
    # file with multiple tabs
    if (!dir.exists(file.path(path,"candidates")))
        dir.create(file.path(path,"candidates"),showWarnings=FALSE,
            recursive=TRUE)
    for (m in names(prismaOut$results)) {
        out <- file.path(path,"candidates",paste0(m,"_candidates.xlsx"))
        write.xlsx(prismaOut$results[[m]]$candidates,file=out,keepNA=TRUE)
    }
    
    if (!is.null(lookupOut)) {
       for (m in names(lookupOut)) {
           out <- file.path(path,"candidates",paste0(m,"_evidence.xlsx"))
           write.xlsx(lookupOut[[m]],file=out,keepNA=TRUE)
       }
    }
    
    REP_ENV <- .makeReportEnv(gwe,prismaOut,cvMetricsOut,lookupOut)
    
    file.copy(#from="/media/raid/software/prisma/inst/prisma_report.Rmd",
        from=system.file(package="prisma","prisma_report.Rmd"),
        to=file.path(path,"prisma_report.Rmd"),overwrite=TRUE)
    render(
        #input="/media/raid/tmp/prisma_report/prisma_report.Rmd",
        input=file.path(path,"prisma_report.Rmd"),
        output_file="index.html",
        #output_dir="/media/raid/tmp/prisma_report",
        output_dir=path,
        envir=REP_ENV
    )
    unlink(file.path(path,"prisma_report.Rmd"),force=TRUE)
}

.makeReportEnv <- function(gwe,prismaOut,cvMetricsOut,lookupOut) {
    re <- new.env(parent=globalenv())
    
    re$gwe <- gwe
    re$prismaOut <- prismaOut
    re$cvMetricsOut <- cvMetricsOut
    re$lookupOut <- lookupOut
    
    return(re)
}

prsEvalReport <- function(evalList,path=NULL) {
    if (is.null(path))  # prsWorkspace from prismaOut
        path <- paste0("prs_eval_report_",format(Sys.time(),"%Y%m%d%H%M%S"))
    if (!dir.exists(path))
        dir.create(path,showWarnings=FALSE,recursive=TRUE)
    
    if (is.null(names(evalList))) { # PRSs must be somehow described
        nsnp <- unlist(lapply(evalList,function(x) {
            return(nrow(x$prs))
        }))
        names(evalList) <- paste("PRS ",seq_along(evalList)," - ",nsnp,sep="")
    }
    
    REP_ENV <- new.env(parent=globalenv())
    REP_ENV$evalList <- evalList
    
    file.copy(#from="/media/raid/software/prisma/inst/eval_report.Rmd",
        from=system.file(package="prisma","eval_report.Rmd"),
        to=file.path(path,"eval_report.Rmd"),overwrite=TRUE)
    render(
        input=file.path(path,"eval_report.Rmd"),
        output_file="index.html",
        output_dir=path,
        envir=REP_ENV
    )
    unlink(file.path(path,"eval_report.Rmd"),force=TRUE)
}
