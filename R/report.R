prismaReport <- function(gwe,prismaOut,cvMetricsOut,lookupOut,
    methods=names(prismaOut$results),path=NULL) {
    # Check the desired algorithms to report
    .checkTextArgs("GWAS methods to report (methods)",methods,
        names(prismaOut$results),multiarg=TRUE)

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
    
#~     tmp <- prismaOut
#~     tmp$results <- tmp$results[1:3]
#~     REP_ENV <- .makeReportEnv(gwe,tmp,cvMetricsOut[1:3],lookupOut[1:3])
    
    file.copy(#from="/media/raid/software/prisma/inst/prisma_report.Rmd",
        from=system.file(package="prisma","inst/prisma_report.Rmd"),
        #to="/media/raid/tmp/prisma_report/prisma_report.Rmd",overwrite=TRUE)
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
