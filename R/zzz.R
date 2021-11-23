.onLoad <- function(libname,pkgname) {
    prismaVerbosity("normal")
}

# Package environment to cache paths to external programs (SNPTEST, PRSice etc.)
.EXTERNAL_TOOLS_LIST <- c("snptest","plink","prsice","impute","gtool")
.EXTERNALS <- new.env()
for (e in .EXTERNAL_TOOLS_LIST) {
    .EXTERNALS[[e]]$dir <- NULL
    .EXTERNALS[[e]]$exec <- NULL
    .EXTERNALS[[e]]$version <- NULL
}




