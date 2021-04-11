getAPIBase <- function() {
    base <- getOption("rpgscat_base")
    if (!is.null(base) && .isValidUrl(base))
        return(base)
    else
        return(.defaultUrlBase())
}

setAPIBase <- function(url=NULL) {
    if (is.null(url))
        url <- .defaultUrlBase()
    else {
        if (!.isValidUrl(url)) {
            warning("The provided URL string does not seem to be a valid ",
                "URL! Assuming default PGS Catalog base...",immediate.=TRUE)
            url <- .defaultUrlBase()
        }
    }
    options("rpgscat_base"=url)
}

.defaultUrlBase <- function() {
    return("https://www.pgscatalog.org/rest")
}

.isValidUrl <- function(url) {
    regx <- paste0("((([A-Za-z]{3,9}:(?:\\/\\/)?)(?:[-;:&=\\+\\$,\\w]+@)",
        "?[A-Za-z0-9.-]+|(?:www.|[-;:&=\\+\\$,\\w]+@)[A-Za-z0-9.-]+)",
        "((?:\\/[\\+~%\\/.\\w-_]*)?\\??(?:[-\\+=&;%@.\\w_]*)#?(?:[\\w]*))?)")
    return(is.character(url) && grepl(regx,url,perl=TRUE))
}
