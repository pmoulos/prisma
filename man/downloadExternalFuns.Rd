\name{downloadExternalFuns}
\alias{downloadExternalTools}
\alias{downloadGtool}
\alias{downloadImpute}
\alias{downloadPlink}
\alias{downloadPrsice}
\alias{downloadSnptest}
\alias{downloadQctool}
\alias{download1000GP3}
\title{Functions to download external GWA and PRS tools}
\usage{
    downloadExternalTools(tools = c("snptest","plink","prsice",
        "impute","gtool","qctool"),
        ver = getDefaults("externalTools"), path = NULL)
    
    downloadGtool(ver = c("v0.7.5"), path = NULL)
    
    downloadImpute(ver = c("v2.3.2"), path = NULL)
    
    downloadPlink(ver = c("v1.90"), path = NULL)
    
    downloadPrsice(ver = c("v2.3.3"), path = NULL) 
    
    downloadSnptest( ver = c("v2.5.6", "v2.5.4", "v2.5.2"),
        path = NULL)
    
    downloadQctool(ver = c("v2.3.2"), path = NULL)
    
    download1000GP3(path = NULL)
}
\arguments{
    \item{tools}{The tools required by PRISMA for some of
    its operations. One or more from \code{"snptest"},
    \code{"plink"},\code{"prsice"},\code{"impute"},
    \code{"gtool"}, \code{"qctool"}.}
    
    \item{ver}{One of the supported tool versions. See the
    usage of each function.}
    
    \item{path}{Each tool's installation path or the path
    to download 1000 genomes imputation reference files.}
}
\value{
    These functions do not return a value, only update the
    environment holding paths to these tools. 
}
\description{
    These functions can be used to download 3rd party tools
    used within PRISMA. The user may also download them 
    outside R and make sure that they are findable by the
    system, for example being included in the environmental 
    \code{PATH} variable. If the tool is installed from
    within R and the automatically assigned directories,
    then it will also be findable from PRISMA. If however
    the user chooses a custom path, (e.g. using the 
    \code{path} argument), this will have to be added to the
    \code{PATH} environmental variable in order to be 
    findable from PRISMA.
    After the download and installation process is complete 
    (if any), these funcitons update the environmental variable 
    (an R environment essentially).
}
\examples{
\dontrun{
    downloadExternalTools(c("snptest","plink"))
    downloadQctool()
    download1000GP3()
}    
}
\author{
    Panagiotis Moulos
}

