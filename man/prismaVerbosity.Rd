\name{prismaVerbosity}
\alias{prismaVerbosity}
\title{Verbosity of output messages}
\usage{
    prismaVerbosity(level = NULL)
}
\arguments{
    \item{level}{The verbosity level, one of \code{"silent"},
    \code{"normal"} (default), \code{"full"} or 
    \code{"debug"}. See also Details}
}
\value{
    This function does not return any value. Instead, it 
    sets the R option which controls the verbosity level.
}
\description{
    This function prints (or not) helpful user messages.
}
\details{
    Regarding the \code{level} argument, when \code{"silent"}
    all output messages are suppressed, while when 
    \code{"normal"}, a reasonable informative output is
    produced. When \code{"full"}, more verbose output is
    produced for some fuction while most messages are
    produced in \code{"debug"} mode.
}
\examples{
    prismaVerbosity("silent")
    prismaVerbosity("normal")
}
\author{
    Panagiotis Moulos
}
