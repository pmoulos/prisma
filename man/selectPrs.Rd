\name{selectPrs}
\alias{selectPrs}
\title{Select a PRS given validation metrics}
\usage{
    selectPrs(metrics, snpSelection, gwe, 
        method = c("maxima", "elbow"),
        crit = c("prs_r2", "prs_pvalue", "prs_aic"),
        stat = c("mean", "median", "none"),
        r2type = c("adjusted", "raw"), base = NULL)
}
\arguments{
    \item{metrics}{A list with metrics returned from 
    \code{\link{prsSelection}}.}
    
    \item{snpSelection}{A data frame with candidate SNPs, the
    output, or part of the output of 
    \code{\link{aggregatePrsMarkers}}, or one of the list 
    members of the output of \code{\link{selectPrs}}.}
    
    \item{gwe}{The GWASExperiment object used to derive the
    PRS candidates.}
    
    \item{method}{How to select the number of SNPs
    for the PRS after \code{niter} iterations. One of
    \code{"maxima"} (default) or \code{"elbow"}. See also
    Details.}
    
    \item{crit}{One of \code{"prs_r2"} (default), 
    \code{"prs_pvalue"} or \code{"prs_aic"}. Determines what
    criterion is used for PRS selection. See also Details.}
    
    \item{stat}{One of \code{"mean"} (default),
    \code{"median"} or \code{"none"}. Determines which
    calculated statistic (mean, median) is used to summarize
    the \code{prsSelectCrit} in the evaluation metrics and
    when applicable (e.g. not applicable in p-values).}
    
    \item{r2type}{One of \code{"adjusted"} (default) or 
    \code{"raw"}. The first one uses the PRS contribution to
    the evaluation regression model and the second one the 
    R^2 of the full model.}
    
    \item{base}{A vector of preferred R^2 values to serve as
    baseline R^2 for further result refinement.}
}
\value{
    A named list with names \code{"main"} and 
    \code{"others"}. The former contains the best candidate
    according the selection criteria and the latter is a
    list with the rest of the candidates, named with the
    number of SNPs in each candidate.
}
\description{
    This suggests a list of PRS candidates given a series
    of evaluation metrics and selection criteria. It is
    part of the main PRISMA pipeline.
}
\details{
    Regarding the \code{method} argument, this controls how 
    the final PRS is selected based on the given evaluation
    metrics (typically the PRS R^2). When \code{"maxima"} 
    (default), a distribution of PRS R^2 values is created 
    with the PRS R^2 values (or rarer, other chosen metric) 
    derived from \code{niter} iterations. This distribution 
    (if R^2-related) may further be adjusted according to 
    the number of SNPs comrpising in the PRS, according to 
    the \code{crit} argument. Then, from this distribution, 
    local maxima are calculated. The local maximum with the 
    highest R^2 is considered the best candidate and 
    provides the number of SNPs. The other maxima are 
    returned also as further candidates. When \code{"elbow"} 
    the elbow point method is used for determining the 
    optimal number of SNPs. The elbow method is a heuristic 
    used in determining the number of features in a data set 
    with respect to some optimization criterion. In this case 
    the criterion is the value of PRS R^2 and the method 
    consists of plotting the PRS explained variation (as 
    depicted by R^2) as a function of the number of SNPs 
    in the PRS.
    
    Regarding the \code{crit} argument, this controls the 
    criterion used to select the final PRS after the end of 
    the iterative process and according to the evaluation 
    metrics. When \code{"prs_r2"} (default), the decision 
    making plots and the \code{prsSelectMethod} operate on 
    the PRS R^2 calculated in each PRS extraction iteration. 
    When \code{"prs_pvalue"}, the -log_10(PRS p-valye) is 
    used and when \code{"prs_aic"}, the AIC (Akaike's 
    Information Criterion) is used.
}
\examples{
    # Soon available
}
\author{
    Panagiotis Moulos
}
