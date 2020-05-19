\name{mechAdv}
\alias{mechAdv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
~~ function to calculate mechanical advantages along the tooth row ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mechAdv(bite, musc, labels = NULL, average = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bite}{
%%     ~~Describe \code{bite} here~~
}
  \item{musc}{
%%     ~~Describe \code{musc} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
}
  \item{average}{
%%     ~~Describe \code{average} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (bite, musc, labels = NULL, average = NULL) 
{
    if (is.null(labels)) {
        labels <- c(seq(1:length(musc)))
    }
    if (is.null(average)) {
        average <- "mean"
    }
    bite <- as.numeric(bite)
    musc <- as.numeric(musc)
    ma <- matrix(nrow = length(bite), ncol = length(musc))
    for (i in 1:length(musc)) {
        ma[, i] <- musc[i]/bite
    }
    percent <- (bite - min(bite))/(max(bite) - min(bite)) * 100
    ma <- cbind(percent, ma)
    colnames(ma) <- c("Percent", as.character(labels))
    if (length(musc) > 1) {
        if (average == "mean") {
            ma.avg <- cbind(percent, rowMeans(ma[, -1]))
            colnames(ma.avg) <- c("Percent Toothrow", "Mechanical Advantage")
        }
        if (average == "median") {
            ma.df <- data.frame(t(ma[, -1]))
            ma.med <- unlist(lapply(ma.df, median))
            names(ma.med) <- NULL
            ma.avg <- cbind(percent, ma.med)
            colnames(ma.avg) <- c("Percent Toothrow", "Mechanical Advantage")
        }
    }
    else {
        ma.avg <- cbind(ma)
        colnames(ma.avg)[2] <- c("MA")
    }
    results <- list(ma, ma.avg)
    names(results) <- c("MA", paste(average, "MA", sep = "_"))
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
