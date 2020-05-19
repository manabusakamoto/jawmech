\name{plotAllProfiles}
\alias{plotAllProfiles}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
~~ function to plot multiple profiles ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotAllProfiles(x, col = NULL, pch = NULL, lcol = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{pch}{
%%     ~~Describe \code{pch} here~~
}
  \item{lcol}{
%%     ~~Describe \code{lcol} here~~
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
function (x, col = NULL, pch = NULL, lcol = NULL) 
{
    if (is.null(col)) {
        col = rep("gray50", length(files))
    }
    if (is.null(pch)) {
        pch = 1
    }
    if (is.null(lcol)) {
        lcol = rep("gray50", length(files))
    }
    KK <- x
    rangeMA <- matrix(nrow = length(KK), ncol = 2)
    for (i in 1:length(KK)) {
        X <- as.matrix(KK[[i]])
        musc.i <- X[1:3, ]
        bite.i <- X[-c(1:3), ]
        MA.i <- mechAdv(bite.i, musc.i)
        MA.i <- MA[grep("_", names(MA.i))][[1]]
        rangeMA.i <- c(max(MA.i[, 2]), min(MA.i[, 2]))
        rangeMA[i, ] <- rangeMA.i
    }
    minMA <- min(rangeMA)
    maxMA <- max(rangeMA)
    X <- as.matrix(KK[[1]])
    musc <- X[1:3, ]
    bite <- X[-c(1:3), ]
    MA <- mechAdv(bite, musc)
    MA <- MA[grep("_", names(MA))][[1]]
    predMA <- predictAvgMA(bite, musc)
    plot(MA, ylim = c(minMA, maxMA), pch = pch, col = col[1])
    lines(predMA, col = lcol[1])
    for (i in 2:length(KK)) {
        X <- as.matrix(KK[[i]])
        musc.i <- X[1:3, ]
        bite.i <- X[-c(1:3), ]
        MA.i <- mechAdv(bite.i, musc.i)
        MA.i <- MA.i[grep("_", names(MA.i))][[1]]
        predMA.i <- predictAvgMA(bite.i, musc.i)
        points(MA.i, pch = pch, col = col[i])
        lines(predMA.i, col = lcol[i])
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
