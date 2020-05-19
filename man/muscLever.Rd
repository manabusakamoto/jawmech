\name{muscLever}
\alias{muscLever}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
~~ function to calculate muscle levers ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
muscLever(musc, fulc, method = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{musc}{
%%     ~~Describe \code{musc} here~~
}
  \item{fulc}{
%%     ~~Describe \code{fulc} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (musc, fulc, method = NULL) 
{
    if (is.null(method)) {
        method <- "mean"
    }
    if (class(fulc) == "data.frame") {
        fulc <- unlist(fulc)
    }
    m.lev <- numeric(length(musc))
    for (i in 1:length(musc)) {
        m.i <- musc[[i]]
        n2 <- nrow(m.i)
        P1 <- m.i[which(1:n2\%\%2 == 1), ]
        P2 <- m.i[which(1:n2\%\%2 == 0), ]
        L <- numeric(nrow(P1))
        for (j in 1:nrow(P1)) {
            p1 <- unlist(P1[j, ])
            p2 <- unlist(P2[j, ])
            L[j] <- lineDist(p1, p2, fulc)
        }
        if (method == "mean") 
            m.lev[i] <- mean(L)
        if (method == "median") 
            m.lev[i] <- median(L)
    }
    return(m.lev)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
