\name{lineDist}
\alias{lineDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
~~ function to calculate distance between line and point from coordinates ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
lineDist(p1, p2, fulc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p1}{
%%     ~~Describe \code{p1} here~~
}
  \item{p2}{
%%     ~~Describe \code{p2} here~~
}
  \item{fulc}{
%%     ~~Describe \code{fulc} here~~
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
function (p1, p2, fulc) 
{
    px1 <- p1[1]
    py1 <- p1[2]
    px2 <- p2[1]
    py2 <- p2[2]
    f.x <- fulc[1]
    f.y <- fulc[2]
    A <- py1 - py2
    B <- px2 - px1
    C <- px1 * py2 - px2 * py1
    l <- abs(A * f.x + B * f.y + C)/sqrt(A^2 + B^2)
    return(l)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
