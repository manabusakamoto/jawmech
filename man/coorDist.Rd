\name{coorDist}
\alias{coorDist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  ~~ function to calculate distances between two coordinate points ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
coorDist(p1, p2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p1}{
%%     ~~Describe \code{p1} here~~
}
  \item{p2}{
%%     ~~Describe \code{p2} here~~
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
function (p1, p2) 
{
    p1.x <- p1[1]
    p1.y <- p1[2]
    p2.x <- p2[1]
    p2.y <- p2[2]
    d <- sqrt((p1.x - p2.x)^2 + (p1.y - p2.y)^2)
    return(d)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
