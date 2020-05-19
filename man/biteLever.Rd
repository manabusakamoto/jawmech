\name{biteLever}
\alias{biteLever}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bite Levers
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Calculates bite levers as distances between each biting point and the jaw joint from X,Y coordinates along the entire tooth row
}
\usage{
  biteLever(bite, fulc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bite}{
  data.frame object of X,Y coordinates with X as first column and Y as second column
}
  \item{fulc}{
X,Y coordinates for fulcrum (jaw joint), with X as first element and Y as second element
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
numeric vector of levers for each biting position
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
Manabu Sakamoto
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
function (bite, fulc) 
{
    if (class(bite) == "list") {
        bite <- bite[[1]]
    }
    if (class(fulc) == "data.frame") {
        fulc <- unlist(fulc)
    }
    L <- numeric(nrow(bite))
    for (i in 1:nrow(bite)) {
        b.i <- unlist(bite[i, ])
        L[i] <- coorDist(b.i, fulc)
    }
    return(L)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
