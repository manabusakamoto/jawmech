\name{lineDist_UniLateral}
\alias{lineDist_UniLateral}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to calculate distance between line and point from coordinates for a unilateral bite~~
function to calculate distance between line and point from coordinates for a unilateral bite
}
\description{
%%  ~~ Calculates the distance between a point and a straight line but also returns the scaling factors for the working and balancing sides of the unilateral bite ~~
}
\usage{
lineDist_UniLateral(p1, p2, p0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p1}{
%%     ~~x,y coordinates for jaw joint 1~~
}
  \item{p2}{
%%     ~~x,y coordinates for jaw joint 2~~
}
  \item{p0}{
%%     ~~x,y coordinates for bite point~~
}
}
\details{
%%  ~~ Based on three points, this function will calculate the distance between the bite point and a straight line connecting the two jaw joints for a unilateral bite. Additionally, the proportional weighting (scaling factor) will be returned for the working and balancing sides ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{lever }{Value for the calculated lever = distance between the biting point and the rotational axis, or straight line connecting the two jaw joints}
%%  \item{scale.contralateral }{Scaling factor (weight) for the contralateral side}
%%  \item{scale.ipsilateral }{Scaling factor (weight) for the contralateral side}
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
function(p1, p2, p0){
  px1 <- p1[1]
  py1 <- p1[2]
  px2 <- p2[1]
  py2 <- p2[2]
  px0 <- p0[1]
  py0 <- p0[2]
  A <- py1 - py2
  B <- px2 - px1
  C <- px1*py2 - px2*py1
  l <- abs(A*px0 + B*py0 + C)/sqrt(A^2 + B^2)
  
  d01 <- sqrt({px1 - px0}^2 + {py1 - py0}^2)
  d02 <- sqrt({px2 - px0}^2 + {py2 - py0}^2)
  
  d.p1_m <- sqrt(abs(d01^2 - l^2))
  d.p2_m <- sqrt(abs(d02^2 - l^2))
  
  d.p1_p2 <- sqrt({px2 - px1}^2 + {py2 - py1}^2)
  
  d <- c(d.p1_m, d.p2_m)
  contra.lateral <- 1 - max(d)/d.p1_p2
  ipsi.lateral <- 1 - min(d)/d.p1_p2
  
  output <- list(lever=l, scale.contralateral=contra.lateral, scale.ipsilateral=ipsi.lateral)
  return(output)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
