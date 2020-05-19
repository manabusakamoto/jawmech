\name{plotMA}
\alias{plotMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
~~ function to plot mechanical advantage profile ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotMA(bite, musc, col = NULL, title = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bite}{
%%     ~~Describe \code{bite} here~~
}
  \item{musc}{
%%     ~~Describe \code{musc} here~~
}
  \item{col}{
%%     ~~Describe \code{col} here~~
}
  \item{title}{
%%     ~~Describe \code{title} here~~
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
function (bite, musc, col = NULL, title = NULL) 
{
    if (is.null(col)) {
        col <- "red"
    }
    MA <- mechAdv(bite, musc)
    MA <- MA[grep("_", names(MA))][[1]]
    if (is.null(title)) {
        plot(MA, axes = F)
        box()
        axis(1)
        axis(2, las = 2)
    }
    plot(MA, main = title, axes = F)
    box()
    axis(1)
    axis(2, las = 2)
    x <- MA[, 1]
    y <- MA[, 2]
    xp <- seq(0, 100, 1)
    lm.4 <- lm(y ~ x + I(x^2) + I(x^3) + I(x^4))
    lm.3 <- lm(y ~ x + I(x^2) + I(x^3))
    lm.2 <- lm(y ~ x + I(x^2))
    yp.4 <- predict(lm.4, list(x = xp))
    yp.3 <- predict(lm.3, list(x = xp))
    yp.2 <- predict(lm.2, list(x = xp))
    aic <- cbind(AIC(lm.4), AIC(lm.3), AIC(lm.2))
    colnames(aic) <- c("4th", "3rd", "2nd")
    wi.akaike <- AkaikeWeights(aic)
    if (wi.akaike[1] > wi.akaike[2] & wi.akaike[1] > wi.akaike[3]) {
        lines(xp, yp.4, col = col)
    }
    else {
        if (wi.akaike[2] > wi.akaike[1] & wi.akaike[2] > wi.akaike[3]) {
            lines(xp, yp.3, col = col)
        }
        else {
            if (wi.akaike[3] > wi.akaike[1] & wi.akaike[3] > 
                wi.akaike[2]) {
                lines(xp, yp.2, col = col)
            }
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
