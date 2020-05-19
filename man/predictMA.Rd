\name{predictMA}
\alias{predictMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
~~ function to predict MA ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
predictMA(bite, musc)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bite}{
%%     ~~Describe \code{bite} here~~
}
  \item{musc}{
%%     ~~Describe \code{musc} here~~
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
function (bite, musc) 
{
    MA <- mechAdv(bite, musc)$MA
    matMA <- matrix(nrow = 101, ncol = length(MA[1, -1]))
    for (i in 1:length(MA[1, -1])) {
        x <- MA[, 1]
        y <- MA[, 1 + i]
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
            p.ma <- yp.4
        }
        else {
            if (wi.akaike[2] > wi.akaike[1] & wi.akaike[2] > 
                wi.akaike[3]) {
                p.ma <- yp.3
            }
            else {
                if (wi.akaike[3] > wi.akaike[1] & wi.akaike[3] > 
                  wi.akaike[2]) {
                  p.ma <- yp.2
                }
            }
        }
        matMA[, i] <- p.ma
    }
    colnames(matMA) <- colnames(MA)[-1]
    results <- cbind(xp, matMA)
    colnames(results)[1] <- "Percent Toothrow"
    return(results)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
