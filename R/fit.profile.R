
#####  Function mechAdv  #####
mechAdv <- function(bite, musc, lab=NULL){
  if(is.null(lab)){lab <- paste("muscle", seq(1:length(musc)), sep="_")}
  bite <- as.numeric(bite)
  musc <- as.numeric(musc)
  ma <- matrix(nrow=length(bite),ncol=length(musc))
  for(i in 1:length(musc)){
    ma[,i]<-musc[i]/bite
  }
  percent <- (bite-min(bite))/(max(bite)-min(bite))*100
  ma <- cbind(percent,ma)
  colnames(ma) <- c("Percent",as.character(lab))
  if(length(musc)>1){
    ma <- cbind(ma,rowMeans(ma[,-1]))
    colnames(ma)[ncol(ma)] <- "mean"
    ma.df <- data.frame(t(ma[,-1]))
    ma.med <- unlist(lapply(ma.df, median))
    names(ma.med) <- NULL
    ma <- cbind(ma, ma.med)
    colnames(ma)[ncol(ma)] <- "median"
  }
  results <- ma
  return(results)
}


#####  function fit.profile  #####
fit.profile <- function(x, y){
  require(evoldiver)
  require(AICcmodavg)
  xp <- seq(0,100,1)
  lm.4 <- lm(y~x+I(x^2)+I(x^3)+I(x^4))
  lm.3 <- lm(y~x+I(x^2)+I(x^3))
  lm.2 <- lm(y~x+I(x^2))
  yp.4 <- predict(lm.4,list(x=xp))
  yp.3 <- predict(lm.3,list(x=xp))
  yp.2 <- predict(lm.2,list(x=xp))
  aic <- cbind(AIC(lm.4),AIC(lm.3),AIC(lm.2))
  aicc <- cbind(AICc(lm.4),AICc(lm.3),AICc(lm.2))
  colnames(aic) <- colnames(aicc) <- c("4th_Order","3rd_Order","2nd_Order")
  wi.akaike <- AkaikeWeights(aicc)
  if(wi.akaike[1]>wi.akaike[2] & wi.akaike[1]>wi.akaike[3]){
    #coeff <- lm.4$coeff[1:3]
    coeff <- lm.4$coeff
    pred <- yp.4
    model <- lm.4
  } else {
    if(wi.akaike[2]>wi.akaike[1] & wi.akaike[2]>wi.akaike[3]){
      #coeff<-lm.3$coeff[1:3]
      coeff<-lm.3$coeff
      pred <- yp.3
      model <- lm.3
    } else {
      if(wi.akaike[3]>wi.akaike[1] & wi.akaike[3]>wi.akaike[2]){
        #coeff <- lm.2$coeff[1:3]
        coeff <- lm.2$coeff
        pred <- yp.2
        model <- lm.2
      }
    }
  }
  n <- seq(0, {length(coeff)-1})
  names(coeff) <- paste("B", n , sep="")
  #names(coeff) <- c("B0","B1","B2")
  return(list(model=model, data=data.frame(Percent=x, MA=y), Coefficients=coeff, fitted.values=data.frame(Percent=xp, fitted.values=pred), AIC=aic, AICc=aicc, AkaikeWeights=wi.akaike))
}


#####  function profileMA  #####
profileMA <- function(bite, musc, lab=NULL){
  require(evoldiver)
  if(is.null(lab)) {lab <- paste("muscle", seq(1:length(musc)), sep="_")}
  MA <- mechAdv(bite, musc, lab)
  x <- MA[,1]
  Z <- list(MA)
  for(i in 2:ncol(MA)){
    y <- MA[,i]
    Z[[i]] <- fit.profile(x, y)
  }
  names(Z) <- c("MechanicalAdvantage", colnames(MA)[2:ncol(MA)])
  return(Z)
}


#####  function summary.profileMA  #####
summary.profileMA <- function(x, type=NULL){
  if(is.null(type)) type <- "mean"
  if(type=="mean"){
    X <- x$mean
  }
  if(type=="median"){
    X <- x$median
  }
  return(X)
}

plot.profileMA <- function(x, type=NULL, col=NULL, pch=NULL, ...){
  if(is.null(type)) type <- "all"
  if(is.null(pch)) pch <- 1
  if(type=="all"){
    n.mean <- grep("mean", names(x))
    n.med <- grep("median", names(x))
    mean <- x[n.mean][[1]]
    med <- x[n.med][[1]]
    rng <- range(x$MechanicalAdvantage[,-1])
    m <- x[-c(1,n.mean,n.med)]
    if(is.null(col)) {
      col <- c("lightseagreen", "maroon", "salmon", "palevioletred4", "turquoise4")
      col <- col[1:length(m)]
    }
    for(i in 1:length(m)){
      X <- m[[i]]
      dt <- X$data
      fit <- X$fitted.values
      if(i==1){
        plot(dt, xlab="Percent Toothrow", ylab="Mechanical advantage", col=col[i], pch=pch, ylim=rng, ...)
        lines(fit$Percent, fit$fitted.values, col=col[i])
      }else{
        points(dt, xlab="Percent Toothrow", ylab="Mechanical advantage", col=col[i], pch=pch)
        lines(fit$Percent, fit$fitted.values, col=col[i])
      }
    }
    lines(mean$fitted.values[,1], mean$fitted.values[,2], xlab="Percent Toothrow", ylab="Mechanical advantage", col="red", lwd=3, lty=2)
    lines(med$fitted.values[,1], med$fitted.values[,2], xlab="Percent Toothrow", ylab="Mechanical advantage", col="blue", lwd=3, lty=2)
    legend("topright", c(names(m),"mean","median"), col=c(col,"red", "blue"), lwd=c(rep(1, length(col)),3,3), lty=c(rep(1, length(col)),2,2))
  }
  if(type=="mean"){
    if(is.null(col)) col <- "red"
    X <- x$mean
    dt <- X$data
    fit <- X$fitted.values
    plot(dt, xlab="Percent Toothrow", ylab="Mechanical advantage", col=col, pch=pch, ylim=rng, ...)
    lines(fit$Percent, fit$fitted.values, col=col)
  }
  if(type=="median"){
    if(is.null(col)) col <- "blue"
    X <- x$median
    dt <- X$data
    fit <- X$fitted.values
    plot(dt, xlab="Percent Toothrow", ylab="Mechanical advantage", col=col, pch=pch, ylim=rng, ...)
    lines(fit$Percent, fit$fitted.values, col=col)
  }
}

