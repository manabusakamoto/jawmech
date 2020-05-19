# calculate distance from line using three coordinates
# two for line and a third for the point of interest
# https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#An_algebraic_proof
lineDist <- function(p1, p2, p0){
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
  return(l)
}

# calculate distance from line using three coordinates
# two for line and a third for the point of interest
# and calculates scaling factors for ipsilateral and contralateral positions
lineDist_UniLateral <- function(p1, p2, p0){
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
  
  
  contra.lateral <- 1 - d.p1_m/d.p1_p2
  ipsi.lateral <- 1 - d.p2_m/d.p1_p2
  
  output <- list(lever=l, scale.contralateral=contra.lateral, scale.ipsilateral=ipsi.lateral)
  return(output)
}

# calculate distance between two points
coorDist <- function(p1, p2){
  p1.x <- p1[1]
  p1.y <- p1[2]
  p2.x <- p2[1]
  p2.y <- p2[2]
  d <- sqrt((p1.x-p2.x)^2 + (p1.y - p2.y)^2)
  return(d)
}

# calculate average muscle levers for each muscle group from X,Y coordinates
muscLever <- function(musc, fulc, method=NULL){
  if(is.null(method)){method <- "mean"}
  if(class(fulc)=="data.frame"){fulc <- unlist(fulc)}
  m.lev <- numeric(length(musc))
  for(i in 1:length(musc)){
    m.i <- musc[[i]]
    n2 <- nrow(m.i)
    
    P1 <- m.i[which(1:n2 %% 2==1), ]
    P2 <- m.i[which(1:n2 %% 2==0), ]
    
    L <- numeric(nrow(P1))
    for(j in 1:nrow(P1)){
      p1 <- unlist(P1[j,])
      p2 <- unlist(P2[j,])
      L[j] <- lineDist(p1, p2, fulc)
    }
    if(method=="mean") m.lev[i] <- mean(L)
    if(method=="median") m.lev[i] <- median(L)
  }
  return(m.lev)
}

# calculate bite levers for each toothrow position
biteLever <- function(bite, fulc){
  if(class(bite)=="list"){bite <- bite[[1]]}
  if(class(fulc)=="data.frame"){fulc <- unlist(fulc)}
  L <- numeric(nrow(bite))
  for(i in 1:nrow(bite)){
    b.i <- unlist(bite[i,])
    L[i] <- coorDist(b.i, fulc)
  }
  return(L)
}
