##################################################
##    FUNCTIONS TO CALCULATE B1, SE(B1) AND T   ##
##################################################

tbeta <- function(volume,design){

  #volume <- sim.data ; array with smoothed signal
  #design <- pred    ; design 
  
  time <- length(design) #length of time series: 400
  dimens <- dim(volume) #dimension of signal array (in the example: 40 40 40 400)

  varact <- c(var(design)) # variance of x
  meanact <- c(mean(design)) # average of x
  # for every voxel, we estimate the covariance between the design (x) and the signal (y). Then divided by variance of design (x)
  b1 <- apply(volume,c(1,2,3),cov,y=design)/varact 

}

tsebeta <- function(b1, volume, design){

  time <- length(design) #length of time series: 400
  dimens <- dim(volume) #dimension of signal array (in the example: 40 40 40 400)

  varact <- c(var(design)) # variance of x
  meanact <- c(mean(design)) # average of x

  # mean signal for every voxel - b1 * average of x, the design
  b0 <- apply(volume,c(1,2,3),mean)-b1*meanact 
  # in every voxel: b1 * value of the design to calculate predicted value under the model
  predact <- array(array(b1,dim=c(dimens[1],dimens[2],dimens[3],1))%*%array(design,dim=c(1,time)),dim=dimens) 
  # same, but with b0 added
  pred <- array(rep(b0,time),dim=c(dimens[1],dimens[2],dimens[3],time)) + predact 
  # calculation of SSE (for all voxels together)
  help <- (pred-volume)^2 
  # se2: calculated for every voxel through SSE/(n - 2)
  se2 <- apply(help,c(1,2,3),sum)/(time-2)
  # denominator for variance of b1
  help2 <- sum((design-mean(design))^2)
  # standard error of b1
  sb1 <- sqrt(se2/help2)

}

tvolume <- function(b1, sb1) {  

  tmap <- b1/sb1 # t-value for b1

}

glmresid <- function(b1, volume, design){

  time <- length(design) #length of time series: 400
  dimens <- dim(volume) #dimension of signal array (in the example: 40 40 40 400)

  varact <- c(var(design)) # variance of x
  meanact <- c(mean(design)) # average of x

  # mean signal in every voxel - b1 * average of x, the design
  b0 <- apply(volume,c(1,2,3),mean)-b1*meanact 
  # in every voxel: b1 * value of the design to calculate predicted value under the model
  predact <- array(array(b1,dim=c(dimens[1],dimens[2],dimens[3],1))%*%array(design,dim=c(1,time)),dim=dimens) 
  # Same but with b0 added
  pred <- array(rep(b0,time),dim=c(dimens[1],dimens[2],dimens[3],time)) + predact
  # calculation of SSE (for all voxels together)
  help <- (volume-pred) 

}

