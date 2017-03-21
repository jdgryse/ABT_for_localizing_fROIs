##################################################
##    FUNCTIES OM B1, SE(B1) EN T TE BEREKENEN    ##
##################################################

tbeta <- function(volume,design){

  #volume <- sim.data ; array met het gesmoothe signaal
  #design <- pred    ; design wordt meegegeven
  
  time <- length(design) #lengte van de tijdsreeksen: 400
  dimens <- dim(volume) #dimensie van de signaal array (in het voorbeeld: 40 40 40 400)

  varact <- c(var(design)) # variantie van x
  meanact <- c(mean(design)) # gemiddelde van x
  b1 <- apply(volume,c(1,2,3),cov,y=design)/varact # op elke cel van de drie dimensies (dus voor elk voxel), wordt de covariantie tussen x (design) en y (signaal) geschat. Die wordt dan gedeeld door de variantie van het design (x)

}

tsebeta <- function(b1, volume, design){

  time <- length(design) #lengte van de tijdsreeksen: 400
  dimens <- dim(volume) #dimensie van de signaal array (in het voorbeeld: 40 40 40 400)

  varact <- c(var(design)) # variantie van x
  meanact <- c(mean(design)) # gemiddelde van x

  b0 <- apply(volume,c(1,2,3),mean)-b1*meanact # het gemiddelde signaal voor elke voxel - b1 * het gemiddelde voor x, het design
  predact <- array(array(b1,dim=c(dimens[1],dimens[2],dimens[3],1))%*%array(design,dim=c(1,time)),dim=dimens) # voor elke voxel wordt de b! * de designwaarde gedaan om de voorspelde waarden onder het model te berekenen
  pred <- array(rep(b0,time),dim=c(dimens[1],dimens[2],dimens[3],time)) + predact #idem, maar nu met b0 erbij gedaan
  help <- (pred-volume)^2 # berekenen van SSE (voor alle voxels samen, weer een array)
  se2 <- apply(help,c(1,2,3),sum)/(time-2) #se2 wordt voor elke voxel apart berekend door SSE te delen door n - 2
  help2 <- sum((design-mean(design))^2) # noemer voor de variantie van b1 berekenen
  sb1 <- sqrt(se2/help2) #standaardfout van b1

}

tvolume <- function(b1, sb1) {  

  tmap <- b1/sb1 # t-value voor b1

}

glmresid <- function(b1, volume, design){

  time <- length(design) #lengte van de tijdsreeksen: 400
  dimens <- dim(volume) #dimensie van de signaal array (in het voorbeeld: 40 40 40 400)

  varact <- c(var(design)) # variantie van x
  meanact <- c(mean(design)) # gemiddelde van x

  b0 <- apply(volume,c(1,2,3),mean)-b1*meanact # het gemiddelde signaal voor elke voxel - b1 * het gemiddelde voor x, het design
  predact <- array(array(b1,dim=c(dimens[1],dimens[2],dimens[3],1))%*%array(design,dim=c(1,time)),dim=dimens) # voor elke voxel wordt de b! * de designwaarde gedaan om de voorspelde waarden onder het model te berekenen
  pred <- array(rep(b0,time),dim=c(dimens[1],dimens[2],dimens[3],time)) + predact #idem, maar nu met b0 erbij gedaan
  help <- (volume-pred) # berekenen van SSE (voor alle voxels samen, weer een array)  

}

