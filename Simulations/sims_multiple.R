#####################################################
################################################################
## FOUR REGIONS, GAUSSIAN WHITE NOISE ##########################
################################################################
#####################################################

#es: 0.005 0.015
#sigma noise: 1.5, 2
#cnr: 0.5/1.5 (0.333), 0.5/2 (0.25), 1.5/1.5 (1), 1.5/2 (0.75)


## HPC

PATH1 <- "path where your source files are located"
PATH2 <- "path where your maps will be stored"
setwd(PATH1)

args <- commandArgs(TRUE)
arrayid <- as.numeric(as.character(args))[1]
simnr <- as.numeric(as.character(args))[2]
s <- simnr + ((arrayid - 1)*100)
sigmnoise <- as.numeric(as.character(args))[3]
nscan <- as.numeric(as.character(args))[4]

set.seed(171990 + s)



## PRELIMINARY

#Libraries
library(neuRosim)
library(AnalyzeFMRI)
#library(lattice)

# Sourcing functions to calculate B1, SE(B1) AND T 
source("functions.GLM.R")

#constant parameters
es1 <- 0.01
es2 <- 0.02
base <- 100
actes1 <- es1*base
actes2 <- es2*base

fwhm <- 8
sigma <- fwhm/sqrt(8*log(2))
width <- 5


#Image characteristics
dim1 <- 64
dim2 <- 64
dim3 <- 40
imdim <- c(dim1,dim2,dim3) #image dimensions
voxdim <- c(3.5,3.5,3.51) #voxelsize
ext <- 4
TR <- 2
nregio <- 2 
coord1 <- c(22,17,11)
coord2 <- c(40,18,11)
coord <- list(coord1,coord2)

#Mask
noisemask <- read.table("rmask.txt")$x
noisemask <- array(noisemask, dim=imdim)
nobrain <- which(noisemask == 0)
nvoxels <- sum(noisemask==1)
  


#CREATING GROUND TRUTH

regions1 <- simprepSpatial(regions=1, coord=coord1, radius=ext, form="sphere", fading=0)
truthdesign <- simprepTemporal(1,1,onsets=1,effectsize = 1, durations=1,TR=1,acc=0.1)
truth1 <- simVOLfmri(design=truthdesign,image=regions1,dim=imdim,SNR=1,noise="none", template=noisemask)[,,,1]
truth1 <- ifelse(truth1 > 0, 1, 0)

regions2 <- simprepSpatial(regions=1, coord=coord2, radius=ext, form="sphere", fading=0)
truthdesign <- simprepTemporal(1,1,onsets=1,effectsize = 1, durations=1,TR=1,acc=0.1)
truth2 <- simVOLfmri(design=truthdesign,image=regions2,dim=imdim,SNR=1,noise="none", template=noisemask)[,,,1]
truth2 <- ifelse(truth2 > 0, 1, 0)

truth <- truth1 + truth2
truth <- ifelse(truth > 0, 1, 0)
notruth <- which(truth==0)
actvox <- sum(truth==1)
#levelplot(truth)



## CREATING DESIGN AND SIMULATED fMRI TIME SERIES 

#putting in the temporal signal: block design 20s ON/OFF
total.time <- nscan * TR
dur <- 20
onsets <- seq(1, total.time, 40)
# Generating a design matrix
X <- simprepTemporal(total.time,1,onsets=onsets,effectsize = 100,durations=dur,TR=TR,acc=0.1, hrf="double-gamma") 
# Generate time series for ONE active voxel
pred <- simTSfmri(design=X, base=100, SNR=1, noise="none", verbose=FALSE) 

design1 <- es1 * (pred-base) + base
design2 <- es2 * (pred-base) + base


#creating the signal in the anatomical mask
smtruth1 <- GaussSmoothArray(truth1,voxdim,ksize=width,sigma=diag(sigma,3))
smtruth2 <- GaussSmoothArray(truth2,voxdim,ksize=width,sigma=diag(sigma,3))
smtruth <- smtruth1 + smtruth2
rawsignal1 <- smtruth1 %o% design1
rawsignal2 <- smtruth2 %o% design2
rawsignal <- rawsignal1 + rawsignal2

signal <- array(NA,dim=c(imdim,nscan))

for (p in 1:nscan) { # Make sure the smoothed signal remains within the anatomical region

    slice <- rawsignal[,,,p]
    slice[notruth] <- 0
    signal[,,,p] <- slice
    rm(slice)
        
}


# creating gaussian noise in the brain mask
noisim <- array(rnorm(prod(imdim)*nscan,0,sigmnoise),dim=c(imdim,nscan))
snoisim <- GaussSmoothArray(noisim, voxdim=voxdim, ksize = width, sigma = diag(sigma,3))
	
gaussnoise <- array(NA,dim=c(imdim,nscan))
for (p in 1:nscan) { # Make sure the smoothed signal remains within the anatomical region
  
slice <- snoisim[,,,p]
slice[nobrain] <- 0
gaussnoise[,,,p] <- slice
rm(slice)
	
}
	
	
#creating the final image
data <- gaussnoise + signal

	

## CALCULATING B1, S(B1), T and P FOR EACH VOXEL
b1 <- tbeta(data,pred)
sb1 <- tsebeta(b1,data,pred)
tmap <- b1/sb1
tmap[is.na(tmap)] <- 0

	
## WRITING AWAY THE DATA
setwd(PATH2)
output.b1 <- array(NA,c(imdim))
output.sb1 <- array(NA,c(imdim))
output.tmap.null <- array(NA,c(imdim))

output.b1  <- b1
output.sb1  <- sb1
output.tmap.null  <- tmap
	
f.write.nifti(output.b1, paste("B1_Multiple_simulation_",s,"_noise_",sigmnoise,"_scans_",nscan,".nii",sep=""),size="float",nii=TRUE)
f.write.nifti(output.sb1, paste("SB1_Multiple_simulation_",s,"_noise_",sigmnoise,"_scans_",nscan,".nii",sep=""),size="float",nii=TRUE)
f.write.nifti(output.tmap.null, paste("TMAP_Multiple_null_simulation_",s,"_noise_",sigmnoise,"_scans_",nscan,".nii",sep=""),size="float",nii=TRUE)

rm("b1","sb1","tmap", "output.tmap.null","output.sb1","output.b1","data","snoisim", "noisim","signal","rawsignal")





